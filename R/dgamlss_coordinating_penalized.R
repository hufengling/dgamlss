#' Simulation: End-to-end coordinating function for distributed GAMLSS model fitting with automated penalty selection.
#'
#' @param mu.formula Formula for the mu component of the GAMLSS model.
#' @param sigma.formula Formula for the sigma component of the GAMLSS model.
#' @param nu.formula Formula for the nu component of the GAMLSS model.
#' @param tau.formula Formula for the tau component of the GAMLSS model.
#' @param family GAMLSS family object specifying the distribution and link functions.
#' @param weights Optional weights for likelihood analysis.
#' @param contrasts Optional contrasts for model parameters.
#' @param site_data List of data frames containing local datasets.
#' @param all_inits Optional named list containing initial values for all coefficients.
#' @param coef_crit Convergence criterion for coefficient changes.
#' @param local_site_data Data from central site. Used to set up initial values for some GAMLSS distributions.
#' @param lambda_list Named list of fixed lambdas, if fixed penalty smooth terms are desired. Names must be "mu", "sigma", "nu", and "tau".
#' @param penalty_matrix_list Named list of penalty matrices, output from dgamlss_bs(). Names must be "mu", "sigma", "nu", and "tau".
#' @param method Option to use either GAIC or GCV for automated penalty selection.
#' @param k GAIC multiplier for automated penalty selection. AIC uses k = 2. BIC uses k = log(n). Higher k corresponds to stronger penalty.
#' @param is_orthogonal TRUE/FALSE indicator of whether the dgamlss_bs() output was orthogonalized or not.
#' @param basis_sizes Vector of number of basis functions for each spline. Must be length 4.
#' @param verbose TRUE/FALSE for verbosity.
#' @param ... Additional arguments to be passed to the underlying GAMLSS fitting function.
#'
#' Simulates dgamlss fitting between a central site and all local sites. Only summary statistics are passed from local sites. The central site then aggregates local summary statistics and sends back pooled updates.
#'
#' @import gamlss
#' @importFrom stats formula
#' @return A list containing pooled coefficients, number of communications, number of reduced updates, and global deviance.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data(abdom)
#' site1_data <- abdom[1:110, ]
#' site2_data <- abdom[111:610, ]
#' site_data <- list(site1_data, site2_data)
#' site1_gamlss <- gamlss(mu.formula = y ~ x1 + x2, sigma.formula = ~ x1,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1_data)
#' result <- dgamlss_coordinating(mu.formula = y ~ x1 + x2,
#' sigma.formula = ~ x1, nu.formula = ~ 1, tau.formula = ~ 1,
#' family = BCT, site_data, all_inits = dgamlss_get_inits(site1_gamlss))
#' }
dgamlss_coordinating_penalized <- function(mu.formula,
                                           sigma.formula = "~ 1",
                                           nu.formula = "~ 1",
                                           tau.formula = "~ 1",
                                           family = NO(),
                                           weights = NULL, # for weighted likelihood analysis
                                           contrasts = NULL, # one type of contrasts for all parameters
                                           local_site_data,
                                           site_data,
                                           all_inits = NULL,
                                           coef_crit = 0.01,
                                           lambda_list = NULL,
                                           penalty_matrix_list = NULL,
                                           method = c("GAIC", "GCV"),
                                           k = 2,
                                           is_orthogonal = FALSE,
                                           basis_sizes = NULL,
                                           verbose = FALSE,
                                           ...) {
  # Return one site ==========
  return_one_site <- function(site, update,
                              get_penalty_metric = FALSE,
                              proposed_coef_list = NULL) {
    dgamlss_RS(formula = formula(mu.formula),
               sigma.formula = formula(sigma.formula),
               nu.formula = formula(nu.formula),
               tau.formula = formula(tau.formula),
               family = family,
               weights = weights,
               data = site,
               dgamlss.update = update,
               get_penalty_metric = get_penalty_metric,
               proposed_coef_list = proposed_coef_list,
               ...)
  }

  method <- match.arg(method)

  outer_iter <- 0
  n_communications <- 0
  n_reduced <- -1

  inits <- dgamlss_get_inits(mu.formula,
                             sigma.formula,
                             nu.formula,
                             tau.formula,
                             family,
                             weights,
                             contrasts,
                             local_site_data,
                             all_inits,
                             penalty_matrix_list,
                             is_orthogonal,
                             basis_sizes,
                             verbose, ...)

  mu_coefs <- inits$mu
  sigma_coefs <- inits$sigma
  nu_coefs <- inits$nu
  tau_coefs <- inits$tau

  new_update <- dgamlss_send_coefs("mu",
                                   mu_coefs,
                                   sigma_coefs,
                                   nu_coefs,
                                   tau_coefs)

  if (verbose) {
    print("Initial coefficients below:")
    print(new_update)
  }

  old_update <- new_update
  old_outer_update <- new_update

  # Get XTX matrix =================================================
  n_communications <- n_communications + 1
  site_xtx_list <- lapply(site_data, function(site) {
    dgamlss_get_inference(gamlss_mock_fit(formula = formula(mu.formula),
                                          sigma.formula = formula(sigma.formula),
                                          nu.formula = formula(nu.formula),
                                          tau.formula = formula(tau.formula),
                                          family = BCPE(),
                                          data = site),
                          penalized_parameters = names(penalty_matrix_list[-length(penalty_matrix_list)]))
  })
  pooled_n <- site_xtx_list$n

  xtx_list <- list(xtx_mu = lapply(site_xtx_list, function(site) site$xtx_list$mu),
                   xtx_sigma = lapply(site_xtx_list, function(site) site$xtx_list$sigma),
                   xtx_nu = lapply(site_xtx_list, function(site) site$xtx_list$nu),
                   xtx_tau = lapply(site_xtx_list, function(site) site$xtx_list$tau))
  n <- sum(sapply(site_xtx_list, function(site) site$n))

  # Outer iteration ===========================================================
  while(TRUE) {
    outer_iter <- outer_iter + 1
    if (verbose) {
      print(paste("Outer", outer_iter))
    }
    # Update mu ================================================================================
    n_reduced <- n_reduced + 1
    is_updated <- FALSE
    new_update$to_update <- "mu"
    while(!dgamlss_check_convergence(new_update,
                                     old_update,
                                     coef_crit = coef_crit) | !is_updated) {
      if (verbose) {
        print("Updating mu")
      }

      n_communications <- n_communications + 1
      site_info <- lapply(site_data, return_one_site, update = new_update)
      if ("mu" %in% names(penalty_matrix_list)) {
        proposed_mu <- dgamlss_aggregate_coef(site_info,
                                              penalty_matrix = penalty_matrix_list$mu)
        n_communications <- n_communications + 1
        site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                get_penalty_metric = TRUE, proposed_coef_list = proposed_mu$proposed_coef_list)
        new_mu <- dgamlss_select_lambda(proposed_mu$lambda_vec, proposed_mu$proposed_coef_list,
                                        site_ssr_list, xtx_list$xtx_mu, penalty_matrix_list$mu,
                                        method = method, pooled_n = n,
                                        k = k)
        mu_coefs <- proposed_mu$proposed_coef_list[[new_mu$index]]
      } else {
        new_mu <- dgamlss_aggregate_coef(site_info)
        mu_coefs <- new_mu$coef
      }

      if (verbose) {
        print(mu_coefs)
        print(new_mu$deviance)
      }

      old_update <- new_update
      new_update <- dgamlss_send_coefs("mu",
                                       mu_coefs,
                                       sigma_coefs,
                                       nu_coefs,
                                       tau_coefs)
      is_updated <- TRUE
    }
    global_deviance <- new_mu$deviance

    if (dgamlss_check_convergence(new_update,
                                  old_outer_update,
                                  coef_crit = coef_crit) & (outer_iter > 1)) {
      break
    }

    # Update sigma ============================================================================================
    if (!is.null(inits$sigma)) {
      n_reduced <- n_reduced + 1
      is_updated <- FALSE
      new_update$to_update <- "sigma"
      while(!dgamlss_check_convergence(new_update,
                                       old_update,
                                       coef_crit = coef_crit) | !is_updated) {
        if (verbose) {
          print("Updating sigma")
        }

        n_communications <- n_communications + 1
        site_info <- lapply(site_data, return_one_site, update = new_update)
        if ("sigma" %in% names(penalty_matrix_list)) {
          proposed_sigma <- dgamlss_aggregate_coef(site_info,
                                                   penalty_matrix = penalty_matrix_list$sigma)
          n_communications <- n_communications + 1
          site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                  get_penalty_metric = TRUE, proposed_coef_list = proposed_sigma$proposed_coef_list)
          new_sigma <- dgamlss_select_lambda(proposed_sigma$lambda_vec, proposed_sigma$proposed_coef_list,
                                             site_ssr_list, xtx_list$xtx_sigma, penalty_matrix_list$sigma,
                                             method = method, pooled_n = n,
                                             k = k)
          sigma_coefs <- proposed_sigma$proposed_coef_list[[new_sigma$index]]
        } else {
          new_sigma <- dgamlss_aggregate_coef(site_info)
          sigma_coefs <- new_sigma$coef
        }
        if (verbose) {
          print(sigma_coefs)
          print(new_sigma$deviance)
        }

        old_update <- new_update
        new_update <- dgamlss_send_coefs("sigma",
                                         mu_coefs,
                                         sigma_coefs,
                                         nu_coefs,
                                         tau_coefs)
        is_updated <- TRUE
      }
      global_deviance <- new_sigma$deviance
    }

    # Update nu ==============================================================================================================
    if (!is.null(inits$nu)) {
      n_reduced <- n_reduced + 1
      is_updated <- FALSE
      new_update$to_update <- "nu"
      while(!dgamlss_check_convergence(new_update,
                                       old_update,
                                       coef_crit = coef_crit) | !is_updated) {
        if (verbose) {
          print("Updating nu")
        }

        n_communications <- n_communications + 1
        site_info <- lapply(site_data, return_one_site, update = new_update)
        if ("nu" %in% names(penalty_matrix_list)) {
          proposed_nu <- dgamlss_aggregate_coef(site_info,
                                                penalty_matrix = penalty_matrix_list$nu)
          n_communications <- n_communications + 1
          site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                  get_penalty_metric = TRUE, proposed_coef_list = proposed_nu$proposed_coef_list)
          new_nu <- dgamlss_select_lambda(proposed_nu$lambda_vec, proposed_nu$proposed_coef_list,
                                          site_ssr_list, xtx_list$xtx_nu, penalty_matrix_list$nu,
                                          method = method, pooled_n = n,
                                          k = k)
          nu_coefs <- proposed_nu$proposed_coef_list[[new_nu$index]]
        } else {
          new_nu <- dgamlss_aggregate_coef(site_info)
          nu_coefs <- new_nu$coef
        }
        if (verbose) {
          print(nu_coefs)
          print(new_nu$deviance)
        }

        old_update <- new_update
        new_update <- dgamlss_send_coefs("nu",
                                         mu_coefs,
                                         sigma_coefs,
                                         nu_coefs,
                                         tau_coefs)
        is_updated <- TRUE
      }
      global_deviance <- new_nu$deviance
    }

    # Update tau ====================================================================================================
    if (!is.null(inits$tau)) {
      n_reduced <- n_reduced + 1
      is_updated <- FALSE
      new_update$to_update <- "tau"
      while(!dgamlss_check_convergence(new_update,
                                       old_update,
                                       coef_crit = coef_crit) | !is_updated) {
        if (verbose) {
          print("Updating tau")
        }

        n_communications <- n_communications + 1
        site_info <- lapply(site_data, return_one_site, update = new_update)
        if ("tau" %in% names(penalty_matrix_list)) {
          proposed_tau <- dgamlss_aggregate_coef(site_info,
                                                 penalty_matrix = penalty_matrix_list$tau)
          n_communications <- n_communications + 1
          site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                  get_penalty_metric = TRUE, proposed_coef_list = proposed_tau$proposed_coef_list)
          new_tau <- dgamlss_select_lambda(proposed_tau$lambda_vec, proposed_tau$proposed_coef_list,
                                           site_ssr_list, xtx_list$xtx_tau, penalty_matrix_list$tau,
                                           method = method, pooled_n = n,
                                           k = k)
          tau_coefs <- proposed_tau$proposed_coef_list[[new_tau$index]]
        } else {
          new_tau <- dgamlss_aggregate_coef(site_info)
          tau_coefs <- new_tau$coef
        }
        if (verbose) {
          print(tau_coefs)
          print(new_tau$deviance)
        }

        old_update <- new_update
        new_update <- dgamlss_send_coefs("tau",
                                         mu_coefs,
                                         sigma_coefs,
                                         nu_coefs,
                                         tau_coefs)
        is_updated <- TRUE
      }
      global_deviance <- new_tau$deviance
    }

    old_outer_update <- new_update
  }

  return(list(mu_coefs = new_update$mu.coef,
              sigma_coefs = new_update$sigma.coef,
              nu_coefs = new_update$nu.coef,
              tau_coefs = new_update$tau.coef,
              n_communications = n_communications,
              n_reduced = n_reduced,
              global_deviance = global_deviance,
              edf_vec = c(new_mu$edf_opt, new_sigma$edf_opt,
                          new_nu$edf_opt, new_tau$edf_opt)))
}
