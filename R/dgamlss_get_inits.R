#' Central site: Get initial values based on running gamlss on the central site data.
#'
#'
#' @param mu.formula Formula for the mu component of the GAMLSS model.
#' @param sigma.formula Formula for the sigma component of the GAMLSS model.
#' @param nu.formula Formula for the nu component of the GAMLSS model.
#' @param tau.formula Formula for the tau component of the GAMLSS model.
#' @param family GAMLSS family object specifying the distribution and link functions.
#' @param weights Optional weights for likelihood analysis.
#' @param contrasts Optional contrasts for model parameters.
#' @param local_site_data Data from central site. Used to set up initial values for some GAMLSS distributions.
#' #' @param all_inits Optional named list containing initial values for all coefficients.
#' @param penalty_matrix_list Named list of penalty matrices, output from dgamlss_bs(). Names must be "mu", "sigma", "nu", and "tau".
#' @param is_orthogonal TRUE/FALSE indicator of whether the dgamlss_bs() output was orthogonalized or not.
#' @param verbose TRUE/FALSE for verbosity.
#' @param ... Additional arguments to be passed to the underlying GAMLSS fitting function.
#'
#' @return
#' @export
#'
#' @examples
dgamlss_get_inits <- function(mu.formula,
                              sigma.formula = "~ 1",
                              nu.formula = "~ 1",
                              tau.formula = "~ 1",
                              family = NO(),
                              weights = NULL, # for weighted likelihood analysis
                              contrasts = NULL, # one type of contrasts for all parameters
                              local_site_data,
                              all_inits = NULL,
                              penalty_matrix_list = NULL,
                              is_orthogonal = FALSE,
                              basis_sizes = NULL,
                              verbose = FALSE,
                              ...) {
  if (verbose) {
    print("Running local GAMLSS to initialize coefficients")
  }

  local_gamlss <- gamlss(formula = formula(mu.formula),
                         sigma.formula = formula(sigma.formula),
                         nu.formula = formula(nu.formula),
                         tau.formula = formula(tau.formula),
                         family = family,
                         weights = weights,
                         data = local_site_data,
                         control = gamlss.control(trace = verbose),
                         ...)
  # Initialize coefficients ====================================================
  if (!is.null(all_inits)) {
    mu_coef_init <- all_inits$mu
    sigma_coef_init <- all_inits$sigma
    nu_coef_init <- all_inits$nu
    tau_coef_init <- all_inits$tau
  } else {
    if (!is_orthogonal) {
      if (!is.null(basis_sizes)) {
        repeat_init <- list(mu = basis_sizes[1], sigma = basis_sizes[2], nu = basis_sizes[3], tau = basis_sizes[4])
      } else {
        repeat_init <- lapply(penalty_matrix_list$smooth_index_list, sum) # If basis is not orthogonal, every basis element contributes to Intercept term
      }
    } else {
      repeat_init <- lapply(penalty_matrix_list$smooth_index_list, \(x) {1}) # If basis is orthogonal, initial value only pertains to Intercept
    }

    gamlss_family <- as.gamlss.family(family)
    mu_coef_init <- c(rep(gamlss_family$mu.linkfun(mean(local_gamlss$y)), repeat_init$mu),
                      rep(0, length(local_gamlss$mu.coefficients) - repeat_init$mu))
    y <- local_gamlss$y
    if ("sigma" %in% names(gamlss_family$parameters)) {
      sigma_coef_init <- c(rep(gamlss_family$sigma.linkfun(eval(gamlss_family$sigma.initial))[1], repeat_init$sigma),
                           rep(0, length(local_gamlss$sigma.coefficients) - repeat_init$sigma))
    } else {
      sigma_coef_init <- NULL
    }
    if ("nu" %in% names(gamlss_family$parameters)) {
      nu_coef_init <- c(rep(gamlss_family$nu.linkfun(eval(gamlss_family$nu.initial))[1], repeat_init$nu),
                        rep(0, length(local_gamlss$nu.coefficients) - repeat_init$nu))
    } else {
      nu_coef_init <- NULL
    }
    if ("tau" %in% names(gamlss_family$parameters)) {
      tau_coef_init <- c(rep(gamlss_family$tau.linkfun(eval(gamlss_family$tau.initial))[1], repeat_init$tau),
                         rep(0, length(local_gamlss$tau.coefficients) - repeat_init$tau))
    } else {
      tau_coef_init <- NULL
    }
  }

  return(list(mu = mu_coef_init,
              sigma = sigma_coef_init,
              nu = nu_coef_init,
              tau = tau_coef_init))
}
