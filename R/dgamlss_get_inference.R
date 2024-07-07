#' Local site: Get local Hessian, X^TX matrices, and local sample size for distributed inference tasks.
#'
#' @param local_object Local fit of GAMLSS object. Used to get formatted matrices and coefficients to obtain X^TX matrix and formatted Hessian.
#' @param penalized_parameters String vector containing names of parameters that are to be penalized. X^TX matrices from relevant parameters will be generated.
#' @param pooled_coefs Named list of final pooled coefficients. Names must be "mu", "sigma", "nu", "tau".
#'
#' @return Communication object: List containing local Hessian, local X^TX matrix for relevant penalized parameters, and local sample size.
#' @export
#'
#' @examples
#' \dontrun{
#' local_object = gamlss(y ~ x + bs_1 + bs_2 + bs_3 + bs_4, ~ x, ~ 1, data = df)
#' dgamlss_get_inference(local_object, penalized_parameters = "mu",
#' pooled_coefs = NULL) # Only gives X^TX matrix.
#' dgamlss_get_inference(local_object, penalized_parameters = "mu",
#' pooled_coefs = pooled_coefs)
#' }
dgamlss_get_inference <- function(local_object,
                                  penalized_parameters = NULL,
                                  pooled_coefs = NULL) {
  dgamlss_get_hessian <- function(local_object, pooled_coefs) {
    replace_local_coefs <- function(local_gamlss, pooled_coefs) {
      names(pooled_coefs$mu) <- names(local_gamlss$mu.coefficients)
      local_gamlss$mu.coefficients <- pooled_coefs$mu

      if (!is.null(pooled_coefs$sigma)) {
        names(pooled_coefs$sigma) <- names(local_gamlss$sigma.coefficients)
        local_gamlss$sigma.coefficients <- pooled_coefs$sigma
      }

      if (!is.null(pooled_coefs$nu)) {
        names(pooled_coefs$nu) <- names(local_gamlss$nu.coefficients)
        local_gamlss$nu.coefficients <- pooled_coefs$nu
      }

      if (!is.null(pooled_coefs$tau)) {
        names(pooled_coefs$tau) <- names(local_gamlss$tau.coefficients)
        local_gamlss$tau.coefficients <- pooled_coefs$tau
      }

      return(local_gamlss)
    }

    ## end of local ----------------------------------------------------------------
    local_object <- replace_local_coefs(local_object, pooled_coefs)

    if (!is.gamlss(local_object)) {
      stop(paste("This is not an gamlss object", "\n", ""))
    }

    coefBeta <- list()
    for (i in local_object$par)
    {
      if (length(eval(parse(text = paste(paste("local_object$", i, sep = ""), ".fix==TRUE", sep = "")))) != 0) {
        ff <- eval(parse(text = paste(paste(paste(local_object$family[1], "()$", sep = ""), i, sep = ""), ".linkfun", sep = "")))
        fixvalue <- ff(fitted(local_object, i)[1])
        names(fixvalue) <- paste("fixed", i, sep = " ")
        coefBeta <- c(coefBeta, fixvalue)
      } else {
        if (!is.null(unlist(attr(terms(formula(local_object, i), specials = .gamlss.sm.list), "specials")))) {
          warning(paste("Additive terms exists in the ", i, "formula. \n ",
                        "Standard errors for the linear terms maybe are not appropriate"))
        }
        nonNAcoef <- !is.na(coef(local_object, i))
        coefBeta <- c(coefBeta, coef(local_object, i)[nonNAcoef])
      }
    }

    betaCoef <- unlist(coefBeta)
    like.fun <- dgamlss_gen_likelihood(local_object)
    hess <- optimHess(betaCoef, like.fun)
    return(hess)
  }

  if (is.null(pooled_coefs) & is.null(penalized_parameters)) {
    warning("Both pooled_coefs and penalized_parameters are NULL. Only n = number of subjects will be returned.")
  }

  if (!is.null(pooled_coefs)) {
    hessian <- dgamlss_get_hessian(local_object, pooled_coefs)
  } else {
    hessian <- NULL
  }

  if (!is.null(penalized_parameters)) {
    xtx_list <- list()
    # rss_list <- list()
    if ("mu" %in% penalized_parameters & "mu" %in% local_object$parameters) {
      xtx_list$mu <- t(local_object$mu.x) %*% local_object$mu.x
      # rss_list$mu <- sum((family$mu.linkfun(local_object$y) -
      #                       local_object$mu.x %*% pooled_coefs$mu)^2)
    }
    if ("sigma" %in% penalized_parameters & "sigma" %in% local_object$parameters) {
      xtx_list$sigma <- t(local_object$sigma.x) %*% local_object$sigma.x
      # rss_list$sigma <- sum((family$sigma.linkfun(local_object$y) -
      #                          local_object$sigma.x %*% pooled_coefs$sigma)^2)
    }
    if ("nu" %in% penalized_parameters & "nu" %in% local_object$parameters) {
      xtx_list$nu <- t(local_object$nu.x) %*% local_object$nu.x
      # rss_list$nu <- sum((family$nu.linkfun(local_object$y) -
      #                       local_object$nu.x %*% pooled_coefs$nu)^2)
    }
    if ("tau" %in% penalized_parameters & "tau" %in% local_object$parameters) {
      xtx_list$tau <- t(local_object$tau.x) %*% local_object$tau.x
      # rss_list$tau <- sum((family$tau.linkfun(local_object$y) -
      #                        local_object$tau.x %*% pooled_coefs$tau)^2)
    }
  } else {
    xtx_list <- NULL
  }

  list(hessian = hessian,
       xtx_list = xtx_list,
       n = length(local_object$y))
}
