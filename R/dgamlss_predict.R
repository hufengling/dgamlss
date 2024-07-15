#' dGAMLSS predictions
#'
#' @param dgamlss_object Output from dgamlss_create_summary()
#' @param local_object Local GAMLSS object or gamlss_mock_fit() output containing data matrices.
#'
#' @return Predictions from dgamlss coefficients using local_object data
#' @export
#'
#' @examples
#' \dontrun{
#' dgamlss_predict(dgamlss_output, local_object)
#' }
dgamlss_predict <- function(dgamlss_object, local_object) {
  gamlss_family <- as.gamlss.family(dgamlss_object$family[1])
  y <- local_object$y
  new_mu.x <- local_object$mu.x
  new_sigma.x <- local_object$sigma.x
  new_nu.x <- local_object$nu.x
  new_tau.x <- local_object$tau.x

  dgamlss_prediction <- list(
    mu.fv = if ("mu" %in% dgamlss_object$parameters) {
      gamlss_family$mu.linkinv(new_mu.x %*% dgamlss_object$mu.coefficients)},
    sigma.fv = if ("sigma" %in% dgamlss_object$parameters) {
      gamlss_family$sigma.linkinv(new_sigma.x %*% dgamlss_object$sigma.coefficients)},
    nu.fv = if ("nu" %in% dgamlss_object$parameters) {
      gamlss_family$nu.linkinv(new_nu.x %*% dgamlss_object$nu.coefficients)},
    tau.fv = if ("tau" %in% dgamlss_object$parameters) {
      gamlss_family$tau.linkinv(new_tau.x %*% dgamlss_object$tau.coefficients)})

  pFamily <- match.fun(paste0("p", gamlss_family$family[1]))
  qFamily <- match.fun(paste0("q", gamlss_family$family[1]))

  if (length(dgamlss_object$parameters) == 1) {
    dgamlss_prediction$median <- qFamily(0.5, dgamlss_prediction$mu.fv)
    dgamlss_prediction$quantile <- pFamily(y, dgamlss_prediction$mu.fv)
  }
  if (length(dgamlss_object$parameters) == 2) {
    dgamlss_prediction$median <- qFamily(0.5, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv)
    dgamlss_prediction$quantile <- pFamily(y, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv)
  }
  if (length(dgamlss_object$parameters) == 3) {
    dgamlss_prediction$median <- qFamily(0.5, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                                         dgamlss_prediction$nu.fv)
    dgamlss_prediction$quantile <- pFamily(y, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                                           dgamlss_prediction$nu.fv)
  }
  if (length(dgamlss_object$parameters) == 4) {
    dgamlss_prediction$median <- qFamily(0.5, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                                         dgamlss_prediction$nu.fv, dgamlss_prediction$tau.fv)
    dgamlss_prediction$quantile <- pFamily(y, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                                           dgamlss_prediction$nu.fv, dgamlss_prediction$tau.fv)
  }

  return(dgamlss_prediction)
}

#' dGAMLSS predictions with new data
#'
#' @param dgamlss_object Output from dgamlss_create_summary()
#' @param local_object Local GAMLSS object or gamlss_mock_fit() output containing data matrices.
#' @param continuous_var Continuous variable that should be plotted on the x axis. The continuous_var (such as age) should correspond one-to-one with rows of new_mu.x, new_sigma.x, new_nu.x, new_tau.x.
#' @param new_mu.x Matrix or data frame containing columns corresponding to dgamlss coefficients for mu parameter. Categorical variables should all be "fixed" at one level to ensure smooth centile fitting (ie one set of centiles should be fit for males and another set for females, even if the same underlying model is used to fit both sets of centiles.) Can build matrices using dgamlss_bs() with appropriate knots.
#' @param new_sigma.x Matrix or data frame containing columns corresponding to dgamlss coefficients for sigma parameter. See above.
#' @param new_nu.x Matrix or data frame containing columns corresponding to dgamlss coefficients for nu parameter. See above.
#' @param new_tau.x Matrix or data frame containing columns corresponding to dgamlss coefficients for tau parameter. See above.
#' @param centiles Desired percentiles to plot. Default (0.4, 2, 10, 25, 50, 75, 90, 98, 99.6)
#'
#' @import ggplot2
#' @return Predicted centiles from dgamlss coefficients using new data with fixed categorical variables.
#' @export
#'
#' @examples
#' \dontrun{
#' dgamlss_centiles(dgamlss_output, age, new_mu.x, new_sigma.x, new_nu.x, new_tau.x)
#' }
dgamlss_centiles <- function(dgamlss_object,
                             continuous_var,
                             local_object = NULL,
                             new_mu.x = NULL,
                             new_sigma.x = NULL,
                             new_nu.x = NULL,
                             new_tau.x = NULL,
                             centiles = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6)) {
  gamlss_family <- as.gamlss.family(dgamlss_object$family[1])

  if (length(centiles) == 0) {
    stop("Must specify at least one centile")
  }

  if (!is.null(local_object)) {
    new_mu.x <- local_object$mu.x
    new_sigma.x <- local_object$sigma.x
    new_nu.x <- local_object$nu.x
    new_tau.x <- local_object$tau.x
  }

  if (length(unique(length(continuous_var),
                    nrow(new_mu.x), nrow(new_sigma.x),
                    nrow(new_nu.x), nrow(new_tau.x))) != 1) {
    stop("Length of x_axis_variable must be equivalent to number of rows in new_mu.x, new_sigma.x, new_nu.x, and new_tau.x.")
  }

  centile_order <- order(continuous_var, decreasing = FALSE)

  dgamlss_prediction <- list(
    mu.fv = if ("mu" %in% dgamlss_object$parameters) {
      gamlss_family$mu.linkinv(new_mu.x %*% dgamlss_object$mu.coefficients)},
    sigma.fv = if ("sigma" %in% dgamlss_object$parameters) {
      gamlss_family$sigma.linkinv(new_sigma.x %*% dgamlss_object$sigma.coefficients)},
    nu.fv = if ("nu" %in% dgamlss_object$parameters) {
      gamlss_family$nu.linkinv(new_nu.x %*% dgamlss_object$nu.coefficients)},
    tau.fv = if ("tau" %in% dgamlss_object$parameters) {
      gamlss_family$tau.linkinv(new_tau.x %*% dgamlss_object$tau.coefficients)})

  qFamily <- match.fun(paste0("q", gamlss_family$family[1]))

  dgamlss_prediction$centiles <- vector("list", length(centiles))
  for (i in 1:length(centiles)) {
    if (length(dgamlss_object$parameter) == 1) {
      dgamlss_prediction$centiles[[i]] <- cbind(
        x = continuous_var[centile_order],
        pred = qFamily(centiles[i] / 100, dgamlss_prediction$mu.fv)[centile_order],
        centile = centiles[i])
    }
    if (length(dgamlss_object$parameter) == 2) {
      dgamlss_prediction$centiles[[i]] <- cbind(
        x = continuous_var[centile_order],
        pred = qFamily(centiles[i] / 100, dgamlss_prediction$mu.fv,
                       dgamlss_prediction$sigma.fv)[centile_order],
        centile = centiles[i])
    }
    if (length(dgamlss_object$parameter) == 3) {
      dgamlss_prediction$centiles[[i]] <- cbind(
        x = continuous_var[centile_order],
        pred = qFamily(centiles[i] / 100, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                       dgamlss_prediction$nu.fv)[centile_order],
        centile = centiles[i])
    }
    if (length(dgamlss_object$parameter) == 4) {
      dgamlss_prediction$centiles[[i]] <- cbind(
        x = continuous_var[centile_order],
        pred = qFamily(centiles[i] / 100, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                       dgamlss_prediction$nu.fv, dgamlss_prediction$tau.fv)[centile_order],
        centile = centiles[i])
    }
  }

  centiles_df <- data.frame(do.call(rbind, dgamlss_prediction$centiles))
  centiles_df$centile <- as.factor(centiles_df$centile)

  centile_plot <- ggplot(centiles_df, aes(x, pred, color = fct_rev(centile))) +
    geom_line()
  return(centile_plot)
}
