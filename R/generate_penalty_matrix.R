#' Central Site: Generate penalty matrix based for use with fixed or automated penalized splines
#'
#' @param local_object Output from gamlss() or gamlss_mock_fit() from central site
#' @param smooth_penalty_list Named list of penalty matrices contained in "P" attribute of dgamlss_bs() output.
#' @param smooth_index_list Named list of 0/1 vectors indicating which coefficients correspond to the penalty. Length of each smooth_index_list item must be equal to number of coefficients for that parameter. Sum of 1s for each vector must be equal to size of penalty matrix for that parameter.
#'
#' @return List of penalty matrices formatted for direct use with dgamlss_aggregate_coefs()
#' @export
#'
#' @examples
#' \dontrun{
#' generate_penalty_matrix(local_object,
#' smooth_penalty_list = list(mu = attr(age_spline, "P"),
#' sigma = attr(age_spline, "P"),
#' nu = attr(age_spline, "P"),
#' tau = attr(age_spline, "P"),
#' smooth_index_list = list(mu = c(rep(1, 10), numeric(5)),
#' sigma = c(rep(1, 10), numeric(5)),
#' nu = c(rep(1, 10), 0),
#' tau = c(rep(1, 10), 0)))
#' }
generate_penalty_matrix <- function(local_object,
                                    smooth_penalty_list,
                                    smooth_index_list) {
  penalty_matrix_list <- list()
  if ("mu" %in% local_object$parameters) {
    if (("mu" %in% names(smooth_penalty_list) & !("mu" %in% names(smooth_index_list))) |
        (!("mu" %in% names(smooth_penalty_list)) & "mu" %in% names(smooth_index_list))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("mu" %in% names(smooth_penalty_list) & "mu" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$mu
      tmp_index <- smooth_index_list$mu
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$mu <- tmp_combined
    }
  }

  if ("sigma" %in% local_object$parameters) {
    if (("sigma" %in% names(smooth_penalty_list) & !("sigma" %in% names(smooth_index_list))) |
        (!("sigma" %in% names(smooth_penalty_list)) & ("sigma" %in% names(smooth_index_list)))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("sigma" %in% names(smooth_penalty_list) & "sigma" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$sigma
      tmp_index <- smooth_index_list$sigma
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$sigma <- tmp_combined
    }
  }

  if ("nu" %in% local_object$parameters) {
    if (("nu" %in% names(smooth_penalty_list) & !("nu" %in% names(smooth_index_list))) |
        (!("nu" %in% names(smooth_penalty_list)) & "nu" %in% names(smooth_index_list))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("nu" %in% names(smooth_penalty_list) & "nu" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$nu
      tmp_index <- smooth_index_list$nu
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$nu <- tmp_combined
    }
  }

  if ("tau" %in% local_object$parameters) {
    if (("tau" %in% names(smooth_penalty_list) & !("tau" %in% names(smooth_index_list))) |
        (!("tau" %in% names(smooth_penalty_list)) & "tau" %in% names(smooth_index_list))) {
      stop("Both smooth_penalty_list and smooth_index_list must be named lists with the same penalized parameters as the names.")
    }

    if ("tau" %in% names(smooth_penalty_list) & "tau" %in% names(smooth_index_list)) {
      tmp_matrix <- smooth_penalty_list$tau
      tmp_index <- smooth_index_list$tau
      if ((sum(tmp_index) != nrow(tmp_matrix)) | !all(tmp_index == 0 | tmp_index == 1)) {
        stop("smooth_index_list must be a list of vectors of 0s and 1s indicating which coefficient indices correspond to the spline. Sum of smooth_index_list vectors must be equal to sides of corresponding smooth_penalty_list matrices.")
      }
      tmp_combined <- tmp_index %*% t(tmp_index)
      tmp_combined[tmp_combined == 1] <- tmp_matrix
      penalty_matrix_list$tau <- tmp_combined
    }
  }

  penalty_matrix_list$smooth_index_list <- smooth_index_list
  penalty_matrix_list
}
