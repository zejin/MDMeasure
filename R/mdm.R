#' Mutual Dependence Measures
#'
#' \code{mdm} computes mutual dependence measures.
#'
#' @param X A matrix.
#' @param dim_comp The dimensions of components in X \code{X}.
#' @param dist_comp Logical. If \code{TRUE}, the distances between all components
#'   in X will be returned.
#' @param type The type of mutual dependence measures, including
#'   - complete measure based on complete V-statistics
#'   - simplified complete measure based on incomplete V-statistics
#'   - asymmetric measure based on distance covariance
#'   - symmectric measure based on distance covariance
#'
#' @return \code{mdm} returns a list containing the following components:
#' \item{stat}{The value of mutual dependence measure.}
#' \item{dist}{The distances between all components.}
#'
#' @references Jin, Z. and Matteson, D. S. (2017).
#'   Generalizing Distance Covariance to Measure and Test Multivariate Mutual Dependence.
#'   arXiv preprint arXiv:1709.02532.
#'   \url{https://arxiv.org/abs/1709.02532}.
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(30), 10, 3)
#' mdm(X)

mdm <- function(X, dim_comp = NULL, dist_comp = FALSE, type = c('complete')) {
  if (!is.numeric(X)) {
    stop('X must be numeric.')
  }

  X <- as.matrix(X)
  num_obs <- nrow(X)
  num_dim <- ncol(X)

  if (is.null(dim_comp)) {
    dim_comp <- rep(1, num_dim)
  }

  if (num_dim != sum(dim_comp)) {
    stop('The dimensions of X and components do not agree.')
  }

  X <- as.vector(t(X))
  num_comp <- length(dim_comp)
  index_comp <- cumsum(c(1, dim_comp)) - 1

  if (type == 'complete') {
    out <- .C("est_complete",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")
    return(list(stat = out$Q, dist = out$D))
  }
}

