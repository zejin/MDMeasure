#' Mutual Dependence Measures
#'
#' \code{mdm} computes mutual dependence measures.
#'
#' @param X A matrix.
#' @param dim_comp The dimensions of components in X \code{X}.
#' @param dist_comp Logical. If \code{TRUE}, the distances between all components
#'   in X will be returned.
#' @param type The type of mutual dependence measures, including
#'   - \code{asym_dcov}: asymmetric measure based on distance covariance
#'   - \code{sym_dcov}: symmectric measure based on distance covariance
#'   - \code{comp}: complete measure based on complete V-statistics
#'   - \code{comp_simp}: simplified complete measure based on incomplete V-statistics
#'   - \code{asym_comp}: asymmetric measure based on complete measure
#'   - \code{asym_comp_simp}: simplified asymmetric measure based on simplified complete measure
#'   - \code{sym_comp}: symmectric measure based on complete measure
#'   - \code{sym_comp_simp}: simplified symmectric measure based on simplified complete measure
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
#' mdm(X, type = 'asym_dcov')
#' mdm(X, type = 'sym_dcov')
#' mdm(X, type = 'comp')
#' mdm(X, type = 'comp_simp')
#' mdm(X, type = 'asym_comp')
#' mdm(X, type = 'asym_comp_simp')
#' mdm(X, type = 'sym_comp')
#' mdm(X, type = 'sym_comp_simp')

mdm <- function(X, dim_comp = NULL, dist_comp = FALSE, type = c('comp_simp')) {
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

  if (type == 'asym_dcov') {
    out <- .C("dCov_asymmetric",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'sym_dcov') {
    out <- .C("dCov_symmetric",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'comp') {
    out <- .C("MDM_complete",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'comp_simp') {
    out <- .C("MDM_complete_simple",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'asym_comp') {
    out <- .C("MDM_asymmetric",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'asym_comp_simp') {
    out <- .C("MDM_asymmetric_simple",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'sym_comp') {
    out <- .C("MDM_symmetric",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else if (type == 'sym_comp_simp') {
    out <- .C("MDM_symmetric_simple",
              X = as.double(X),
              D = as.double(numeric(num_comp * num_obs * num_obs)),
              Q = as.double(numeric(1)),
              NOBS = as.integer(num_obs),
              NDIM = as.integer(num_dim),
              NCOMP = as.integer(num_comp),
              ICOMP = as.integer(index_comp),
              PACKAGE = "MDMeasure")

  } else {
    stop('Invalid type. Read ?mdm for proper syntax.')
  }

  if (dist_comp) {
    return(list(stat = out$Q, dist = out$D))
  } else {
    return(list(stat = out$Q))
  }

}

