#' Mutual Independence Tests
#'
#' \code{mdm_test} tests mutual independence of all components in \code{X},
#' where each component contains one variable (univariate) or more variables (multivariate).
#' All tests are implemented as permutation tests.
#'
#' @param X A matrix or data frame, where rows represent samples, and columns represent variables.
#' @param dim_comp The numbers of variables included by all components in \code{X}.
#'   If omitted, each component is assumed to contain exactly one variable.
#' @param num_perm The number of permutation samples drawn to approximate the asymptotic distributions
#'   of mutual dependence measures.
#' @param type The type of mutual dependence measures, including
#'   (1) \code{asym_dcov}: asymmetric measure based on distance covariance;
#'   (2) \code{sym_dcov}: symmectric measure based on distance covariance;
#'   (3) \code{comp}: complete measure based on complete V-statistics;
#'   (4) \code{comp_simp}: simplified complete measure based on incomplete V-statistics;
#'   (5) \code{asym_comp}: asymmetric measure based on complete measure;
#'   (6) \code{asym_comp_simp}: simplified asymmetric measure based on simplified complete measure;
#'   (7) \code{sym_comp}: symmectric measure based on complete measure;
#'   (8) \code{sym_comp_simp}: simplified symmectric measure based on simplified complete measure.
#'
#' @return \code{mdm_test} returns a list containing the following components:
#' \item{stat}{The value of mutual dependence measure.}
#' \item{pval}{The p-value of mutual independence test.}
#'
#' @references Jin, Z. and Matteson, D. S. (2017).
#'   Generalizing Distance Covariance to Measure and Test Multivariate Mutual Dependence.
#'   arXiv preprint arXiv:1709.02532.
#'   \url{https://arxiv.org/abs/1709.02532}.
#'
#' @include mdm.R
#'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(30), 10, 3)
#' mdm_test(X, type = 'comp_simp')

mdm_test <- function(X, dim_comp = NULL, num_perm = NULL, type = c('comp_simp')) {
  X <- as.matrix(X)
  num_obs <- nrow(X)
  num_dim <- ncol(X)

  if (is.null(dim_comp)) {
    dim_comp <- rep(1, num_dim)
  }

  if (num_dim != sum(dim_comp)) {
    stop("The dimensions of X and components do not agree.")
  }

  num_comp <- length(dim_comp)

  out_sample <- mdm(X, dim_comp = dim_comp, dist_comp = TRUE, type = type)
  mdm_sample <- out_sample$stat
  dist_sample <- out_sample$dist

  if (is.null(num_perm)) {
    num_perm <- 200 + 5000 %/% num_obs
  }

  count <- 0
  for (i in 1:num_perm) {
    index_perm <- matrix(NA, num_comp, num_obs)

    # keep the order of the first component
    index_perm[1, ] <- 1:num_obs

    # permute the orders of other components
    for (j in 2:num_comp) {
      index_perm[j, ] <- sample(num_obs)
    }

    index_perm <- as.vector(index_perm) - 1

    if (type == "asym_dcov") {
      out_perm <- .C("dCov_asymmetric_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "sym_dcov") {
      out_perm <- .C("dCov_symmetric_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "comp") {
      out_perm <- .C("MDM_complete_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "comp_simp") {
      out_perm <- .C("MDM_complete_simple_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "asym_comp") {
      out_perm <- .C("MDM_asymmetric_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "asym_comp_simp") {
      out_perm <- .C("MDM_asymmetric_simple_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "sym_comp") {
      out_perm <- .C("MDM_symmetric_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else if (type == "sym_comp_simp") {
      out_perm <- .C("MDM_symmetric_simple_perm",
                     D = as.double(dist_sample),
                     Q = as.double(numeric(1)),
                     NOBS = as.integer(num_obs),
                     NCOMP = as.integer(num_comp),
                     IPERM = as.integer(index_perm),
                     PACKAGE = "MDMeasure")
    } else {
      stop("Invalid type. Read ?mdm_test for proper syntax.")
    }

    mdm_perm <- out_perm$Q

    if (mdm_perm >= mdm_sample) {
      count <- count + 1
    }
  }

  return(list(stat = mdm_sample, pval = count / num_perm))
}

