#' Complete measure.
#' 
#' @param X A matrix.
#' @param D A matrix.
#' @param Q A scalar.
#' @param num_obs The number of observations.
#' @param num_dim The number of dimensions.
#' @param num_comp The number of components.
#' @param index_comp The index of components.
#'
#' @return The sum of \code{x} and \code{y}.
#'
#' @export

est_complete <- function(X, D, Q, num_obs, num_dim, num_comp, index_comp) {
  if (!(is.numeric(X) && is.numeric(D) && is.numeric(Q) && 
        is.numeric(num_obs) && is.numeric(num_dim) && 
        is.numeric(num_comp) && is.numeric(index_comp)))
    stop("arguments must be numeric")
  out <- .C("est_complete",
            X=as.double(X),
            D=as.double(D),
            Q=as.double(Q),
            NOBS=as.integer(num_obs),
            NDIM=as.integer(num_dim),
            NCOMP=as.integer(num_comp),
            ICOMP=as.integer(index_comp),
	    PACKAGE = "MDMeasure")
  return(list(D = out$D, Q = out$Q))
}