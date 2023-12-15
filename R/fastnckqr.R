#' Solve the non-crossing kernel quantile regression
#'
#' Trains the kernel quantile regression
#'
#' @param x A numerical input matrix. The dimension is \eqn{n+1} by \eqn{ntau} by \eqn{L1} by \eqn{L2}.
#' @param y Response variable. The length is \eqn{n}.
#' @param lambda1 A user-supplied \code{lambda1} sequence. The length is \eqn{L1}.
#' @param lambda2 A user-supplied \code{lambda2} sequence. The length is \eqn{L2}.
#' @param tau A user-supplied \code{tau} sequence for quantile levels. The length is \eqn{ntau}.
#' @param delta The smoothing index for \code{method='huber'}. Default is 0.125.
#' @param eps Stopping criterion.
#' @param maxit Maximum number of iterates.
#' @param gam A small number for numerical stability.
#' @param sigma Kernel bandwidth.
#' @param kernel Name of kernel function. Default is "Gaussian".
#' @param is_exact Exact or approximated solutions.
#'
#' @details
#' The function implements the majorization-minimization method to solve
#' non-crossing kernel quantile regression.
#'
#' @return
#' An object with S3 class \code{fastnckqr}
#' \item{alpha}{An \eqn{n+1} by \eqn{L} matrix of coefficients, where \eqn{n} represents the number of observations,
#' \eqn{ntau} represents the number of quantile levels, and \eqn{L} denotes the number of tuning parameters.}
#' \item{tau}{The \code{tau} sequence that was actually used.}
#' \item{lambda1}{The \code{lambda1} sequence that was actually used.}
#' \item{lambda2}{The \code{lambda2} sequence that was actually used.}
#' \item{delta}{The smoothing index.}
#' \item{npass}{The total number of iterates used to train the classifier.}
#' \item{jerr}{Warnings and errors; 0 if none.}
#' \item{info}{A list includes some settings used to fit this object: \code{eps}, \code{maxit}}.
#'
#' @keywords quantile regression
#' @useDynLib QuikQuan, .registration=TRUE
#' @export
#'
#' @examples
#' lambda2 <- 1e-4
#' tau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' lambda1 <- 10^seq(-8, 2, length.out=10)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' fit <- fastnckqr(x ,y, lambda1 = lambda1 , lambda2 = lambda2,
#'                  tau = tau)
fastnckqr <- function(x ,y, lambda1, lambda2, tau, delta=.125, eps=1e-08,
                      maxit=5e+06, gam=1e-07, sigma=NULL, kernel="rbfdot",
                     is_exact=FALSE){
  nobs <- nrow(x)
  p <- ncol(x)
  tau <- sort(tau)
  ntau <- length(tau)
  if (is.null(sigma)) sigma <- sigest(x)
  Kmat <- kernelMat(x, x, sigma = sigma, kernel = kernel)
  jerr <- c(0,0)
  if (missing(lambda1)) {
    stop("Users have to provide a lambda1 sequence.")
  } else {
    ulam1 = sort(lambda1, decreasing=FALSE)
    nlam1 = length(lambda1)
  }
  if (missing(lambda2)) {
    stop("Users have to provide a lambda2 sequence or a lambda2 value.")
  } else {
    ulam2 = sort(lambda2, decreasing=TRUE)
    nlam2 = length(lambda2)
  }
  npass <- array(0, c(nlam1, nlam2))
  alpmat <- array(0, c((nobs+1), ntau, nlam1, nlam2))
  
  fit <- dotCall64::.C64("fast_nckqr",
     SIGNATURE = c("double", "double", "integer", "double",
                   "integer", "double", "integer", "double",
                   "integer", "double", "double", "integer",
                   "double", "integer", "integer", "integer",
                   "integer", "integer", "double"),
    delta = as.double(delta), Kmat = as.double(Kmat), nobs = as.integer(nobs),
    y=as.double(y), nlam1=as.integer(nlam1), ulam1=as.double(ulam1),
    nlam2=as.integer(nlam2), ulam2=as.double(ulam2), ntau=as.integer(ntau),
    utau=as.double(tau), eps=as.double(eps), maxit=as.integer(maxit),
    gamma=as.double(gam), is_exact=as.integer(is_exact), anlam1=integer(1),
    anlam2=integer(1), npass=as.integer(npass), jerr=as.integer(jerr),
    alpmat=as.double(alpmat),
    PACKAGE = "QuikQuan")
   alpha <- array(fit$alpmat[seq((nobs+1) * ntau * nlam1 * nlam2)],
               c((nobs+1), ntau, nlam1, nlam2))
   ############################################################################
   ## wrap up output
   info <- list(eps=eps, maxit=signif(maxit))
   outlist <- list(alpha=alpha, tau=tau, lambda1=ulam1, lambda2=ulam2,
     sigma=sigma, tau=tau, delta=fit$delta, npass=fit$npass,
     jerr=fit$jerr, info=info)
   class(outlist) <- c("fastnckqr")
   outlist
}
