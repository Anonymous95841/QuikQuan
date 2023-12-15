#' Solve the kernel quantile regression. The solution path is computed
#' at a grid of values of tuning parameter \code{lambda}.
#'
#' @param x A numerical input matrix. The dimension is \eqn{n} rows and \eqn{p} columns.
#' @param y Response variable. The length is \eqn{n}.
#' @param lambda A user-supplied \code{lambda} sequence.
#' @param tau A user-supplied \code{tau} value for a quantile level.
#' @param delta The smoothing index for \code{method='huber'}. Default is 0.125.
#' @param eps Stopping criterion.
#' @param maxit Maximum number of iterates.
#' @param gam A small number for numerical stability.
#' @param sigma Kernel bandwidth.
#' @param kernel A kernel function. Default is "Gaussian".
#' @param is_exact Exact or approximated solutions. Default is \code{FALSE}.
#'
#' @details
#' The function implements an accelerated proximal gradient descent to solve
#' kernel quantile regression.
#'
#' @return
#' An object with S3 class \code{fastkqr}
#' \item{alpha}{An \eqn{n+1} by \eqn{L} matrix of coefficients, where \eqn{n} is the number of observations
#' and \eqn{L} is the number of tuning parameters. The first row of \code{alpha} contains the intercepts.}
#' \item{lambda}{The \code{lambda} sequence that was actually used.}
#' \item{delta}{The smoothing index.}
#' \item{npass}{The total number of iterates used to train the classifier.}
#' \item{jerr}{Warnings and errors; 0 if none.}
#' \item{info}{A list includes some settings used to fit this object: \code{eps}, \code{maxit}}.
#'
#' @keywords quantile regression
#' @useDynLib QuikQuan, .registration=TRUE
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' lambda <- 10^(seq(1, -4, length.out=30))
#' fit <- fastkqr(x, y, lambda=lambda, tau=0.1, is_exact=TRUE)

fastkqr <- function(x ,y, lambda, tau, delta=.125, eps=1e-08, maxit=1e+06,
                       gam=1e-07, sigma = NULL, kernel="rbfdot",
                       is_exact=FALSE){
    x <- as.matrix(x)
    nrows <- as.integer(NROW(x))
    np <- as.integer(NCOL(x))
    if (length(y) != nrows)
      stop("x and y have different number of observations.")
    if (missing(lambda)) {
      stop("Users have to provide a lambda sequence.")
    } else {
      ulam = sort(lambda, decreasing=TRUE)
      nlam = length(lambda)
    }
    if (is.null(sigma)) sigma <- sigest(x)
    Kmat <- kernelMat(x,x, sigma=sigma, kernel = kernel)
    npass <- rep(0, nlam)
    alpmat <- array(0, c((nrows + 1), nlam))

    fit <- dotCall64::.C64("fast_kqr",
       SIGNATURE = c("double", "double", "integer", "double", "integer",
                     "double", "double", "double", "integer", "double",
                     "integer", "integer", "integer", "integer", "double"),
      delta = as.double(delta), Kmat = as.double(Kmat), nobs = nrows,
      y=as.double(y), nlam=as.integer(nlam), ulam=as.double(ulam),
      tau=as.double(tau), eps=as.double(eps), maxit=as.integer(maxit),
      gamma=as.double(gam), is_exact=as.integer(is_exact), anlam=integer(1),
      npass=as.integer(npass), jerr=integer(1), alpmat=as.double(alpmat),
      PACKAGE = "QuikQuan")
    alpha <- matrix(fit$alpmat[seq((nrows + 1) * nlam)], (nrows + 1), nlam)
    ############################################################################
    ## wrap up output
    info <- list(eps=eps, maxit=signif(maxit))
    outlist <- list(alpha=alpha, lambda=ulam, delta=fit$delta, sigma=sigma,
      tau=tau, npass=fit$npass, jerr=fit$jerr, info=info)
    class(outlist) <- c("fastkqr")
    outlist
}
