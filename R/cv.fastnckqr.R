#' cross-validation for selecting the tuning parameter 'lambda2' of non-crossing
#' kernel quantile regression
#'
#' Performs k-fold cross-validation for [fastnckqr()].
#' This function is largely similar [glmnet::cv.glmnet()].
#'
#' @param x A numerical input matrix. The dimension is \eqn{n} rows and \eqn{p} columns.
#' @param y Response variable.
#' @param kernel A kernel function.
#' @param lambda1 A user-supplied \code{lambda1} sequence.
#' @param lambda2 A user-supplied \code{lambda2} sequence.
#' @param tau A user-supplied \code{tau} sequence.
#' @param nfolds The number of folds in cross-validation. Default is 5.
#' @param foldid An optional vector which indexed the observations into each
#'   cross-validation fold. If supplied, \code{nfolds} is overridden.
#' @param delta The smoothing index for [fastnckqr()]. Default is .125.
#' @param ... Additional arguments passed into \code{fastnckqr}
#'
#' @details
#' The function computes the average cross-validation error and reports the standard error.
#'
#' @return
#' An object of class [cv.fastnckqr()] is returned, which is a
#'   list with the components describing the cross-validation error.
#' \item{lambda2}{The \code{lambda2} candidate values.}
#' \item{cvm}{Mean cross-validation error.}
#' \item{cvsd}{Estimates of standard error of cross-validation error.}
#' \item{cvup}{The upper curve: \code{cvm + cvsd}.}
#' \item{cvlo}{The lower curve: \code{cvm - cvsd}.}
#' \item{lambda.min}{The \code{lambda2} incurring the minimum cross-validation error.}
#' \item{lambda.1se}{The largest \code{lambda} whose error is within one standard error of the minimum.}
#' \item{cv.min}{The cross-validation error at \code{lambda.min}.}
#' \item{cv.1se}{The cross-validation error at \code{lambda.1se}.}
#' @keywords kernel quantile regression
#' @export
#'
#' @examples
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' ttau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' l1_list <- 10^seq(-8, 2, length.out=10)
#' l2_list <- c(0.01, 0.0001)
#' cvres <- cv.fastnckqr(x, y, ttau, lambda1 = l1_list, lambda2 = l2_list)


cv.fastnckqr <- function(x, y, tau, lambda1 = NULL, lambda2 = NULL,
  sigma = NULL, nfolds=5L, foldid, delta=.125, kernel = "rbfdot", ...){
  ############################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  ntau <- length(tau)
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ###Now fit the nfold models and store them
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = x.row))
  } else nfolds <- max(foldid)
  if (nfolds < 3L)
    stop("nfolds must be at least 3; nfolds = 5 recommended")
  lambda2 <- sort(lambda2, decreasing=TRUE)
  if (is.null(sigma)) sigma <- sigest(x)
  cvraw <- matrix(0, x.row, length(lambda2))
  for(t in seq(ntau)) {
    outlist <- as.list(seq(nfolds))
    for (i in seq(nfolds)) {
      which <- foldid == i
      outlist[[i]] <- fastkqr(x=x[!which, , drop=FALSE], y=y[!which],
        lambda=lambda2, tau=tau[t], sigma = sigma)
      if (outlist[[i]]$jerr != 0)
        stop(paste("Error occurs when fitting the", i, "th folder."))
    }
    cvstuff <- cvpath(outlist, x, y, tau[t], lambda2,
      foldid, x.row, delta, ... )
    cvraw <- cvraw + cvstuff$cvraw
  }
  N <- length(y) - apply(is.na(cvraw), 2, sum)
  cvm <- colMeans(cvraw, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(colMeans(scaled, na.rm = TRUE) / (N - 1))
  ## wrap up output 
  out <- list(lambda2=lambda2, cvm=cvm, cvsd=cvsd,
    cvup=cvm + cvsd, cvlo=cvm - cvsd)
  obj <- c(out, as.list(getmin(lambda2, cvm, cvsd)))
  class(obj) <- "cv.fastnckqr"
  obj
}

