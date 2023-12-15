#' cross-validation for selecting the tuning parameter of kernel quantile regression
#'
#' Performs k-fold cross-validation for [fastkqr()].
#' This function is largely similar [glmnet::cv.glmnet()].
#'
#' @param x A numerical input matrix. The dimension is \eqn{n} rows and \eqn{p} columns.
#' @param y Response variable.
#' @param kern A kernel function.
#' @param lambda A user-supplied \code{lambda} sequence.
#' @param nfolds The number of folds in cross-validation. Default is 5.
#' @param foldid An optional vector which indexed the observations into each
#'   cross-validation fold. If supplied, \code{nfolds} is overridden.
#' @param delta The smoothing index for [fastkqr()]. Default is .125.
#' @param ... Additional arguments passed into \code{fastkqr}
#'
#' @details
#' The function computes the average cross-validation error and reports the standard error.
#'
#' @return
#' An object of class [cv.fastkqr()] is returned, which is a
#'   list with the components describing the cross-validation error.
#' \item{lambda}{The \code{lambda} candidate values.}
#' \item{cvm}{Mean cross-validation error.}
#' \item{cvsd}{Estimates of standard error of cross-validation error.}
#' \item{cvup}{The upper curve: \code{cvm + cvsd}.}
#' \item{cvlo}{The lower curve: \code{cvm - cvsd}.}
#' \item{lambda.min}{The \code{lambda} incurring the minimum cross-validation error.}
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
#' lambda <- 10^(seq(1, -4, length.out=10)
#' cv.fit <- cv.fastkqr(x, y, lambda=lambda, tau=0.1)

cv.fastkqr <- function(x, y, tau, lambda = NULL, sigma = NULL, nfolds=5L,
                       foldid, delta=.125, kernel = "rbfdot", ...){
  ############################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ###Now fit the nfold models and store them
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = x.row))
  } else nfolds <- max(foldid)
  if (nfolds < 3L)
    stop("nfolds must be at least 3; nfolds = 5 recommended")
  lambda <- sort(lambda, decreasing=TRUE)
  if (is.null(sigma)) sigma <- sigest(x)
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- fastkqr(x=x[!which, , drop=FALSE], y=y[!which],
      lambda=lambda, tau=tau, eps=1e-05)
    if (outlist[[i]]$jerr != 0)
      stop(paste("Error occurs when fitting the", i, "th folder."))
  }
  cvstuff <- cvpath(outlist, x, y, tau, lambda, foldid, x.row, delta, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  ## wrap up output
  out <- list(lambda=lambda, cvm=cvm, cvsd=cvsd,
    cvup=cvm + cvsd, cvlo=cvm - cvsd)
  obj <- c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) <- "cv.fastkqr"
  obj
}

cvpath <- function(outlist, x, y, tau, lambda, foldid, x.row,
                       delta, kernel = "rbfdot", ...){
  nfolds <- max(foldid)
  predmat <- matrix(NA, x.row, length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[!whichfold, , drop = FALSE],
      x[whichfold, , drop = FALSE])
    nlami <- length(fitobj$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- check_loss(y-predmat, tau)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- colMeans(cvraw, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(colMeans(scaled, na.rm = TRUE) / (N - 1))
  out <- list(cvm=cvm, cvsd=cvsd, cvraw=cvraw)
  out
}
