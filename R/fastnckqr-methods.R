#' Extract model coefficients from a `fastnckqr` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda1` for a given
#' 'lambda2' from a [fastnckqr()] object.
#'
#' `s1` is the new vector of `lambda1` values at which predictions are requested.
#' If `s1` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#' @param object A fitted \code{fastnckqr} object.
#' @param s2 Value of the penalty parameter `lambda2` at which
#'   coefficients are required.
#' @param s1 Value(s) of the penalty parameter `lambda1` at which
#'   coefficients are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#' @seealso [fastnckqr()] and [predict.fastnckqr()].
#'
#' @return The coefficients for the non-crossing kernel quantile regression model.
#'
#' @method coef fastnckqr
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' l2 <- 1e-4
#' ttau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' l1_list <- 10^seq(-8, 2, length.out=10)
#' fit <- fastnckqr(x,y, lambda1=l1_list, lambda2=l2, tau=ttau)
#' coef(fit, s1=l1_list[1:3], s2=1e-4)

coef.fastnckqr <- function(object, s1=NULL, s2, ...) {
  rlang::check_dots_empty()
  if (length(s2) != 1)
    stop("s2 must be assigned a single value.")
  lamlist2 <- lambda.interp(object$lambda2, s2)
  alpha <- object$alpha
  nobs <- dim(alpha)[1]
  alpha <- alpha[, ,,lamlist2$left, drop = FALSE] * lamlist2$frac +
    alpha[, ,,lamlist2$right, drop = FALSE] * (1 - lamlist2$frac)
  tau <- object$tau
  ntau <- length(tau)
  l1 <- length(object$lambda1)
  if(is.null(s1)) s1 <- object$lambda1
  ls <- length(s1)
  outlist <- array(NA, c(nobs, ntau, ls))
  for (i in seq(ntau)){
    alpha_tau <- matrix(alpha[,i,,1, drop = FALSE], ncol=l1)
    b0 <- matrix(alpha_tau[1,], nrow = 1)
    rownames(b0) <- "(Intercept)"
    alpha_tau <- rbind2(b0, alpha_tau[-1,, drop = FALSE])
    if (ls!=l1) {
      vnames <- dimnames(alpha_tau)[[1]]
      dimnames(alpha_tau) <- list(NULL, NULL)
      lambda1 <- object$lambda1
      lamlist1 <- lambda.interp(lambda1, s1)
      if (ls == 1) {
        alpha_tau = alpha_tau[, lamlist1$left, drop = FALSE] * lamlist1$frac +
          alpha_tau[, lamlist1$right, drop = FALSE] * (1 - lamlist1$frac)
      } else {
        alpha_tau = alpha_tau[, lamlist1$left, drop = FALSE] %*%
          Matrix::Diagonal(ls, lamlist1$frac) +
          alpha_tau[, lamlist1$right, drop = FALSE] %*%
          Matrix::Diagonal(ls, 1 - lamlist1$frac)
      }
    }
    outlist[,i,] <- as.matrix(alpha_tau)
  }
  outlist
}

#' Predict the fitted values for a \code{fastnckqr} object.
#'
#' @param object A fitted \code{fastnckqr} object.
#' @param x The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{fastnckqr} object.
#' @param newx A matrix of new values for \code{x} at which predictions are to be made. Note
#' that \code{newx} must be of a matrix form, predict function does not accept a vector or other
#' formats of \code{newx}.
#' @param s2 Value of the penalty parameter `lambda2` at which
#'   predictions are required.
#' @param s1 Value(s) of the penalty parameter `lambda1` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#'
#' @return
#' Returns the fitted values for the non-crossing kernel quantile regression model.
#' @keywords regression kernel
#' @useDynLib QuikQuan, .registration=TRUE
#' @method predict fastnckqr
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' l2 <- 1e-4
#' ttau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' l1_list <- 10^seq(-8, 2, length.out=10)
#' fit <- fastnckqr(x,y, lambda1=l1_list, lambda2=l2, tau=ttau)
#' predict(fit, x, tail(x), s1=l1_list[1:3], s2=1e-4)

predict.fastnckqr <- function(object, x, newx=NULL, s2, s1=NULL, ...) {
  rlang::check_dots_empty()
  if (is.null(newx)) newx <- x
  if (is.null(dim(newx))) dim(newx) <- c(length(newx), 1)
  sigma <- object$sigma
  newK <- kernelMat(newx, x, sigma=sigma)
  tau <- object$tau
  ntau <- length(tau)
  out <- coef(object, s1, s2, ...)
  if(is.null(s1)) s1 <- object$lambda1
  ls <- length(s1)
  nobs <- as.integer(NROW(newx))
  outlist <- array(NA, c(nobs, ntau, ls))
  for(i in seq(ntau)){
    fit <- as.matrix(cbind2(1, newK) %*% as.matrix(out[,i,],
      nrow=nobs, byrow=TRUE))
    outlist[,i,] <- fit
  }
  outlist
}

#' @export
fitted.fastnckqr <- function(object, x, s2, s1=NULL, ...){
  fit <- predict(object, x, s2=s2, s1=s1, ...)
  fit
}

#' Calculate the crossing numbers for a \code{fastnckqr} object.
#'
#' @param object A fitted \code{fastnckqr} object.
#' @param x The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{fastnckqr} object.
#' @param newx A matrix of new values for \code{x} at which predictions are to be made.
#' @param s2 Value of the penalty parameter `lambda2` at which
#'   predictions are required.
#' @param s1 Value(s) of the penalty parameter `lambda1` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#'
#' @return
#' Returns the crossing numbers for the non-crossing kernel quantile regression model.
#' @keywords regression kernel
#' @useDynLib QuikQuan, .registration=TRUE
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' l2 <- 1e-4
#' ttau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' l1_list <- 10^seq(-8, 2, length.out=10)
#' fit <- fastnckqr(x,y, lambda1=l1_list, lambda2=l2, tau=ttau)
#' cross.fastnckqr(fit, x, s1=l1_list[1:3], s2=1e-4)

cross.fastnckqr <- function(object, x, newx = NULL, s2, s1=NULL, ...) {
  if(is.null(s1)) s1 <- object$lambda1
  fit <- predict(object, x, newx=newx, s1=s1, s2=s2, ...)
  ls <- length(s1)
  x.row <- as.integer(NROW(x))
  count <- rep(NA, ls)
  for (i in seq(ls)){
    mat <- matrix(fit[,,i], nrow = x.row)
    count[i] <- sum(mat[, -ncol(mat), drop = FALSE] > mat[, -1, drop = FALSE])
  }
  count
}

#' Calculate the objective values for a \code{fastnckqr} object.
#'
#' @param object A fitted \code{fastnckqr} object.
#' @param x The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{fastnckqr} object.
#' @param y Response variable. The length is \eqn{n}.
#' @param s2 Value of the penalty parameter `lambda2` at which
#'   predictions are required.
#' @param s1 Value(s) of the penalty parameter `lambda1` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#'
#' @return
#' Returns the objective values for the non-crossing kernel quantile regression model.
#' @keywords regression kernel
#' @useDynLib QuikQuan, .registration=TRUE
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' l2 <- 1e-4
#' ttau <- c(0.1, 0.3, 0.5, 0.7, 0.9)
#' l1_list <- 10^seq(-8, 2, length.out=10)
#' fit <- fastnckqr(x,y, lambda1=l1_list, lambda2=l2, tau=ttau)
#' objective.fastnckqr(fit, x, y, l2)
objective.fastnckqr <- function(object, x, y, s2, s1=NULL){
  if(is.null(s1)) s1 <- object$lambda1
  out <- predict(object, x, s1=s1, s2=s2)
  ls <- length(s1)
  x.row <- as.integer(NROW(x))
  err <- rep(NA, ls)
  tau <- object$tau
  K <- kernelMat(x, x, sigma = object$sigma)
  coefs <- coef(object, s1, s2)
  for (i in seq(ls)){
    mat <- matrix(out[,,i], nrow = x.row)
    err1 <- s1[i] * sum(smooth_relu(mat[, -ncol(mat), drop = FALSE] 
      - mat[, -1, drop = FALSE]))
    err2 <- sum(colMeans(sapply(1:ncol(mat),
        function(i) check_loss(y-mat[, i], tau[i]))))
    alpha <- matrix(coefs[,,i], nrow = (x.row + 1))
    alpha <- alpha[-1, , drop = FALSE]
    b <- alpha[1, ]
    Kalpha <- K %*% alpha

    err3 <- .5 * s2 * sum(sapply(1:length(tau), function(j) {
         t(alpha[, j]) %*% Kalpha[, j]
       })) + (1e-8) * sum(b^2)
    err[i] <- err1+err2+err3
  }
  err
}
