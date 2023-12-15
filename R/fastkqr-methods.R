#' Extract model coefficients from a `fastkqr` object.
#'
#' Computes the coefficients at the requested value(s) for `lambda` from a
#' [fastkqr()] object.
#'
#' `s` is the new vector of `lambda` values at which predictions are requested.
#' If `s` is not in the lambda sequence used for fitting the model, the `coef`
#' function will use linear interpolation to make predictions. The new values
#' are interpolated using a fraction of coefficients from both left and right
#' `lambda` indices.
#'
#' @param object Fitted [fastkqr()] object.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'  coefficients are required. Default is the entire sequence.
#' @param ... Not used.
#' @seealso [fastkqr()] and [predict.fastkqr()].
#'
#' @return The coefficients at the requested values for `lambda`.
#'
#' @method coef fastkqr
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' lambda <- 10^(seq(1, -4, length.out=10)
#' fit <- fastkqr(x, y, lambda=lambda, tau=0.1)
#' coef(fit)

coef.fastkqr <- function(object, s = NULL, ...) {
  rlang::check_dots_empty()
  b0 <- matrix(object$alpha[1,], nrow = 1)
  rownames(b0) <- "(Intercept)"
  alpha <- rbind2(b0, object$alpha[-1,,drop=FALSE])
  if (!is.null(s)) {
    vnames <- dimnames(alpha)[[1]]
    dimnames(alpha) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    ls <- length(s)
    if (ls == 1) {
      alpha = alpha[, lamlist$left, drop = FALSE] * lamlist$frac +
        alpha[, lamlist$right, drop = FALSE] * (1 - lamlist$frac)
    } else {
      alpha = alpha[, lamlist$left, drop = FALSE] %*%
        Matrix::Diagonal(ls, lamlist$frac) +
        alpha[, lamlist$right, drop = FALSE] %*%
        Matrix::Diagonal(ls, 1 - lamlist$frac)
    }
    if (is.null(names(s))) names(s) <- paste0("s", seq_along(s))
    dimnames(alpha) <- list(vnames, names(s))
  }
  return(alpha)
}



#' Predict the fitted values for a \code{fastkqr} object.
#'
#' @param object A fitted \code{fastkqr} object.
#' @param x The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{fastkqr} object.
#' @param newx A matrix of new values for \code{x} at which predictions are to be made. Note
#' that \code{newx} must be of a matrix form, predict function does not accept a vector or other
#' formats of \code{newx}.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#'
#' @details
#' The result is \eqn{\beta_0 + K_i' \alpha} where \eqn{\beta_0} and \eqn{\alpha} are from the
#' \code{fastkqr} object and \eqn{K_i} is the ith row of the kernel matrix.
#'
#' @return
#' Returns the fitted values.
#' @keywords classification kernel
#' @useDynLib QuikQuan, .registration=TRUE
#' @method predict fastkqr
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' lambda <- 10^(seq(1, -4, length.out=30))
#' fit <- fastkqr(x, y, lambda=lambda, tau=0.1, is_exact=TRUE)
#' predict(fit, x, tail(x))

predict.fastkqr <- function(object, x, newx=NULL, s=NULL,...) {
  rlang::check_dots_empty()
  if (missing(newx)) newx <- x
  if (is.null(dim(newx))) dim(newx) <- c(length(newx), 1)
  sigma <- object$sigma
  alpha <- coef(object, s)
  newK <- kernelMat(newx, x, sigma=sigma)
  fit <- as.matrix(cbind2(1, newK) %*% alpha)
  fit
}


#' Return the objective values for a \code{fastkqr} object.
#'
#' @param object A fitted \code{fastkqr} object.
#' @param x The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{fastkqr} object.
#' @param y Response variable. The length is \eqn{n}.
#' @param s Value(s) of the penalty parameter `lambda` at which
#'   predictions are required. Default is the entire sequence used to create the
#'   model.
#' @param ... Not used.
#'
#' @return
#' Returns the objective values.
#' @useDynLib QuikQuan, .registration=TRUE
#' @export
#' @examples
#' library(MASS)
#' data(GAGurine)
#' x <- as.matrix(GAGurine$Age)
#' y <- GAGurine$GAG
#' lambda <- 10^(seq(1, -4, length.out=10))
#' fit <- fastkqr(x, y, lambda=lambda, tau=0.1)
#' objective.fastkqr(fit, x, y, lambda[1])

objective.fastkqr <- function(object, x, y, s = NULL) {
  if (is.null(s)) s <- object$lambda
  predmat <- predict(object, x, s = s)
  K <- kernelMat(x, x, sigma = object$sigma)
  err1 <- colMeans(check_loss(y - predmat, tau = object$tau))
  alpha <- object$alpha[-1, , drop = FALSE]
  b <- object$alpha[1, ]
  Kalpha <- K %*% alpha
  err2 <- sapply(1:length(s), function(i) {
    .5 * s[i] * (t(alpha[, i]) %*% Kalpha[, i]) + (1e-8) * b[i] * b[i]
  })
  err1 + err2
}
