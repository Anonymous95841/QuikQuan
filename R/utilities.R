err <- function(n, maxit) {
# This function is adapted from glmnet package.
  if (n == 0) msg <- ""
  if (n < 0) {
    msg <- paste0("convergence for ", -n,
      "th lambda value not reached after maxit=", maxit,
      " iterations; solutions for larger lambdas returned")
    n <- -1
    msg <- paste("From kerneltool fortran code:", msg)
  }
  list(n=n, msg=msg)
}

error.bars <- function(x, upper, lower, width=0.02, ...) {
# This function is adapted from glmnet package.
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

getmin <- function(lambda, cvm, cvsd) {
# This function is adapted from glmnet package.
  cvmin <- min(cvm, na.rm=TRUE)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin], na.rm=TRUE)
  cvmin2 <- min(cvm[!is.na(cvsd)])
  lambda.min2 <- max(lambda[cvm[!is.na(cvsd)] <= cvmin2], na.rm=TRUE)
  idmin <- match(lambda.min2, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm[!is.na(cvsd)] <= semin
  lambda.1se <- max(lambda[idmin])
  id1se <- match(lambda.1se, lambda)
  cv.1se <- cvm[id1se]
  list(lambda.min=lambda.min, lambda.1se=lambda.1se,
	  cvm.min=cvmin, cvm.1se=cv.1se)
}

lambda.interp <- function(lambda, s){
# lambda is the index sequence that is produced by the model
# s is the new vector at which evaluations are required.
# the value is a vector of left and right indices, and a vector of fractions.
# the new values are interpolated bewteen the two using the fraction
# Note: lambda decreases. you take:
# sfrac*left+(1-sfrac*right)
  if (length(lambda) == 1) {
  # degenerate case of only one lambda
    nums  <- length(s)
    left  <- rep(1,nums)
    right <- left
    sfrac <- rep(1,nums)
  } else{
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1]-s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac-lambda[right])/(lambda[left] - lambda[right])
    sfrac[left==right]=1
    }
  list(left=left,right=right,frac=sfrac)
}


lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}


smooth_relu <- function(u, del=3.051758e-05){
  ifelse(u >= del, u,
         ifelse(u <= -del, 0,
                (1/(4*del)) * u^2 + (1/2) * u + (1/4) * del))
}



check_loss <- function(r, tau) {
  0.5 * (abs(r) + (2.0 * tau - 1.0) * r)
}

# kernel matrix
kernelMat <- function(x1, x2=NULL, sigma=NULL, kernel="rbfdot") {
  n1 <- nrow(x1)
  p1 <- ncol(x1)
  if (is.null(x2)) x2 <- x1
  if (is.null(sigma)) sigma <- sigest(x)
  n2 <- nrow(x2)
  p2 <- ncol(x2)
  Kmat <- array(0, c(n1, n2))
  res = dotCall64::.C64("rbfdot",
    SIGNATURE = c("double", "double", "integer", "integer",
                  "integer", "integer", "double", "double"),
    X1=as.double(x1), X2=as.double(x2), nobs1 = as.integer(n1),
    nobs2=as.integer(n2), p1=as.integer(p1), p2= as.integer(p2),
    sigma=as.double(sigma), Kmat=as.double(Kmat),
    PACKAGE = "QuikQuan")

  Kmat = matrix(res$Kmat, n1, n2)
  return(Kmat)
}

sigest = function(x, frac=0.5) {
  m = dim(x)[1]
  n = floor(frac * m)
  index = sample(1:m, n, replace = TRUE)
  index2 = sample(1:m, n, replace = TRUE)
  temp = x[index, , drop = FALSE] - x[index2, , drop = FALSE]
  dist = rowSums(temp^2)
  srange = 1/quantile(dist[dist != 0], probs = c(0.9, 0.5, 0.1))
  mean(srange[c(1, 3)])
}

objfun <- function(ab,tau,y,K,lambda){
  n = length(ab)
  alpha = ab[2:n]
  beta=ab[1]
  mean(check_loss(y-(K%*%alpha+beta),tau =tau)) +
   (lambda/2)*t(alpha)%*%K%*%alpha+(1e-8)*beta*beta
}



objfun_smoothrelu <- function(alp, ttau, y, K, lam1, lam2){
  n <- length(y)
  ntau <- length(ttau)
  alp <- matrix(alp, n+1, ntau)
  val_1 <- K %*% alp[2:(n+1),1]+alp[1,1]
  val_3 <- K %*% alp[2:(n+1),2]+alp[1,2]
  val_5 <- K %*% alp[2:(n+1),3]+alp[1,3]
  val_7 <- K %*% alp[2:(n+1),4]+alp[1,4]
  val_9 <- K %*% alp[2:(n+1),5]+alp[1,5]
  objfun(alp[,1],y,K,lam2,tau=0.1) + objfun(alp[,2],y,K,lam2,tau=0.3) +
        objfun(alp[,3],y,K,lam2,tau=0.5) + objfun(alp[,4],y,K,lam2,tau=0.7) +
        objfun(alp[,5],y,K,lam2,tau=0.9) + lam1 * sum(smooth_relu(val_1-val_3)) +
         lam1 * sum(smooth_relu(val_3-val_5)) + lam1 * sum(smooth_relu(val_5-val_7)) +
         lam1 * sum(smooth_relu(val_7-val_9))
}

