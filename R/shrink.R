SURE_cost <- function(n,p,
                      method.ada,
                      sigma,
                      gamma,
                      lambda,
                      svdXd){
  
  if(method.ada == "qut") {
    gamma <- exp(gamma) + 1
  } else { 
    lambda <- exp(lambda)
  }  
  
  DD  <- svdXd
  DD2 <- DD^2
  lDD <- length(DD)
  D <-  matrix(0, lDD, lDD)
  dhat <- DD * pmax(1-(lambda/DD)^gamma, 0) 
  temp <- DD * dhat
  for (i in 1:lDD){
    DD2i <- DD2[i]
    diff2i <- DD2[i]-DD2
    diff2i[i] <- Inf
    D[i, ] <- temp[i]/diff2i
  }
  
  gradd <- (1+(gamma-1)*(lambda/svdXd)^gamma) * (svdXd >= lambda)
  DIV <- sum(gradd + abs(n-p)*dhat/svdXd) + 2*sum(D)
  
  if(method.ada == "gsure") {
    mse.ada <- sum((dhat-svdXd)^2)/(1-DIV/n/p)^2
  } else{       
    mse.ada <- -n*p*sigma^2 + sum((dhat-svdXd)^2)  + 2*sigma^2 *DIV
  }
  return(mse.ada)
}


#' @examples 
#' X <- matrix(rnorm(10*10), 10, 10)
#' x <- pca(X, ncomp=10)
shrink.pca <- function(x, method=c("GSURE", "QUT", "SURE"), sigma, lambda=NULL, gamma=NULL,
                       nbsim = 500, method.optim = "BFGS", gamma.seq = seq(1, 5, by=.1), lambda0=NA) {
  method <- match.arg(method, c("GSURE", "Gsure", "gsure", 
                                "gSure", "GSure", "QUT", "Qut", "qut", "SURE", "sure", 
                                "Sure"), several.ok = T)[1]
  
  if (!is.na(sigma) & (sigma <= 0)) {
    stop("sigma must be positive")
  }
  
  if (sum(gamma.seq < 1) > 0) {
    stop("the gammas may be greater or equal to 1")
  }
  if (length(gamma.seq) == 1) {
    if (gamma.seq == 1) 
      warning("you are performing soft singular values thresholding")
  }
  
  svdXd <- x$d
  
  if(is.na(lambda0)){
    lambda0 <- log(median(svdXd))
  }

  method <- tolower(method)
  
  val.optglob <- Inf
  n <- nrow(scores(x))
  p <- nrow(loadings(x))
  
  if (method == "qut") {
    nbsim <- nbsim
    maxd <- rep(1, nbsim)
    for (i in 1:nbsim) {
      maxd[i] <- max(RSpectra::svds(matrix(rnorm(n * p, sd = sigma), n, p), k=1)$d)
    }
    lambda.o <- quantile(maxd, 1 - 1/sqrt(log(max(n, p))))
    if (length(gamma.seq) != 1) {
      res.opti <- optim(par=0, SURE_cost, n=n, p=p, method.ada = method, 
                        sigma = sigma,
                        lambda = lambda.o, 
                        svdXd = svdXd,
                        method = method.optim)
      gamma.o <- exp(res.opti$par) + 1
    } else {
      gamma.o <- gamma.seq
    }
  } else {
    if (is.na(lambda0)) {
      lambda0 <- log(median(svdXd))
    }
    for (gamma in gamma.seq) {
      res.opti <- optim(par=lambda0, fn=SURE_cost, n=n, p=p, sigma = sigma, 
                        gamma = gamma, svdXd = svdXd, method.ada = method, 
                        method = method.optim)
      lambda.temp <- exp(res.opti$par)
      val.opt <- res.opti$val
      if (val.opt < val.optglob) {
        gamma.o <- gamma
        lambda.o <- lambda.temp
        val.optglob <- val.opt
      }
    }
  }
  
  dhat <- svdXd * pmax( 1-(lambda.o/svdXd)^gamma.o, 0) 
  list(dhat=dhat,
       lambda=lambda.o,
       gamma=gamma.o)
  
}

#' shrink_pca
#' 
#' adaptive shrinkage pca from the \code{denoiseR} package
#' 
#'   
#' @param X
#' @param center
#' @param scale
#' @importFrom denoiseR adashrink
#' @export
shrink_pca <- function(X, preproc=center(), method = c("GSURE", "QUT", "SURE"), ...) {
  assert_that(is.matrix(X) || inherits(X, "Matrix"))
  
  procres <- prep(preproc, X)
  Xp <- procres$init(X)
  
  res <- denoiseR::adashrink(Xp, method=method, center=FALSE, ...)
  
  keep <- res$singval > 1e-06
  
  if (sum(keep) == 0) {
    keep <- 1
  } 
  
  v=res$low.rank$v[,keep,drop=FALSE]
  u=res$low.rank$u[,keep,drop=FALSE]
  d=res$low.rank$d[keep]
  
  ret <- bi_projector(
    preproc=procres,
    ncomp=length(d),
    v=v, 
    u=u,
    d=d,
    scores=t(t(as.matrix(u)) * d),
    classes=c("shrink_pca", "pca"))
  
}


