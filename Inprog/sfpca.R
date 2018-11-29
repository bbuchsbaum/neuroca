

soft_threshold <- function(x,d) return(sign(x)*pmax(0, abs(x)-d))


prox <- function(X, u, Lu, v, Su, lambda, conv=1e-5) {
  u0 <- u
  crit <- Inf
  while (crit > conv) {

    ut <- u + 1/Lu * (X %*% v - Su %*% u)
    ut <- soft_threshold(ut, lambda/Lu)
    
    if (sum(ut^2) > 0) {
      ut <- ut/sqrt(t(ut) %*% Su %*% ut)[1,1]
    } else {
      return(ut)
    }
    crit <- sum( (ut - u)^2)
    u <- ut
  }
  
  u
}

#' sparse and smooth functional principal components analysis
#' 
#' @param X the data matrix
#' @param ncomp number of components to estimate
#' @param lu lambda sparsity parameter for the rows
#' @param lv lambda sparsity parameter for the columns
#' @param au smoothing parameter for the rows
#' @param av smoothing parameter for the columns
#' @param Pu optional penalty matrix for rows
#' @param Pv optional penalty matrix for columns
#' 
#'
#' gen_mat <- function(sd=.2) {
#'   X <- matrix(0,64,64)
#'   X[32,32] <- 1
#'   X[16,16] <- 1
#'   S <- neighborweights::spatial_smoother(as.matrix(expand.grid(1:64, 1:64)), sigma=8,nnk=80)
#'   x=S %*% as.vector(X) 
#'   x + rnorm(length(x), sd=sd)
#' }
#' 
#' X <- replicate(50, gen_mat(.02))
#' X <- do.call(cbind,X)
#' X <- t(X)
#' 
#' res <- sfpca(X, ncomp=1, av=.5, lv=.5, coords=as.matrix(expand.grid(1:64, 1:64)))
sfpca <- function(X, ncomp=2, lu=.1, lv=1, au=.1, av=.1, coords, Pu=NULL, Pv=NULL, conv=1e-5) {


  if (is.null(Pv)) {
    Pv <- neighborweights::spatial_laplacian(coords, nnk=9, weight_mode="binary")
  }
  

  Sv <- Matrix::Diagonal(n=ncol(X)) + av*Pv
  Lv <- RSpectra::eigs(Sv,1)$values + .01

  if (is.null(Pu)) {
    Pu <- neighborweights::temporal_laplacian(1:nrow(X))
  }
  
  Su <- Matrix::Diagonal(n=nrow(X)) + au*Pu
  Lu <- RSpectra::eigs(Su,1)$values + .01
  
  init_svd <- RSpectra::svds(X, 1)
  v <- init_svd$v[,1]
  u <- init_svd$u[,1]
  criterion <- Inf
  Xhat <- X
  
  vout <- matrix(0, length(v), ncomp)
  uout <- matrix(0, length(u), ncomp)
  dout <- length(ncomp)
  for (i in 1:ncomp) {

    while (criterion > conv) {
      u_new <- prox(Xhat, u, Lu, v, Su, lu)
      while (all(u_new == 0)) {
        lu <- lu * .8
        message("new lu: ", lu)
        u_new <- prox(Xhat, u, Lu, v, Su, lu)
      }
      
      v_new <- prox(t(Xhat), v, Lv, u, Sv, lv)
      while (all(v_new == 0)) {
        lv <- lv * .8
        message("new lv: ", lv)
        v_new <- prox(t(Xhat), v, Lv, u, Sv, lv)
      }
      
      #print(sum(u_new^2))
      #print(sum(v_new^2))
    
      #v_new <- v_new/sqrt(t(v_new) %*% Sv %*% v_new)[1,1]
      #u_new <- u_new/sqrt(t(u_new) %*% Su %*% u_new)[1,1]
    
      criterion <- max(sum( (v_new - v)^2), sum((u_new - u)^2))
      print(criterion)
      u <- u_new
      v <- v_new
    }
  
    v <- v/sum(v^2)
    u <- u/sum(u^2)
    d <- t(u) %*% Xhat %*% v
    
    vout[,i] <- as.vector(v)
    uout[,i] <- as.vector(u)
    dout[i] <- d[i,1]
    #Xhat = Xhat - d %*% u  %*% t(v)
    dout[i] * (u  %*% t(v))
  }
  
  
  list(d=dout, v=vout, u=uout)
}