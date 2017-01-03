

cobec <- function(Xlist, ncomp=3, pcadim=nrow(Xlist[[1]]), maxiter=50, ctol=1e-3) {
  
  N <- length(Xlist)
  NR <- nrow(Xlist[[1]])
  
  J <- rep(0,N)
  U <- vector(N, mode="list")
  
  for (n in 1:N) {
    edecomp <- eig(t(Xlist[[n]]) %*% Xlist[[n]])
    D <- diag(edecomp$values)
    flag <- D > 1e-3
    U[[n]]= Xlist[[n]] %*% edecomp$values[,flag] * diag(D[flag] ^-.5)
    
    J[n] <- ncol(U[[n]])
    
  }
  
  idx <- order(J)
  
  Ac <- ccak_init(U[[idx[[1]]]],U[[idx[[2]]]],ncomp)
  x <- vector(N, mode="list")
  
  for (it in 1:maxiter) {
    c0 <- Ac
    c2 <- matrix(0, NRows, ncomp)
    for (n in 1:N) {
      x[[n]] <- t(U[[n]]) %*% Ac
      c2 <- c2 + U[[n]] %*% x[[n]]
    }
    
    if (mean(abs(diag(t(Ac)  %*% c0))) >1-ctol) {
        break
    }
  }
  
  Ac
  
}

#' @importFrom pracma orth
#' @importFrom RSpectra eigs
ccak_init(X,Y,K) <- function(X, Y, K=NULL) {
  #CCA calculate canonical correlations with K largest canonical correlation
  dimX=dim(X)
  dimY=dim(Y)
  
  Ix <- dimX[1]
  Iy <- dimY[1]
  
  Nx <- dimX[2]
  Ny <- dimY[2]
  
  if (is.null(K)) {
    K=min(Nx,Ny)
  }

  canres <- cancor(X,Y)
  #[Wx,Wy,r]=canoncorr(X,Y);
  # [r,pos]=sort(r,'descend');
  Wx <- canres$xcoef
  Wy <- canres$ycoef
  
  Wx=Wx[,1:K]
  Wy=Wy[,1:K]
  
  Wx=pracma::orth(X %*% Wx)
  Wy=pracma::orth(Y %*% Wy)
  
  C <- diag(2*K)
  
  C[1:K,(K+1):(2*K)] <- t(Wx) %*% Wy
  C[(K+1):(2*K),1:K] <- t(Wy) %*% Wx
  
  eres <- eigs(C,k=K,'LM')
  u <- eres$vectors %*% diag(eres$values)^.5
  proj=Wx %*% u[1:K,] + Wy %*% u[K+1:2*K,,]
  
  proj=pracma::orth(proj)
}
  
  
  