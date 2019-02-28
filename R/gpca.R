
.msqrt <- function(a) {
  a.eig <- eigen(a)
  a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
}

#' Generalized Principal Components Analysis
#' 
#' Compute a PCA in a inner-product space defined by row and coulmn constraint matrices.
#' 
#' @inheritParams pca
#' @param A the column constraints. Can be a \code{vector}, symmetric \code{matrix}, or symmetric sparse matrix with \code{ncol(X)} rows and columns.
#' @param M the row constraints. Can be a \code{vector}, symmetric \code{matrix}, or symmetric sparse matrix with \code{nrow(X)} rows and columns.
#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseMatrix t isDiagonal
#' @export
#' 
#' @references 
#' 
#' Abdi, H. (2007). Singular value decomposition (SVD) and generalized singular value decomposition. \emph{Encyclopedia of measurement and statistics}, 907-912.
#' 
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014). A generalized least-square matrix decomposition. \emph{Journal of the American Statistical Association}, 109(505), 145-159.
#' 
#' 
#' @examples 
#' 
#' coords <- expand.grid(x=seq(1,100), y=seq(1,100))
#' img <- apply(coords, 1, function(x) {
#'   x1 <- 1 - pnorm(abs(x[1] - 50), sd=8)
#'   x2 <- 1 - pnorm(abs(x[2] - 50), sd=8)
#'   x1*x2
#' })
#' 
#' mat <- matrix(img, 100,100)
#' mlist <- replicate(50, as.vector(mat + rnorm(length(mat))*.8), simplify=FALSE)
#' X <- do.call(rbind, mlist)
#' 
#' ## spatial smoother
#' S <- neighborweights:::spatial_smoother(coords, sigma=8, nnk=27)
#' gp1 <- genpca(X, A=S, ncomp=2)
#' 
#' Xs <- do.call(rbind, lapply(1:nrow(X), function(i) X[i,,drop=FALSE] %*% S))
#' gp2 <- genpca(as.matrix(Xs), ncomp=2, center=FALSE)
#' 
#' ## use an adjacency matrix to weight items sharing an index.
#' 
#' X <- matrix(rnorm(50*100), 50, 100)
#' colind <- rep(1:10, 10)
#' S <- neighborweights:::spatial_adjacency(as.matrix(colind), dthresh=4, sigma=1, nnk=27, normalized=TRUE, include_diagonal=TRUE, weight_mode="heat")
#' diag(S) <- 1
#' S <- S/RSpectra::svds(S,k=1)$d
#' gp1 <- genpca(X, A=S, ncomp=2)
genpca <- function(X, A=NULL, M=NULL, ncomp=min(dim(X)), 
                   preproc=center(), deflation=FALSE) {
  
 
  if (is.null(A)) {
    A <- sparseMatrix(i=1:ncol(X), j=1:ncol(X),x=rep(1, ncol(X)))
  }
  
  if (is.null(M)) {
    M <- sparseMatrix(i=1:nrow(X), j=1:nrow(X),x=rep(1,nrow(X)))
  }
  
  if (is.vector(A)) {
    assert_that(length(A) == ncol(X))
    A <- sparseMatrix(i=1:length(A), j=1:length(A),x=A)
  } else {
    assert_that(nrow(A) == ncol(X), msg=paste("nrow(A) != ncol(X) -- ", nrow(A), " != ", ncol(X)))
    assert_that(ncol(A) == ncol(X), msg=paste("ncol(A) != ncol(X) -- ", ncol(A), " != ", ncol(X)))
  }
  
  if (is.vector(M)) {
    assert_that(length(M) == nrow(X))
    M <- sparseMatrix(i=1:length(M), j=1:length(M),x=M)
  } else {
    assert_that(nrow(M) == nrow(X))
    assert_that(ncol(M) == nrow(X))
  }
  
  procres <- prep(preproc, X)
  Xp <- procres$Xp
  
  assert_that(ncomp > 0)
  ncomp <- min(min(dim(Xp)), ncomp)
  
  n = nrow(Xp)
  p = ncol(Xp)
  
  
  if (deflation) {
    if (n < p) {
      svdfit <- gmd_deflation_cpp(t(X), A, M, ncomp)
      svdfit$d <- svdfit$d[,1]
    } else {
      svdfit <- gmd_deflation_cpp(X, M, A, ncomp)
      svdfit$d <- svdfit$d[,1]
    }
  } else { 
    if(n < p){
      ret = gmdLA(t(Xp), A,M, ncomp,p,n)
      svdfit = list(u=ret$v, v=ret$u,d=ret$d, cumv=ret$cumv,propv=ret$propv)
    } else{
      svdfit = gmdLA(Xp, M,A,ncomp,n,p)
    }
  }
  
  scores <- t(t(as.matrix(M %*% svdfit$u)) * svdfit$d)
  col_scores <- t(t(as.matrix(A %*% svdfit$v)) * svdfit$d)
  
  #scores <- t(t(as.matrix(svdfit$u)) * svdfit$d)
  row.names(scores) <- row.names(X)
  #norm_loadings <- t(t(as.matrix(svdfit$v)) * svdfit$d)
  
  ret <- bi_projector(
              preproc=procres,
              ncomp=length(svdfit$d),
              v=svdfit$v, 
              u=svdfit$u, 
              d=svdfit$d, 
              scores=scores, 
              classes=c("genpca", "pca"),
              col_scores=col_scores,
              A=A,
              M=M)

  ret
}


#' @export
loadings.genpca <- function(obj) {
  obj$A %*% obj$v
}

project_xav <- function(X, A, V) {
  if (is.vector(A)) {
    t(t(X) * A) %*% V
  } else {
    X %*% A %*% V
  }
}
  
#' @export
predict.genpca <- function(x, newdata, comp=1:x$ncomp, pre_process=TRUE) {
  Xsup <- if (pre_process) {
    reprocess(x, newdata)
  } else {
    newdata
  }
  
  project_xav(Xsup, x$A, x$v[,comp,drop=FALSE])
  
}

reconstruct.genpca <- function(x, newdata=NULL, comp=1:x$ncomp, subind=NULL) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  
  if (is.null(subind)) {
    reverse_pre_process(x$preproc, newdata %*% t(loadings(x)[,comp,drop=FALSE]))
  } else {
    reverse_pre_process(x$preproc, newdata %*% t(loadings(x)[,comp,drop=FALSE])[,subind], subind=subind)
  }
}


#' @export
contributions.genpca <- function(x, type=c("column", "row")) {
  type <- match.arg(type)
  if (type == "column") {
    t(x$v^2) %*% x$A
  } else {
    t(x$u^2) %*% x$M
    
  }
}

#' @export
singular_values.genpca <- function(x) {
  x$d
}

#' @export
project.genpca <- function(x, newdata, comp=1:x$ncomp, pre_process=TRUE, colind=NULL) {
  if (is.null(newdata)) {
    return(scores(x)[,comp])
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(colind)) {
    Xsup <- if (pre_process) reprocess(x, newdata) else newdata
    project_xav(Xsup, x$A, x$v[,comp,drop=FALSE])
  } else {
    assertthat::assert_that(length(colind) == ncol(newdata))
    Xsup <- if (pre_process) {
      reprocess(x, newdata, colind)
    } else {
      newdata
    }
    
    comp <- comp[which(comp <= x$ncomp)]
    project_xav(Xsup, x$A[colind,colind,drop=FALSE], x$v[colind,comp,drop=FALSE])
  }
}


#' @export
project_cols.genpca <- function(x, newdata, comp=1:ncomp(x)) {
  
  ## if no newdata, then simply return the loadings
  if (missing(newdata)) {
    return(loadings(x)[,comp, drop=FALSE])
  } else {
    if (is.vector(newdata)) {
      newdata <- as.matrix(newdata, ncol=1)
    }
    assert_that(nrow(newdata) == nrow(x$u))
    t(newdata) %*% x$M %*% (x$u[,comp,drop=FALSE]) %*% diag(1/x$d[comp], nrow=length(comp), ncol=length(comp))
  }
}


#' @keywords internal
gmdLA <- function(X, Q, R, k=min(n,p), n, p) {
  ##computation

  if (isDiagonal(R)) {
    Rtilde <- Matrix::Diagonal(x=sqrt(Matrix::diag(R)))
    Rtilde.inv = Matrix::Diagonal(x=1/sqrt(Matrix::diag(R)))
  } else {
    decomp <- eigen(R)
    
    if (length(decomp$values) > 1) {
      v <- decomp$values
      if (sum(v < 0) > 1) {
        warning(paste("genpca: removing ", sum(v<0), 
                      " negative eigenvalues when computing inverse of constraint matrix"))
      }
    }
    
    keep <- which(decomp$values > 0 & abs(decomp$values) > 1e-7)
  
    decomp$vectors <- decomp$vectors[, keep]
    decomp$values <- decomp$values[keep]
  
    Rtilde <- decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
  
    inv.values = 1 / sqrt(decomp$values)
    Rtilde.inv = decomp$vectors %*% diag(inv.values) %*% t(decomp$vectors)
  }
  
  
  inmat <- Matrix::crossprod(X, Q) %*% X
  
  RtinRt <- Rtilde %*% inmat %*% Rtilde
  
  XR <- X %*% R
  ## nXn * n*n * pXp
  RnR <- R %*% inmat %*% R
  
  xtilde.decomp <- if (k == min(n,p)) {
    eigen(RtinRt)
  } else {
    #RSpectra::eigs_sym(RtinRt, k=k)
    ret <- RSpectra::eigs_sym(RtinRt, k=k+1)
    list(vectors=ret$vectors, values=ret$values)
  }

  keep <- which(abs(xtilde.decomp$values) > 1e-7)
  k <- min(k, length(keep))
  xtilde.decomp$vectors <- xtilde.decomp$vectors[, 1:k]
  xtilde.decomp$values <- xtilde.decomp$values[1:k]
  
  #Rtilde.inv %*% xtilde.decomp$vectors
  
  vgmd <- Rtilde.inv %*% xtilde.decomp$vectors
  dgmd <- sqrt(xtilde.decomp$values[1:k])
  ugmd <- matrix(nrow = n, ncol = k)
  cumv <- rep(0, k)
  propv <- dgmd ^ 2 / sum(diag(as.matrix(inmat %*% R)))
  normalizing.number <- 1
  
  for (i in 1:k) {
    normalizing.number = sqrt(vgmd[, i] %*% RnR %*% vgmd[, i])
    ugmd[, i] = as.matrix(XR %*% vgmd[, i]) / as.double(normalizing.number)
    cumv[i] = sum(propv[1:i])
  }
  
  list(
    u = ugmd[, 1:k, drop=FALSE],
    v = vgmd[, 1:k, drop=FALSE],
    d = dgmd[1:k],
    k = k,
    cumv = cumv,
    propv = propv
  )
  
}


#' @export
truncate.genpca <- function(obj, ncomp) {
  if (ncomp >= obj$ncomp) {
    warning("number of components to keep is greater than or equal to rank of pca fit, returning original model")
    ret <- obj
  } else {
    ret <- list(v=obj$v[,1:ncomp], 
                u=obj$u[,1:ncomp], 
                d=obj$d[1:ncomp], 
                scores=obj$scores[,1:ncomp], 
                ncomp=ncomp, svd.method=obj$svd.method, 
                pre_process=obj$pre_process)
    class(ret) <- c("genpca", "pca", "projector", "list")
  }
  
  ret
}


# mmult <- function(X, q) {
#   if (is.vector(q)) {
#     t(t(X) * q)
#   } else {
#     X %*% q
#   }
# }
# 
# mmult2 <- function(q, X) {
#   if (is.vector(q)) {
#     t(q * t(X))
#   } else {
#     q %*% X
#   }
#   
# }
# 
# cprod <- function(X, q) {
#   if (is.vector(q)) {
#     t(X * q)
#   } else {
#     Matrix::crossprod(X, q)
#   }
# }
