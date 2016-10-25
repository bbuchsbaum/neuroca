
library(RSpectra)
library(FNN)
library(assertthat)

psparse <- function(..., fun=max, na.rm=FALSE) {
  # check that all matrices have conforming sizes
  num.rows <- unique(sapply(list(...), nrow))
  num.cols <- unique(sapply(list(...), ncol))
  stopifnot(length(num.rows) == 1)
  stopifnot(length(num.cols) == 1)
  
  cat.summary <- do.call(rbind, lapply(list(...), summary))
  out.summary <- aggregate(x ~ i + j, data = cat.summary, FUN=fun, na.rm)
  
  sparseMatrix(i = out.summary$i,
               j = out.summary$j,
               x = out.summary$x,
               dims = c(num.rows, num.cols))
}

pmin.sparse <- function(..., na.rm = FALSE) { psparse(..., fun=min) }

pmax.sparse <- function(..., na.rm = FALSE) { psparse(..., fun=max ) }

simgraph_knn_from_adj <- function(A, k=5, type=c("normal", "mutual"), ncores=1) {
  assert_that(k > 0 && k <= nrow(X))
  
  type <- match.arg(type)
  
  jind <- 1:nrow(X)
  A2 <- do.call(rbind, mclapply(1:nrow(A), function(i) {
    ord <- order(A[i,], decreasing=TRUE)
    cbind(i=i, j=jind[ord[1:k]],x=A[i, ord[1:k]])
  }, mc.cores=ncores))
  
  m <- sparseMatrix(i=A2[,1], j=A2[,2], x=A2[,3], dims=c(nrow(A), nrow(A)))
  m <- if (type == "normal") {
    pmax.sparse(m, t(m))
  } else {
    pmin.sparse(m, t(m))
  }
}  
  
simgraph_knn <- function(X, k=5, sigma=1, type=c("normal", "mutual"), ncores=1) {
  assert_that(k > 0 && k <= nrow(X))
  
  type <- match.arg(type)
  nn <- get.knnx(X, X, k=k+1)
  
  jind <- 1:nrow(X)
  A <- do.call(rbind, mclapply(1:nrow(X), function(i) {
    cbind(i=i, j=jind[nn$nn.index[i,2:k]],x=exp(-nn$nn.dist[i,2:k]/sigma^2))
  }, mc.cores=ncores))
  
  m <- sparseMatrix(i=A[,1], j=A[,2], x=A[,3], dims=c(nrow(X), nrow(X)))
  if (type == "normal") {
    pmax.sparse(m, t(m))
  } else {
    pmin.sparse(m, t(m))
  }
}


commute_time <- function(A, ncomp=nrow(A)) {
  D <- rowSums(A)
  Dtilde <- 1/(sqrt(D))
  
  #P <- diag(1/D) %*% A
  
  M <- diag(Dtilde) %*% A %*% diag(Dtilde)
  
  decomp <- eigs(M, k=ncomp)
  pii <- D/sum(A)
  v <- decomp$vectors[, 2:ncomp]
  
  cds <- sweep(v, 2, sqrt(1 - decomp$values[2:ncomp]), "/")
  cds <- sweep(cds, 2, 1/pi)
  
}