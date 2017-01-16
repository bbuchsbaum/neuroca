

#' @importFrom rowr rollApply
temporal_adjacency <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  weight_mode <- match.arg(weight_mode)
  len <- length(time)
  
  wfun <- if (weight_mode == "binary") {
    function(x) 1
  } else {
    function(x) heat_kernel(x, sigma)
  }
  
  m <- do.call(rbind, rollApply(time, window=window, function(t) {
    if (length(t) > 1) {
      h <- wfun( (t[1] - t[-1])^2)
      cbind(t[1], t[-1], h)
    } else {
      NULL
    }
  }))
  
  
  sm <- sparseMatrix(i=m[,1], j=m[,2], x=m[,3], dims=c(len, len))
  (sm + t(sm))/2

}

temporal_laplacian <- function(time, weight_mode = c("heat", "binary"), sigma=1, window=2) {
  adj <- temporal_adjacency(time, weight_mode, sigma, window)
  Diagonal(x=rowSums(adj))  - adj
}


spatial_laplacian <- function(coord_mat, dthresh=1.42, nnk=27) {
  adj <- spatial_adjacency(coord_mat=coord_mat, dthresh=dthresh, nnk=nnk)
  Diagonal(x=rowSums(adj))  - adj
}


spatial_adjacency <- function(coord_mat, dthresh=1.42, nnk=27) {
  full_nn <- FNN::get.knnx(coord_mat, coord_mat, nnk)
  
  triplet <- do.call(rbind, lapply(1:nrow(coord_mat), function(i) {
    dvals <- full_nn$nn.dist[i,1:nnk]
    keep <- dvals > 0 & dvals < dthresh
    if (keep[2]) {
      cbind(i,full_nn$nn.index[i,keep], 1)
    } else {
      NULL
    }
  }))
  
  sparseMatrix(i=triplet[,1], j=triplet[,2], x=triplet[,3])
  
}


indices_to_sparse <- function(nn.index, hval, return_triplet=FALSE) {
  M <- do.call(rbind, lapply(1:nrow(nn.index), function(i) {
    cbind(i, nn.index[i,], hval[i,])
  }))
  
  if (return_triplet) {
    M
  } else {
    Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(nn.index), nrow(nn.index)))
  }
}



heat_kernel <- function(x, sigma=1) {
  exp(-x/(2*sigma^2))
}

normalized_heat_kernel <- function(x, sigma=1, len) {
  norm_dist <- (x^2)/(2*len)
  exp(-norm_dist/(2*sigma^2))
}

#' @param neighbor_mode
#' @param weight_mode binary (1 if neighbor, 0 otherwise), heat kernel, normalized heat kernel
#' @param k number of neighbors
#' @param sigma parameter for heat kernel \code{exp(-dist/(2*sigma^2))}
#' @param labels the class of the categories when \code{weight_mode} is \code{supervised}, supplied as a \code{factor} with \code{nrow(labels) == nrow(X)}
#' @export
construct_weight_matrix <- function(X, neighbor_mode=c("knn", "supervised"), weight_mode=c("heat", "normalized", "binary"), k=5, sigma=1, labels=NULL) {
  neighbor_mode = match.arg(neighbor_mode)
  weight_mode = match.arg(weight_mode)
  
  if (weight_mode == "normalized") {
    X <- t(scale(t(X), center=TRUE, scale=TRUE))
  }
  
  
  
  wfun <- if (weight_mode == "heat") {
    function(x) heat_kernel(x, sigma)
  } else if (weight_mode == "binary") {
    function(x) (x>0) * 1
  } else if (weight_mode == "normalized") {
    function(x) normalized_heat_kernel(x,sigma,ncol(X))
  }
  
  if (neighbor_mode == "knn") {
      W <- weighted_knn(X, k, FUN=wfun)
  } else if (neighbor_mode == "supervised") {
    assertthat::assert_that(!is.null(labels))
    labels <- as.factor(labels)
    if (k == 0) {
      W <- label_matrix(labels, labels)
    } else if (k > 0) {
      M <- do.call(rbind, lapply(levels(labels), function(lev) {
        idx <- which(labels == lev)
        M <- weighted_knn(X[idx,], k, FUN=wfun, return_triplet=TRUE)
        cbind(idx[M[,1]], idx[M[,2]], M[,3])
      }))
      
      Matrix::sparseMatrix(i=M[,1], j=M[,2], x=M[,3], dims=c(nrow(X), nrow(X)))
      
    }
    
  }
}




psparse <- function(..., fun=c("max", "min"), na.rm=FALSE, return_triplet=FALSE) {
  fun <- match.arg(fun)
  # check that all matrices have conforming sizes
  num.rows <- unique(sapply(list(...), nrow))
  num.cols <- unique(sapply(list(...), ncol))
  stopifnot(length(num.rows) == 1)
  stopifnot(length(num.cols) == 1)
  
  cat.summary <- do.call(rbind, lapply(list(...), summary))
  

  out.summary <- if (fun == "min") {
    aggregate(x ~ i + j, data = cat.summary, FUN=function(x) {
      if (length(x) == 1) 0 else x[1]
    })
  } else {
    aggregate(x ~ i + j, data = cat.summary, FUN=max, na.rm)
  }
  
  if (return_triplet) {
    cbind(i=out.summary$i, j=out.summary$j, x=out.summary$x)
  } else {
    sparseMatrix(i = out.summary$i,
                 j = out.summary$j,
                 x = out.summary$x,
                 dims = c(num.rows, num.cols))
  }
}

pmin.sparse <- function(..., na.rm = FALSE, return_triplet=FALSE) { psparse(..., fun="min", return_triplet=return_triplet) }

pmax.sparse <- function(..., na.rm = FALSE, return_triplet=FALSE) { psparse(..., fun="max", return_triplet=return_triplet ) }


#' @export
sim_knn_from_adj <- function(A, k=5, type=c("normal", "mutual"), ncores=1) {
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



#' weighted_knn 
#' 
#' @param X the data matrix, where rows are instances and columns are features
#' @param k the number of nearest neighbors
#' @param FUN the function used to convert euclidean distances to correlations
#' @param type
#' @importFrom assertthat assert_that
#' @importFrom FNN get.knn
#' @importFrom Matrix t
#' @export
weighted_knn <- function(X, k=5, FUN=heat_kernel, type=c("normal", "mutual"), return_triplet=FALSE) {
  assert_that(k > 0 && k <= nrow(X))
  
  type <- match.arg(type)
  nn <- FNN::get.knn(X, k=k)
  
  nnd <- nn$nn.dist

  hval <- FUN(nnd)
  
  W <- indices_to_sparse(nn$nn.index, hval)
  
  if (type == "normal") {
    psparse2(W, pmax, return_triplet=return_triplet)
  } else {
    psparse2(W, pmin, return_triplet=return_triplet)
  }
}

psparse2 <- function(W, FUN, return_triplet=FALSE) {
  ind <- which(W != 0, arr.ind=TRUE)
  x1 <- W[ind]
  x2 <- W[cbind(ind[,2], ind[,1])]
  
  if (return_triplet) {
    cbind(i=c(ind[,1],ind[,2]), j=c(ind[,2],ind[,1]), x=rep(FUN(x1,x2),2))
  } else {
    sm <- sparseMatrix(i=c(ind[,1],ind[,2]), j=c(ind[,2],ind[,1]), x=rep(FUN(x1,x2),2), dims=dim(W), use.last.ij=TRUE)
    #sm <- sparseMatrix(i=ind[,1], j=ind[,2], x=FUN(x1,x2), dims=dim(W), symmetric=TRUE)
    #as(sm, "dgTMatrix")
    sm
  }
}





#' label_matrix
#' 
#' @importFrom Matrix sparseMatrix
label_matrix <- function(a, b, type=c("s", "d"), return_matrix=TRUE) {
  type <- match.arg(type)
  
  a.idx <- which(!is.na(a))
  b.idx <- which(!is.na(b))
  
  sfun <- function(x, y) {
    if (is.na(x) || is.na(y)) {
      0
    } else if (x == y) {
      1
    } else {
      0
    }
  }
  
  dfun <- function(x, y) {
    if (is.na(x) || is.na(y)) {
      0
    } else if (x == y) {
      0
    } else {
      1
    }
  }
  
  fun <- if (type == "s") sfun else dfun
  
  out <- lapply(a.idx, function(i) {
    ret <- lapply(b.idx, function(j) {
      if (fun(a[i], b[j])) {
        c(i,j,1)
      } else {
        NULL
      }
    })
    ret[sapply(ret, function(x) !is.null(x))]
  })
  
  out <- unlist(out, recursive=FALSE)
  out <- do.call(rbind, out)
  
  if (return_matrix) {
    sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(length(a), length(b)))
  } else {
    out
  }
}