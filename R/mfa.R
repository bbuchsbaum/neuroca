compute_sim_mat <- function(block_mat, FUN, ...) {
  pairs <- combn(nblocks(block_mat),2)
  M <- matrix(0, nblocks(block_mat), nblocks(block_mat))
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    sim <- FUN(get_block(block_mat, p1), get_block(block_mat, p2), ...)
    M[p1,p2] <- sim
    M[p2,p1] <- sim
  }
  
  M
}


coefficientRV <- function(X, Y, center=TRUE, scale=FALSE) {
  if (dim(X)[[1]] != dim(Y)[[1]]) 
    stop("no the same dimension for X and Y")
  if (dim(X)[[1]] == 1) {
    print("1 configuration RV is  NA")
    rv = NA
  }
  else {
    if (center || scale) {
      Y <- scale(Y, center=center, scale=scale)
      X <- scale(X, center=center, scale=scale)
    }
    W1 <- X %*% t(X)
    W2 <- Y %*% t(Y)
    rv <- sum(diag(W1 %*% W2))/(sum(diag(W1 %*% W1)) * 
                                  sum(diag(W2 %*% W2)))^0.5
  }
  return(rv)
}

normalization_factors <- function(block_mat, type=c("MFA", "RV", "None")) {
  type <- match.arg(type)
  
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
  } else if (type == "RV" && nblocks(block_mat) > 2) {
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    wts <- abs(svd_wrapper(smat, ncomp=1, method="propack")$u[,1])
  } else {
    rep(1, nblocks(block_mat))
  }
}

#' @export
mfa <- function(X, ncomp=2, center=TRUE, scale=FALSE, normalization=c("MFA", "RV", "None"), rank_k=NULL) {
  assertthat::assert_that(inherits(X, "block_matrix"))
  
  Xr <- if (!is.null(rank_k)) {
    is_reduced <- TRUE
    reducer <- reduce_rank(X, rank_k)
    reducer$x
  } else {
    is_reduced=FALSE
    reducer <- NULL
    X
  }
  
  alpha <- normalization_factors(Xr, type=normalization)
  bind <- attr(Xr, "block_ind")
  
  reprocess <- function(newdat, table_index) {
    ## given a new observation(s), pre-process it in the same way the original observations were processed
    newdat <- Xreduced$pre_process_funs[[table_index]](newdat)
    
    if (is_reduced) {
      newdat <- project(reducer, newdat, table_index)
    }
    
    newdat
    
  }
  
  A <- rep(alpha, block_lengths(Xr))
  
  pca_fit <- genpca(Xr, A=A, 
                    ncomp=ncomp, 
                    center=FALSE, 
                    scale=FALSE)
  
  
  ncomp <- length(pca_fit$d)
  
  
  partial_fscores = 
    lapply(1:nblocks(X), function(i) {
      ind <- attr(Xr, "block_indices")[i,]
      #nblocks(X) * (get_block(Xr, i) * alpha[i]) %*% pca_fit$v[ind[1]:ind[2],]
      project(pca_fit, get_block(Xr, i), subind=ind[1]:ind[2])
    })
  
  sc <- pca_fit$scores
  row.names(sc) <- row.names(X)
  
  
}