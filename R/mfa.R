


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
mfa <- function(X, ncomp=2, center=TRUE, scale=FALSE, 
                normalization=c("MFA", "RV", "None"), 
                rank_k=NULL) {
  
  assertthat::assert_that(inherits(X, "block_matrix"))
  
  normalization <- match.arg(normalization)
  
  X <- pre_processor(X, center=center,scale=scale)
  
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
  
  bind <- block_index_list(X)
  
  reprocess <- function(newdat, table_index) {
    ## given a new observation(s), pre-process it in the same way the original observations were processed
    
    prep <- attr(X, "pre_process")
    newdat <- prep(newdat, bind[[table_index]])
    
    if (is_reduced) {
      newdat <- project(reducer, newdat, table_index)
    }
    
    newdat
    
  }
  
  A <- rep(alpha, block_lengths(Xr))
  
  pca_fit <- genpca(unclass(Xr), 
                    A=A, 
                    ncomp=ncomp, 
                    center=FALSE, 
                    scale=FALSE)
  
  
  ncomp <- length(pca_fit$d)
  
  
  #partial_fscores = 
    #lapply(1:nblocks(X), function(i) {
    #  ind <- attr(Xr, "block_indices")[i,]
      #nblocks(X) * (get_block(Xr, i) * alpha[i]) %*% pca_fit$v[ind[1]:ind[2],]
    #  project(pca_fit, get_block(Xr, i), subind=ind[1]:ind[2])
    #})
  
  sc <- pca_fit$scores
  row.names(sc) <- row.names(X)
  
  result <- list(
    X=X,
    Xr=Xr,

    scores=sc,
    #partial_scores=partial_fscores,
    
    #table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
    #  sapply(1:nblocks(Xr), function(j) {
    #    ind <- bind[j,1]:bind[j,2]
    #    sum(pca_fit$v[ind,i]^2 * alpha[j])
    #  })
    #})),
    
    ntables=nblocks(X),
    pca_fit=pca_fit,
    center=center,
    scale=scale,
    ncomp=ncomp,
    blockIndices=bind,
    alpha=alpha,
    normalization=normalization,
    #refit=refit,
    #table_names=table_names,
    reprocess=reprocess,
    rank_k=rank_k,
    permute_refit=permute_refit
    #center_vec=unlist(Xreduced$centroids),
    #scale_vec=unlist(Xreduced$scales)
  )
  
  class(result) <- c("mfa", "list")
  result
}


#' @export
reconstruct.mfa <- function(x, ncomp=x$ncomp) {
  recon <- reconstruct(x$pca_fit, ncomp)
  x$reverse_pre_process(recon)
}


#' @export
impute_mfa <- function(X, ncomp=min(dim(X)), center=TRUE, scale=FALSE, 
                       normalization=c("MFA", "RV", "None"), 
                       iter_max=100, threshold = 1e-05, 
                       Xmiss_init=NULL) {
  
  
  nas <- is.na(X)
  missing_ind <- which(nas, arr.ind = TRUE)
  missing_1d <- which(nas)
  
  if (is.null(Xmiss_init)) {
    Xbar <- colMeans(X, na.rm=TRUE)
    Ximp <- X 
    Ximp[missing_ind] <- Xbar[missing_ind[,2]]
  } else {
    assertthat::assert_that(all(dim(Xmiss_init) == dim(X)))
    xtmp[missing_ind] <- 0
    Ximp <- xtmp + Xmiss_init
  }
  
  old <- Inf
  iter <- 1
  
  while(iter < iter_max & criterion > threshold) {
    mres <- mfa(Ximp, ncomp=ncomp, center=center, scale=scale, normalization=normalization)
    recon <- reconstruct(mres)
    
    Xnew <- X
    Xnew[missing_ind] <- recon[missing_ind]
    
    obj <- sum((X[-missing_1d] - recon[-missing_1d])^2)
   
    Ximp <- Xnew
    criterion <- abs(1 - obj/old)
    print(criterion)
    old <- obj
    iter <- iter+1
  }
  
  Ximp
}

  
