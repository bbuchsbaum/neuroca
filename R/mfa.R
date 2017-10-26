


normalization_factors <- function(block_mat, type=c("MFA", "RV", "RV-MFA", "None")) {
  type <- match.arg(type)
  
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
  } else if (type == "RV" && nblocks(block_mat) > 2) {
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    wts <- abs(svd_wrapper(smat, ncomp=1, method="propack")$u[,1])
  } else if (type == "RV-MFA") {
    alpha1 <- unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    smat[which(smat < 0)] <- 0
    diag(smat) <- 1
    alpha2 <- abs(svd_wrapper(smat, ncomp=1, method="propack")$u[,1])
    alpha1*alpha2
    
  } else {
    rep(1, nblocks(block_mat))
  }
}


 
#' @export
mfa <- function(X, ncomp=2, center=TRUE, scale=FALSE, 
                normalization=c("MFA", "RV", "None", "RV-MFA"), 
                rank_k=NULL) {
  
  assertthat::assert_that(inherits(X, "block_matrix"))
  
  normalization <- match.arg(normalization)
  
  ## pre-preocess full data matrix
  X <- pre_processor(X, center=center,scale=scale)
  
  ## if required, compute low-rank matrix
  Xr <- if (!is.null(rank_k)) {
    is_reduced <- TRUE
    reducer <- reduce_rank(X, rank_k)
    reducer$x
  } else {
    is_reduced=FALSE
    reducer <- NULL
    X
  }
  
  ## compute nblock normalization factors
  alpha <- normalization_factors(Xr, type=normalization)
  
  ## construct vector of nomalization weights expanded over variables per block
  A <- rep(alpha, block_lengths(Xr))
  
  ## get the block indices 
  bind <- block_index_list(X)

 
  ## compute generalized pca
  fit <- genpca(unclass(Xr), 
                    A=A, 
                    ncomp=ncomp, 
                    center=FALSE, 
                    scale=FALSE)
  
  
  result <- list(
    ## original data matrix
    X=X,
    ## reduced data matrix (equal to X if no data reduction)
    Xr=Xr,
    nvars=ncol(X),
    ncomp=ncomp,
    ntables=nblocks(X),
    fit=fit,
    center=center,
    scale=scale,
    block_indices=bind,
    alpha=alpha,
    normalization=normalization,
    table_names=names(X),
    reduced_rank=rank_k,
    is_reduced=is_reduced,
    reducer=reducer
  )
  
  class(result) <- c("mfa", "multiblock", "list")
  result
}


refit.mfa <- function(x, X, ncomp=x$ncomp) { 
  mfa(X, ncomp, center=x$center, scale=x$scale, normalization=x$normalization, rank_k=x$rank_k) 
}

permute_refit.mfa <- function(ncomp=x$ncomp) {
  Xperm <- block_apply(x$X, function(x) {
    sidx <- sample(1:nrow(x))
    x[sidx,]
  })
  
  refit(x, Xperm, ncomp)
}


#' @export
singular_values.mfa <- function(x) x$fit$d

#' @export
reconstruct.mfa <- function(x, comp=1:ncomp(x)) {
  recon <- reconstruct(x$fit, comp)
}

#' @importFrom vegan procrustes 
procrusteanize.mfa <- function(x, ncomp=2) {
  F <- scores(x)[,1:ncomp]
  
  res <- lapply(1:length(x$Xlist), function(i) {
    xcur <- get_block(x$Xr[[i]])
    xp <- project(x, xcur, comp=1:ncomp, table_index=i)
    pres <- vegan::procrustes(F, xp)
    list(H=pres$rotation,
         scalef=pres$scale)
  })
  
  ret <- list(rot_matrices=res, ncomp=ncomp, musufit=x)
  class(ret) <- c("procrusteanized_mfa", "list")
  ret
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
  criterion <- Inf
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


  
