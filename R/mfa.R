


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
  
  refit <- function(.X, .ncomp=ncomp) { 
    mfa(.X, .ncomp, center=center, scale=scale, normalization=normalization, rank_k=rank_k) 
  }
  
  permute_refit <- function(.ncomp=ncomp) {
    Xperm <- block_apply(X, function(x) {
      sidx <- sample(1:nrow(x))
      x[sidx,]
    })
    
    refit(Xperm, .ncomp)
    
    
  }
  
  A <- rep(alpha, block_lengths(Xr))
  
  pca_fit <- genpca(unclass(Xr), 
                    A=A, 
                    ncomp=ncomp, 
                    center=FALSE, 
                    scale=FALSE)
  
  
  ncomp <- length(pca_fit$d)
  
  
  sc <- pca_fit$scores
  row.names(sc) <- row.names(X)
  
  result <- list(
    X=X,
    Xr=Xr,

    scores=sc,
    ntables=nblocks(X),
    pca_fit=pca_fit,
    center=center,
    scale=scale,
    ncomp=ncomp,
    block_indices=bind,
    alpha=alpha,
    normalization=normalization,
    refit=refit,
    permute_refit=permute_refit,
    table_names=names(X),
    reprocess=reprocess,
    reduced_rank=rank_k,
    is_reduced=is_reduced,
    permute_refit=permute_refit
  )
  
  class(result) <- c("mfa", "list")
  result
}

#' @export
loadings.musu_bada <- function(x, table_index=1:nrow(x$blockIndices), 
                               comp=1:ncol(x$pca_fit$v)) {
  
  lds <- loadings(x$pca_fit)
  block_matrix(lapply(table_index, function(i) {
    ind <- x$block_indices[[i]]
    lds[comp,ind]
  }))
}

#' @export
singular_values.mfa <- function(x) x$pca_fit$d


#' @export
scores.mfa <- function(x) {
  x$scores
}

#' @export
loadings <- function(x) {
  loadings.genpca(x$pca_fit)
}


#' @export
partial_scores.mfa <- function(x, table_index=x$ntables) {
  res <- lapply(table_index, function(i) {
    ind <- attr(Xr, "block_indices")[i,]
    x$ntables * project(pca_fit, get_block(Xr, i), subind=ind[1]:ind[2])
  })
  
  names(res) <- paste0("Block_", table_index)
  res
}

project.mfa <- function(x, newdata, ncomp=x$ncomp, pre_process=TRUE, table_index=1:x$ntables) {
  if (length(table_index) > 1) {
    assert_that(is.list(newdata) && length(newdata) == length(table_index))
  }
  
  if (is.vector(newdata)) {
    assert_that(length(table_index) == 1)
    newdata <- list(matrix(newdata, 1, length(newdata)))
  }
  
  if (is.matrix(newdata)) {
    assert_that(length(table_index) == 1)
    newdata <- list(newdata)
  }
  
  ## project new data-point(s)
  res <- lapply(1:length(table_index), function(i) {
    tbind <- table_index[i]
    xnewdat <- x$reprocess(newdata[[i]], tbind)
   
    x$ntables * project(x$pca_fit, xnewdat, 
                        ncomp=ncomp, subind=x$block_indices[[tbind]])
   
  })
  
  names(res) <- paste0("table_", table_index)
  res
} 

#' @export
reconstruct.mfa <- function(x, ncomp=x$ncomp) {
  recon <- reconstruct(x$pca_fit, ncomp)
}


#' @export
contributions.mfa <- function(x, type=c("column", "row", "table")) {
  contr <- contributions(x$pca_fit)
  type <- match.arg(type)
  
  if (type == "table") {
    out <- do.call(cbind, lapply(1:x$ncomp, function(i) {
      sapply(1:x$ntables, function(j) {
        ind <- x$block_indices[[j]]
        sum(contr[i,ind])
      })
    }))
      
  } else if (type == "row") {
    contributions(x$pca_fit, type="row")
  } else {
    contributions(x$pca_fit, type="column")
  }
  
}

#' @importFrom vegan procrustes 
procrusteanize.mfa <- function(x, ncomp=2) {
  F <- scores(x)[,1:ncomp]
  
  res <- lapply(1:length(x$Xlist), function(i) {
    xcur <- get_block(x$Xr[[i]])
    xp <- project(x, xcur, ncomp=ncomp, table_index=i)
    pres <- vegan::procrustes(F, xp)
    list(H=pres$rotation,
         scalef=pres$scale)
  })
  
  ret <- list(rot_matrices=res, ncomp=ncomp, musufit=x)
  class(ret) <- c("procrusteanized_mfa", "list")
  ret
}



## project from existing table
#' @export
predict.mfa <- function(x, newdata, ncomp=x$ncomp, table_index=1:x$ntables, pre_process=TRUE) {
  type <- match.arg(type)
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  fscores <- if (length(table_index) == x$ntables) {
    assert_that(ncol(newdata) == ncol(x$X))
    Reduce("+", lapply(table_index, function(i) {
      ind <- x$block_indices[[i]]
      
      Xp <- if (pre_process) {
        x$reprocess(newdata[, ind], i)
      } else {
        newdata[, ind]
      }
      
      project(x$pca_fit, Xp, ncomp=ncomp) * x$ntables
    }))
  } else if (length(table_index) == 1) {
    ind <- x$block_indices[[table_index]]
    
    Xp <- if (pre_process) {
      x$reprocess(newdata, table_index)
    } else {
      newdata[, ind]
    }
    
    fscores <- project(x$pca_fit, Xp, ncomp=ncomp) * x$ntables

  }
  
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

  
