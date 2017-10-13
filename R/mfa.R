


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
  
  reprocess <- function(newdat, table_index=NULL) {
    ## given a new observation(s), pre-process it in the same way the original observations were processed
    
    prep <- attr(X, "pre_process")
    
    newdat <- if (is.null(table_index)) {
      prep(newdat)
    } else {
      prep(newdat, bind[[table_index]])
    }
    
    if (is_reduced && !is.null(table_index)) {
      project(reducer, newdat, table_index)
    } else if (is_reduced & is.null(table_index)) {
      project(reducer, newdat)
    } else {
      newdat
    }

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
  
  fit <- genpca(unclass(Xr), 
                    A=A, 
                    ncomp=ncomp, 
                    center=FALSE, 
                    scale=FALSE)
  
  
  result <- list(
    X=X,
    Xr=Xr,
    ntables=nblocks(X),
    fit=fit,
    center=center,
    scale=scale,
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
  
  class(result) <- c("mfa", "multiblock", "list")
  result
}


#' @export
singular_values.mfa <- function(x) x$fit$d

#' @export 
block_index_list.mfa <- function(x) x$block_indices

#' @export
ncomp.mfa <- function(x) ncomp(x$fit)

#' @export
scores.mfa <- function(x) {
  scores(x$fit)
}

#' @export
loadings.mfa <- function(x) {
  loadings(x$fit)
}


#' @export
partial_scores.mfa <- function(x, table_index=x$ntables) {
  bind <- block_indices(x)
  res <- lapply(table_index, function(i) {
    x$ntables * project(x$fit, get_block(Xr, i), subind=bind[[i]])
  })
  
  names(res) <- paste0("Block_", table_index)
  res
}

project.mfa <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, table_index=NULL) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  
  assert_that(is.matrix(newdata))
  
  if (is.null(table_index)) {
    # new data must have same number of columns as original data
    assert_that(ncol(newdata) == ncol(x$X))
    xnewdat <- x$reprocess(newdata)
    project(x$fit, xnewdat, comp=comp)
  } else if (length(table_index) == 1) {

    xnewdat <- x$reprocess(newdata[[i]], table_index)
    x$ntables * project(x$fit, xnewdat, 
                        comp=comp, subind=x$block_indices[[table_index]])
   
  } else {
    stop("table_index must have length of 1")
  }
  
} 

#' @export
reconstruct.mfa <- function(x, comp=1:ncomp(x)) {
  recon <- reconstruct(x$fit, comp)
}


#' @export
contributions.mfa <- function(x, type=c("column", "row", "table")) {
  contr <- contributions(x$fit)
  type <- match.arg(type)
  
  if (type == "table") {
    out <- do.call(cbind, lapply(1:ncomp(x), function(i) {
      sapply(1:x$ntables, function(j) {
        ind <- x$block_indices[[j]]
        sum(contr[i,ind])
      })
    }))
      
  } else if (type == "row") {
    contributions(x$fit, type="row")
  } else {
    contributions(x$fit, type="column")
  }
  
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



## project from existing table
#' @export
predict.mfa <- function(x, newdata, ncomp=x$ncomp, table_index=1:x$ntables, pre_process=TRUE) {
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
      
      project(x$fit, Xp, ncomp=ncomp) * x$ntables
    }))
  } else if (length(table_index) == 1) {
    ind <- x$block_indices[[table_index]]
    
    Xp <- if (pre_process) {
      x$reprocess(newdata, table_index)
    } else {
      newdata
    }

    fscores <- project(x$fit, Xp, ncomp=ncomp, subind=ind) * x$ntables

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


  
