


normalization_factors <- function(block_mat, type=c("MFA", "RV", "RV-MFA", "None")) {
  type <- match.arg(type)
  
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
  } else if (type == "RV" && nblocks(block_mat) > 2) {
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    wts <- abs(svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
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



#' multiple factor analysis
#' 
#' mfa
#' 
#' @param X a \code{block_matrix} or a \code{block_projection_matrix} object
#' @param ncomp the number of components to estimate
#' @param center whether to mean center the variables
#' @param scale whether to standardize the variables
#' @param normalization the normalization method: MFA, RV, RV-MFA, or None (see details).
#' @export
#' @examples 
#' 
#' X <- block_matrix(replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE))
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' p <- partial_scores(res, 1)
#' stopifnot(ncol(scores(res)) == 3)
mfa <- function(X, ncomp=2, center=TRUE, scale=FALSE, 
                normalization=c("MFA", "RV", "None", "RV-MFA", "custom"), A=NULL) {
  
  assertthat::assert_that(inherits(X, "block_matrix"), msg="X must be a 'block_matrix'")
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom") {
    assert_that(!is.null(A))
  }

  ## pre-process the variables.
  preproc <- pre_processor(X, center=center,scale=scale)
  Xr <- pre_process(preproc,X)

  ## normalize the matrices 
  
  if (normalization != "custom") {
    alpha <- normalization_factors(Xr, type=normalization)
    A <- rep(alpha, block_lengths(X))
  } else {
    alpha <- rep(1, nblocks(X))
  }
  
  bind <- block_index_list(X)
  
  fit <- genpca(unclass(Xr), 
                    A=A, 
                    ncomp=ncomp, 
                    center=FALSE, 
                    scale=FALSE)
  
  
  result <- list(
    X=Xr,
    preproc=preproc,
    ntables=nblocks(X),
    fit=fit,
    ncomp=fit$ncomp,
    center=center,
    scale=scale,
    block_indices=bind,
    alpha=alpha,
    normalization=normalization,
    table_names=names(X),
    A=A
  )
  
  class(result) <- c("mfa", "multiblock", "projector", "list")
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
partial_scores.mfa <- function(x, block_index=x$ntables) {
 
  bind <- block_index_list(x)
  res <- lapply(block_index, function(i) {
    ## FIXME
    x$ntables * project(x$fit, get_block(x$Xr, i), subind=bind[[i]], pre_process=FALSE)
  })
  
  names(res) <- paste0("Block_", block_index)
  res
}

#' @export
project.mfa <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, block_index=NULL) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  
  if (!is.null(block_index)) {
    assert_that(length(block_index) == 1, msg="block_index must have length of 1")
    subind <- block_index_list(x$X)[[block_index]]
    newdat <- reprocess(x, newdata, block_index)
    x$ntables * project(x$fit, unclass(newdat), comp=comp, subind=subind)
  } else {
    # new data must have same number of columns as original data
    assert_that(ncol(newdata) == ncol(x$X), msg=paste("ncol(newdata) =  ", ncol(newdata), " ncol(x$X) = ", ncol(x$X)))
    project(x$fit, unclass(reprocess(x, newdata)), comp=comp)
  } 
} 

#' @export
project_cols.mfa <- function(x, newdata=NULL, comp=1:x$ncomp) {
  project_cols(x$fit,newdata,comp=comp)
}


#' @export
reconstruct.mfa <- function(x, newdata=NULL, comp=1:ncomp(x)) {
  recon <- reconstruct(x$fit, newdata, comp)
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
    xp <- project(x, xcur, comp=1:ncomp, block_index=i)
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
predict.mfa <- function(x, newdata, ncomp=x$ncomp, block_index=1:x$ntables, pre_process=TRUE) {
  assert_that(is.matrix(newdata))
  assert_that(length(block_index) == 1 || length(block_index) == x$ntables)
  
  fscores <- if (length(block_index) == x$ntables) {
    assert_that(ncol(newdata) == ncol(x$X))
    Reduce("+", lapply(block_index, function(i) {
      ind <- x$block_indices[[i]]
      
      Xp <- if (pre_process) {
        reprocess(x, newdata[, ind,drop=FALSE], block_index=i)
      } else {
        newdata[, ind]
      }
      
      project(x$fit, Xp, ncomp=ncomp) * x$ntables
    }))
  } else if (length(block_index) == 1) {
    ind <- x$block_indices[[block_index]]
    
    Xp <- if (pre_process) {
      reprocess(x, newdata, block_index=block_index)
    } else {
      newdata
    }

    fscores <- project(x$fit, Xp, comp=1:ncomp, subind=ind) * x$ntables

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
  
  while(iter < iter_max && criterion > threshold) {
    mres <- mfa(Ximp, ncomp=ncomp, center=center, scale=scale, normalization=normalization)
    recon <- reconstruct(mres)
    
    Xnew <- X
    Xnew[missing_ind] <- recon[missing_ind]
    
    obj <- sum((Xnew[-missing_1d] - recon[-missing_1d])^2)
   
    Ximp <- Xnew
    criterion <- abs(1 - obj/old)
    message("iter: ", obj)
    print(criterion)
    old <- obj
    iter <- iter+1
  }
  
  Ximp
}


  
