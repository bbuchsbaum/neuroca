


normalization_factors <- function(block_mat, type=c("MFA", "RV", "RV-MFA", "None", "Frob", "Dual-RV")) {
  type <- match.arg(type)
  assert_that(inherits(block_mat, "block_matrix"))
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
  } else if (type == "RV" && nblocks(block_mat) > 2) {
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    abs(svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
  } else if (type == "RV-MFA") {
    alpha1 <- unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    smat[which(smat < 0)] <- 0
    diag(smat) <- 1
    alpha2 <- abs(svd_wrapper(smat, ncomp=1, method="propack")$u[,1])
    alpha1*alpha2
  } else if (type == "Dual-RV") {
    bl <- block_lengths(block_mat)
    assertthat::assert_that(all(bl[1] == bl), msg="Dual-RV requires that all blocks have same variables.")
    smat1 <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    smat2 <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(t(x1),t(x2)))
    smat <- (smat1 + smat2)/2
    diag(smat) <- 1
    smat[which(smat < 0)] <- 0
    abs(svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
  }else if (type == "Frob") {
    unlist(lapply(as.list(block_mat), function(X) sum(X^2)))
  } else {
    rep(1, nblocks(block_mat))
  }
}



#' multiple factor analysis
#' 
#' mfa
#' 
#' @param X a \code{block_matrix} object
#' @param ncomp the number of components to estimate
#' @param preproc a preprocessing pipeline, default is `center()`
#' @param normalization the normalization method: MFA, RV, RV-MFA, or None (see details).
#' @param A custom weight matrix for the columns
#' @param M custom weight matrix for the rows
#' @export
#' @examples 
#' 
#' X <- block_matrix(replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE))
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' p <- partial_scores(res, 1)
#' stopifnot(ncol(scores(res)) == 3)
#' 
#' labs <- letters[1:10]
#' cfier <- classifier(res, labels=labs, ncomp=3)
#' pred <- predict(cfier, X[1:2,])
#' cfier2 <- classifier(res, labels=labs, ncomp=3, colind=res$block_indices[[2]])
#' pred2 <- predict(cfier2, X[1:2,res$block_indices[[2]]])
mfa <- function(X, ncomp=2, preproc=center(), 
                normalization=c("MFA", "RV", "None", "RV-MFA", "Dual-RV", "Frob", "custom"), M=NULL, A=NULL, ...) {

  
  assertthat::assert_that(inherits(X, "block_matrix"), msg="X must be a 'block_matrix'")
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom") {
    assert_that(!is.null(A))
  }

  ## pre-process the variables.
  procres <- prep(preproc, X)
  Xp <- procres$Xp

  ## normalize the matrices 
  
  if (normalization != "custom") {
    alpha <- normalization_factors(Xp, type=normalization)
    A <- rep(alpha, block_lengths(X))
  } else {
    alpha <- rep(1, nblocks(X))
  }
  
  bind <- block_index_list(X)
  
  fit <- genpca(unclass(Xp), 
                preproc=pass(),
                    A=A, 
                    M=M,
                    ncomp=ncomp,
                    ...)
  

  result <- list(
    X=X,
    preproc=procres,
    ntables=nblocks(X),
    fit=fit,
    ncomp=fit$ncomp,
    block_indices=bind,
    alpha=alpha,
    normalization=normalization,
    table_names=names(X),
    nvars=ncol(X),
    ntables=length(block_lengths(X)),
    A=A,
    M=M
  )
  
  class(result) <- c("mfa", "multiblock", "bi_projector", "projector", "list")
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
vfa.mfa <- function(x, type = c("component", "block"), ncomp=x$ncomp) {
  recon <- reconstruct(x, comp=1)
  block_inds <- block_index_list(x$X)
  
  do.call(rbind, lapply(1:ncomp, function(i) {
    bl <- to_block_matrix(recon, block_lengths(x$X))
    do.call(rbind, lapply(1:length(block_inds), function(j) {
      Xj <- get_block(x$X, j)
      Rj <- get_block(recon, j)
      sserr <- sum((Xj - Rj)^2)
      v <- 1 - (sserr/(Xj^2))
      data.frame(component=i, block=j, vaf=v)
    }))
  }))
    
}

projection_fun.mfa <- function(x, colind=NULL, block=NULL) {
  if (!is.null(colind) && !is.null(block)) {
    stop("projection_fun.mfa: cannot provide both 'colind' or 'block' arguments")
  }
  
  if (!is.null(block)) {
    colind <- res$block_indices[[block]]
  }

  if (is.null(colind)) {
    stop("projection_fun.mfa: muust provide either 'colind' or 'block' argument")
  }
  
  f <- function(newdata) {
    project(x, newdata, colind=colind)
  }
}
  
  
#' @export
partial_scores.mfa <- function(x, block_index=x$ntables) {
 
  bind <- block_index_list(x)
  
  res <- lapply(block_index, function(i) {
    ## FIXME
    x$ntables * project(x$fit, get_block(x$X, i), colind=bind[[i]], pre_process=FALSE)
  })
  
  names(res) <- paste0("Block_", block_index)
  res
}

#' @export
block_project.mfa <- function(x, newdata, block=1, comp=1:ncomp(x)) {
  assert_that(length(block) == 1, msg="block must have length of 1")
  subind <- block_index_list(x$X)[[block]]
  newdat <- reprocess(x, newdata, colind=subind)
  project(x$fit, unclass(newdata), comp=comp, colind=subind)
}

#' @export
project.mfa <- function(x, newdata, comp=1:ncomp(x), colind=NULL) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  
  if (is.null(colind)) {
    assert_that(ncol(newdata) == ncol(x$X), msg=paste("ncol(newdata) =  ", ncol(newdata), " ncol(x$X) = ", ncol(x$X)))
    project(x$fit, unclass(reprocess(x, newdata)), comp=comp)
  } else {
    project(x$fit, unclass(reprocess(x, newdata,colind=colind)), comp=comp, colind=colind)
  }
} 


#' @export
project_cols.mfa <- function(x, newdata=NULL, comp=1:x$ncomp) {
  project_cols(x$fit,newdata,comp=comp)
}


#' @export
reconstruct.mfa <- function(x, newdata=NULL, comp=1:ncomp(x)) {
  recon <- reconstruct(x$fit, newdata, comp)
  x$preproc$reverse_transform(recon)
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
  
  ## TODO this currently rotates the scores. 
  ## could simply rotate the original variables
  res <- lapply(1:length(x$Xlist), function(i) {
    xcur <- get_block(x$Xp[[i]])
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
## not clear that we need this
#' @export
predict.mfa <- function(x, newdata, ncomp=x$ncomp, block_index=1:x$ntables, pre_process=TRUE) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  assert_that(is.matrix(newdata))
  assert_that(length(block_index) == 1 || length(block_index) == x$ntables)
  
  fscores <- if (length(block_index) == x$ntables) {
    assert_that(ncol(newdata) == ncol(x$X))
    Reduce("+", lapply(block_index, function(i) {
      ind <- x$block_indices[[i]]
      
      Xp <- if (pre_process) {
        reprocess(x, newdata[, ind,drop=FALSE], colind=ind)
      } else {
        newdata[, ind]
      }
      
      project(x$fit, Xp, ncomp=ncomp) * x$ntables
    }))
  } else if (length(block_index) == 1) {
    ind <- x$block_indices[[block_index]]
    
    Xp <- if (pre_process) {
      reprocess(x, newdata, colind=ind)
    } else {
      newdata
    }

    fscores <- project(x$fit, Xp, comp=1:ncomp, colind=ind) * x$ntables

  }
  
}



#' @export
impute_mfa <- function(X, ncomp=min(dim(X)), preproc=center(), 
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
    xtmp <- X
    xtmp[missing_ind] <- 0
    xtmp[missing_ind] <- xtmp[missing_ind] + Xmiss_init[missing_ind]
    Ximp <- xtmp
  }
  
  old <- Inf
  criterion <- Inf
  iter <- 1
  
  while(iter < iter_max && criterion > threshold) {
    mres <- mfa(Ximp, ncomp=ncomp, preproc=preproc, normalization=normalization)
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


  
