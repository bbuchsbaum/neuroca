
centering_matrix <- function(vals) {
  p <- cbind(diag(length(vals)), -vals)
  #rbind(p, c(rep(0, length(vals)),1))
}

#' 
#' 
#' @example 
#' 
#' X1 <- matrix(rnorm(20*10), 20, 10)
#' X2 <- matrix(rnorm(20*20), 20, 20)
#' X3 <- matrix(rnorm(20*30), 20, 30)
#' 
#' pc1 <- pca(X1, ncomp=10, preproc=center())
#' pc2 <- pca(X2, ncomp=10, preproc=center())
#' pc3 <- pca(X3, ncomp=10, preproc=center())
#' 
#' Xg <- cbind(X1,X2,X3)
#' 
#' fits <- list(pc1,pc2,pc3)
#' pfit <- metapca(fits, ncomp=15, scale_inner=TRUE)
#' ncol(scores(pfit)) == 10
#' @export
metapca <- function(fits, ncomp=2, scale_inner=FALSE) {
  assert_that(all(sapply(fits,function(f) inherits(f, "projector"))))
  
  X <- block_matrix(lapply(fits, function(x) project(x)))
  
  pres <- if (scale_inner) {
    wts <- (1/sqrt(matrixStats::colSds(X)))
    genpca(unclass(X), A=wts, ncomp=ncomp, preproc=preproc)
  } else {
    pca(X, ncomp=ncomp, preproc=preproc)
  }
  
  
  nvars <- sapply(fits, ncol)
  tot <- sum(nvars)
  offsets <- cumsum(c(1, sapply(fits, ncol)))
  outer_block_indices <- lapply(1:length(nvars), function(i) {
    seq(offsets[i], offsets[i] + nvars[i]-1)
  })
  
  ret <- bi_projector(
    preproc=pass(),
    ncomp=length(pres$d),
    v=pres$v, 
    u=pres$u, 
    d=pres$d, 
    scores=pres$scores,
    metafit=pres,
    fits=fits,
    outer_block_indices=outer_block_indices,
    inner_block_indices=block_index_list(X),
    classes=c("metapca", "pca"))
}



#' @export
block_project.metapca <- function(x, newdata, block=1, comp=1:ncomp(x)) {
  project(x, newdata, comp, colind=x$outer_block_indices[[block]])
}


#' @export
loadings.metapca <- function(x) {
  do.call(rbind, lapply(1:length(fits), function(i) {
    v <- fits[[i]]$v %*% loadings(x$metafit)[x$inner_block_indices[[i]],]
  }))
  
}
  
#' @export
project.metapca <- function(x, newdata, comp=1:ncomp(x), colind=NULL) {
  if (missing(newdata)) {
    scores(x)[,comp,drop=FALSE]
  } else if (!is.null(colind)) {
     
      x0 <- do.call(cbind, lapply(1:length(x$outer_block_indices), function(i) {
        ind <- x$outer_block_indices[[i]]
        keep <- colind %in% ind
        
        if (sum(keep) > 0) {
          idx <- which(colind %in% ind)
          nd <- newdata[,idx,drop=FALSE]
          project(x$fits[[i]], nd,colind=which(keep))
        } else {
          matrix(0, nrow(scores(x$fits[[i]])), x$fits[[i]]$ncomp)
        }
      }))
      project(x$metafit,x0, comp)
  } else {
    assert_that(ncol(newdata) == sum(sapply(x$fits, ncol))) 
    x0 <- do.call(cbind, lapply(1:length(x$outer_block_indices), function(i) {
      ind <- x$outer_block_indices[[i]]
      nd <- newdata[,ind,drop=FALSE]
      project(x$fits[[i]], nd)
    }))
    project(x$metafit,x0, comp)
    
  }
}

#' @export
supplementary_loadings.metapca <- function(x, newdata, block_ind=1, ncomp=x$ncomp) {
  #assert_that(length(block_ind) == 1, msg="can only extract supplementary loadings for one block at a time.")
  if (length(block_ind) == 1) {

    newdata <- reprocess(x$fits[[block_ind]], newdata)
   
    Qsup <- t(newdata) %*% (x$u[,1:ncomp,drop=FALSE]) %*% diag(x=1/x$d[1:ncomp], 
                                                              nrow=length(x$d), ncol=length(x$d))
    t(Qsup)
  } else {
    assert_that(inherits(newdata, "block_matrix"))
    assert_that(nblocks(newdata) == length(block_ind))
    
    ret <- lapply(block_ind, function(bind) {
      nd <- get_block(newdata, bind)
      newdata <- reprocess(x$fits[[bind]], nd)
      #newdata <- reprocess(x$metafit, newdata, colind=x$inner_block_indices[[bind]])
      Qsup <- t(newdata) %*% (x$u[,1:ncomp,drop=FALSE]) %*% diag(x=1/x$d[1:ncomp], 
                                                                 nrow=length(x$d), ncol=length(x$d))
      as.matrix(Qsup)
    })
    
    do.call(rbind, ret)
    
  }
}

#' @export
reconstruct.metapca <- function(x, newdata=NULL, comp=1:length(x$d), 
                                block_ind=seq(1,length(x$fits)), 
                                reverse_pre_process=TRUE) {
  
  recon <- reconstruct(x$metafit, comp=comp, reverse_pre_process=TRUE)
  
  ## reconstruct here means reconstruct the original data
  ret <- lapply(block_ind, function(i) {
    ## reconstruct the inner matrix
    ## recon <- reconstruct(x$metafit, newdata, comp=comp, colind=x$block_indices[[i]], reverse_pre_process=TRUE)
    ## reconstruct the outer matrix
    as.matrix(reconstruct(x$fits[[i]], newdata=recon[,x$inner_block_indices[[i]]], 
                                       reverse_pre_process=reverse_pre_process))
  })
  
  block_matrix(ret)
  
}