

#' block_pca
#' 
#' @param X
#' @param groups
#' @param method
#' @param ncomp
#' @param min_k
#' @param max_k
#' @param center
#' @param scale
#' @importFrom furrr future_map
block_pca <- function(X, est_method=c("gcv", "shrink", "fixed", "nneg"), 
                   ncomp=2, min_k=1, max_k=3, 
                   center=TRUE, scale=FALSE, shrink_method="GSURE", ...) {
  
  
  assert_that(inherits(X, "block_matrix"))
  
  est_method <- match.arg(est_method)
  
  ngroups <- nblocks(X)
  
  if (length(ncomp) == 1) {
    ncomp <- rep(ncomp, ngroups)
  } else {
    assertthat::assert_that(length(ncomp) == ngroups)
  }

  fits <- if (est_method == "gcv") {
    fits <- furrr::future_map(1:nblocks(X), function(i) {
      xb <- get_block(X, i)
      est <- fast_estim_ncomp(xb, ncp.min=min_k, ncp.max=max_k)
      pca(xb, ncomp=max(min_k, est$bestcomp), center=center, scale=scale)
    })
  } else if (est_method == "fixed") {
    ## method is fixed
    furrr::future_map(1:nblocks(X), function(i) {
      xb <- get_block(X, i)
      pca(xb, ncomp=ncomp[i], center=center, scale=scale)
    })
  } else {
    furrr::future_map(1:nblocks(X), function(i) {
      xb <- get_block(X, i)
      shrink_pca(xb, center=center, scale=scale, ...)
    })
  }
  
  bm <- block_projector(fits)

}



multiscale_pca <- function(X, cutmat, est_method=c("fixed", "gcv", "shrink"), ncomp=rep(1, length(cuts)+1), 
                           center=TRUE, scale=FALSE, shrink_method="GSURE") {
  
  est_method <- match.arg(est_method)
  preproc <- pre_processor(X, center=center, scale=scale)
  Xp <- pre_process(preproc, X)
 
  
  #offset <- rowMeans(Xp)
  #Xp <- sweep(Xp, 1, offset, "-")

  do_pca <- function(x, level, cind) {
    isplit <- split(1:length(cind), cind)
    lens <- sapply(isplit, length)
    
    ## rescale x to match number of elements in each cluster
    x <- sweep(x, 2, sqrt(lens), "*")
  
    fit <- if (est_method == "fixed") {
      print(ncomp[level])
      pca(x, ncomp=ncomp[level], center=FALSE, scale=FALSE)
    } else if (est_method == "gcv") {
      #est <- fast_estim_ncomp(x)
      est <- FactoMineR::estim_ncp(x)
      bestcomp <- est$ncp
      print(paste("best", bestcomp))
      if (bestcomp == 0) {
        ## 0 comps, return dummy pca
        pseudo_svd(u=matrix(rep(0, nrow(x))), v=matrix(rep(0,ncol(x))), d=0)
      } else {
        pca(x, ncomp=bestcomp, center=FALSE, scale=FALSE)
      }
    } else {
      shrink_pca(x, center=FALSE, scale=FALSE, method=shrink_method)
    }
    
  
    
    out <- matrix(0, length(cind), ncomp(fit))
    
    for (i in 1:length(isplit)) {
      ind <- isplit[[i]]
      v <- fit$v[i,,drop=FALSE]
      vex <- apply(v, 2, rep, length(ind))
      vex <- apply(vex, 2, function(vals) vals/sqrt(length(ind)))
      out[ind,] <- vex
    }
    
    nout <- apply(out, 2, function(vals) vals/ norm(as.matrix(vals), "F"))
    
    ## expanded X
    #Xex <- x[,cind]
  
    #browser()
    #sc <- fit$scores
    #denom <- sqrt(sapply(isplit, length))
    #sweep(sc, 1, denom, "*")
    proj <- bi_projector(preproc=pre_processor(matrix(), center=FALSE, scale=FALSE),
                 ncomp=ncomp(fit),
                 v=nout,
                 u=fit$u,
                 d=fit$d,
                 scores=fit$scores,
                 classes="expanded_pca")
    
    proj
  }
                 
  fits <- list(ncol(cutmat))
  Xresid <- t(Xp)
  
  for (i in 1:ncol(cutmat)) {
    print(i)
    cind <- cutmat[,i]
    Xbar <- t(group_means(cind, Xresid))
    #offset <- rowMeans(Xbar) 
    #Xbar <- sweep(Xbar, 1, offset, "-")
    
    fits[[i]] <- do_pca(Xbar, i, cind)
    
    supp_lds <- project_cols(fits[[i]], t(Xresid))
    recon <- scores(fits[[i]]) %*% t(supp_lds)
    Xresid <- Xresid - t(recon)
    #Xresid <- t(do.call(rbind, lapply(1:nrow(Xbar), function(j) {
    #  xb <- Xbar[j,]
    #  Xresid[,j] - xb[cind]
    #})))
    
    print(sqrt(sum(Xresid^2)))
   
  }
    
  
}


#' hclust_pca
#' 
#' A type of pca that uses a hierarchical clustering to define a set of nested regions for pca compression.
#'   
#' @param X the input matrix
#' @param hclus a hierarchical clustering instance with as many objects as there are rows in \code{X}
#' @param cuts desired number of clusters at each level of the hierarchy (must be increasing)
#' @param k the number of components to estimate at each level
#' @examples 
#' 
#' grid <- expand.grid(1:10, 1:10)
#' X <- matrix(rnorm(100*50), 100, 50)
#' cuts <- c(4, 8, 16)
#' hclus <- hclust(dist(grid))
#' hres1 <- hpca(X, hclus, cuts, est_method="fixed", ncomp=c(4,1,1,1))
#' ncomp(hres1) == (sum(cuts) +4)
#' 
#' hres2 <- hpca(X, hclus, cuts, est_method="shrink")
#' hres3 <- hpca(X, hclus, cuts, est_method="gcv")
#' 
hpca <- function(X, hclus, cuts, est_method=c("fixed", "gcv", "shrink"), ncomp=rep(1, length(cuts)+1), 
                 center=TRUE, scale=FALSE, shrink_method="GSURE") {
  
  est_method <- match.arg(est_method)
  preproc <- pre_processor(X, center=center, scale=scale)
  Xp <- pre_process(preproc, X)
  
  do_pca <- function(x, level) {
    if (est_method == "fixed") {
      pca(x, ncomp=ncomp[level], center=FALSE, scale=FALSE)
    } else if (est_method == "gcv") {
      est <- fast_estim_ncomp(x)
      if (est$bestcomp == 0) {
        ## 0 comps, return dummy pca
        pseudo_svd(u=matrix(rep(0, nrow(x))), v=matrix(rep(0,ncol(x))), d=0)
      } else {
        pca(x, ncomp=nc, center=FALSE, scale=FALSE)
      }
    } else {
      shrink_pca(x, center=FALSE, scale=FALSE, method=shrink_method)
    }
  }
  
  fit0 <- do_pca(Xp, 1)
  
  if (fit0$d[1] == 0) {
    stop("hpca: degenerate solution, first singular value is zero. For fixed solutions, use  'est_method' = 'fixed'")
  }
  
  
  Xresid <- residuals(fit0,xorig=Xp)
  
  fits <- list()
  fits[[1]] <- fit0
  
  for (i in 1:length(cuts)) {
    print(sum(Xresid^2))
    kind <- cutree(hclus, cuts[i])
    
    pres <- split(1:length(kind), kind) %>% furrr::future_map(function(ind) {
      x <- Xresid[ind,]
    
      fit <- do_pca(x, i+1)
      print(ncomp(fit))
      u <- Matrix::sparseMatrix(i=rep(ind, ncol(fit$u)), 
                                j=rep(1:ncol(fit$u), each=length(ind)), 
                                x=as.vector(fit$u),
                                dims=c(nrow(X), ncol(fit$u)))

      pseudo_svd(u,fit$v,fit$d)
      
    })
    
    ds <- map_dbl(pres, ~ .$d[1])
    
    if (all(ds == 0)) {
      break
    } else if (sum(ds) == 1) {
      ind <- which(ds == 1)
      fits[[i+1]] <- pres[[ind]]
    } else {
      pres <- pres[ds > 0]
      fits[[i+1]] <- do.call(partial(combine, pres[[1]]), pres[2:length(pres)])
    }
    
    Xresid <- residuals(fits[[i+1]], ncomp=ncomp(fits[[i+1]]), xorig=Xresid)
  }
  
  final_fit <- if (length(fits) > 1) {
    do.call(partial(combine, fits[[1]]), fits[2:length(fits)])
  } else {
    fits[[1]]
  }
  
}

