

#' block_pca
#' 
#' @param X the data matrix and instance of class \code{block_matrix}
#' @param est_method the pca method to apply to each block
#' @param ncomp the number of components to extract from each block (can vary by block or be fixed)
#' @param min_k the minimum number of components in each block
#' @param max_k the maximum number of components in each block
#' @param center whether to center columns
#' @param scale whether to scale columns
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

    fits <- furrr::future_map(seq(1,nblocks(X)), function(i) {
      xb <- get_block(X, i)
      est <- fast_estim_ncomp(xb, ncp.min=min_k, ncp.max=max_k)
      pca(xb, ncomp=max(min_k, est$bestcomp), center=center, scale=scale)
    })
  } else if (est_method == "fixed") {
    ## method is fixed
    furrr::future_map(seq(1, nblocks(X)), function(i) {
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
#' @param X the data matrix
#' @param hclus a hierarchical clustering instance with as many objects as there are rows in \code{X}
#' @param cuts desired number of clusters at each level of the hierarchy (must be increasing)
#' @param ncomp the number of components to estimate at each level
#' @param center
#' @param scale
#' @param shrink_method
#' @examples 
#' 
#' grid <- expand.grid(1:10, 1:10)
#' 
#' ## 100 images each with 50 features
#' X <- matrix(rnorm(100*50), 100, 50)
#' cuts <- c(4, 8, 16)
#' hclus <- dclust::dclust(grid, nstart=10)
#' hres1 <- hpca(X, hclus, cuts, est_method="fixed", ncomp=c(4,1,1,1))
#' ncomp(hres1) == (sum(cuts) +4)
#' 
#' hres2 <- hpca(X, hclus, cuts, est_method="shrink")
#' hres3 <- hpca(X, hclus, cuts, est_method="gcv")
#' ## start bottom up?
#' @importFrom dendextend cutree
#' @import dclust
#' 
#' 
# res <- neuroca:::hpca(mat, hclus, cuts=c(2,4,8,16,32), est_method="smooth",
#                comp_method="log",
#                cds=cds, spat_smooth=c(16,8,4,2,2))
hpca <- function(X, hclus, cuts, est_method=c("fixed","shrink", "smooth"), 
                 comp_method=c("fixed", "unit", "log"),
                 ncomp=rep(1, length(cuts)),
                 spat_smooth=rep(0, length(cuts)+1),
                 cds=NULL,
                 preproc=center(), 
                 shrink_method="GSURE") {
  
  est_method <- match.arg(est_method)
  comp_method <- match.arg(comp_method)
  #preproc <- pre_processor(X, center=center, scale=scale)
  #Xp <- pre_process(preproc, X)
  
  if (comp_method == "fixed" && ncomp != length(cuts)) {
    stop("`ncomp` must have same length as `cuts`")
  }
  
  if (est_method == "smooth") {
    assert_that(!is.null(cds) && nrow(cds) == nrow(X))
  }
  
  get_ncomp <- function(level, ind) {
    if (comp_method == "log") {
      max(1, ceiling(log(length(ind))))
    } else if (comp_method == "unit") {
      1
    } else if (comp_method=="fixed") {
      ncomp[level]
    }
  }
  
  do_pca <- function(x, level, preproc, ind) {
    if (est_method == "fixed") {
      pca(x, ncomp=get_ncomp(level, ind), preproc=preproc)
    } else if (est_method == "smooth") {
      message("smooth is", spat_smooth[level])
      coord <- cds[ind,]
      S <- spatial_constraints(cds[ind,], sigma_within=spat_smooth[level],
                               nnk_within=spat_smooth[level] * 10)
      #S <- neighborweights::make_doubly_stochastic(S)
      genpca(x, M=S, ncomp=get_ncomp(level, ind), preproc=preproc)
    } else {
      shrink_pca(x, preproc=preproc, method=shrink_method)
    }
  }
  
  fit0 <- do_pca(X, 1, preproc, 1:nrow(X))
  
  if (fit0$d[1] == 0) {
    stop("hpca: degenerate solution, first singular value is zero. For fixed solutions, use  'est_method' = 'fixed'")
  }
  
  
  Xresid0 <- residuals(fit0, ncomp=fit0$ncomp, xorig=X)
  Xresid <- Xresid0
  
  fits <- vector(length(cuts), mode="list")
  fits[[1]] <- fit0$u
  
  for (i in 1:length(cuts)) {
    #print(i)
    print(sum(Xresid^2))
    kind <- dendextend::cutree(hclus, cuts[i])
    
    pres <- split(1:length(kind), kind) %>% furrr::future_map(function(ind) {
      #print(ind)
      x <- Xresid[ind,]
    
      fit <- do_pca(x, i+1, preproc=center(), ind)
      #print(ncomp(fit))
      
      u <- cbind(rep(1, nrow(fit$u)), fit$u)
      u <- Matrix::sparseMatrix(i=rep(ind, ncol(u)), 
                                j=rep(1:ncol(u), each=length(ind)), 
                                x=as.vector(u),
                                dims=c(nrow(X), ncol(u)))
      
      
      
      #psvd <- pseudo_svd(fit$u,fit$v,fit$d, pre_processor=fit$preproc)
      #residuals(psvd)
      #print(paste("ncomp = ", fit$ncomp, " length(d) = ", length(fit$d)))
      
      list(psvd=fit, u=u, resid=residuals(fit, ncomp=length(fit$d), xorig=x), 
           ind=ind, comp=length(fit$d))
    })
    
  
    
    ## fit u to Xresid -- take resid
    ## regress(u, Xresid, intercept=FALSE)
    ds <- map_dbl(pres, ~ .$psvd$d[1])
    
    if (all(ds == 0)) {
      break
    } else if (sum(ds) == 1) {
      ind <- which(ds == 1)
      pfits <- lapply(pres, "[[", "psvd")
      fits[[i+1]] <- pfits[[ind]]$u
    } else {
      pres <- pres[ds > 0]
      pfits <- lapply(pres, "[[", "psvd")
      fits[[i+1]] <- do.call(cbind, lapply(pres, "[[", "u"))
    }
    
    Xout <- matrix(0, nrow(X), ncol(X))
    
    for (pfit in pres) {
      Xout[pfit$ind,] <- as.matrix(pfit$resid)
    }
    Xresid <- Xout
    
  }
  
  u_final <- do.call(cbind, fits)
  reg <- regress(u_final, X, method="ridge", intercept=TRUE)
  class(reg) <- c("hpca", class(reg))
  reg$levels <- length(cuts)
  reg$hclus <- hclus
  reg$cuts <- cuts
  reg$spat_smooth <- spat_smooth
  reg
}

