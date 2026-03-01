
#' procrustes rotation
proc_rot <- function(X, Y) {
  XY <- crossprod(X, Y)
  sol <- svd(XY)
  A <- sol$v %*% t(sol$u)
  A
}


#' compute standard deviations for a set of bootstrap results
boot_sd <- function(bootlist) {
  boot.mean <- Reduce("+", bootlist)/length(bootlist)
  boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
}

#' compute bootstrap ratios for a list of matrices
boot_ratio <- function(bootlist) {
  boot.mean <- Reduce("+", bootlist)/length(bootlist)
  boot.sd <- boot_sd(bootlist)
  boot.mean/boot.sd	
}

bootci <- function(bootlist, ci=c(.025,.975)) {
  #browser()
  nboot <- length(bootlist)
  nr <- nrow(bootlist[[1]])
  ranks <- ci * nboot
  ret <- lapply(1:nr, function(i) {
    v <- sapply(bootlist, function(x) x[i,])
    apply(v,1, function(vals) {
      svals <- sort(vals)
      c(svals[ranks[1]], svals[ranks[2]])
    })
  })
  
  nc <- ncol(bootlist[[1]])
  
  ret2 <- lapply(1:nc, function(i) {
    m <- do.call(rbind, lapply(ret, function(z) z[,i]))
    colnames(m) <- c("ci_min", "ci_max")
    m
  })
  
  names(ret2) <- paste0("Comp", 1:nc)
  ret2
  
}

# given a matrix and a vector, resample with replacement.
resampleXY <- function(X, Y) {
  ysplit <- split(1:length(Y), Y)
  yidx <- sort(unlist(lapply(ysplit, function(ids) {
      sample(ids, replace=TRUE)
  })))
  
  list(X=X[yidx,], Y=Y[yidx])
}


#' @export
resample.mubada <- function(x) {
  ## check for multiple repetitions per Y levels per block
  Xresam <- x$Xlist
  Yresam <- x$Y
  
  ret <- lapply(1:length(Yresam), function(i) {
    resampleXY(Xresam[[i]], Yresam[[i]])
  })
  
  Xret <- lapply(ret, "[[", "X")
  Yret <- lapply(ret, "[[", "Y")
  list(X=Xret, Y=Yret)
}


#' @export
resample.bada <- function(x) {
  ## split all indices by category
  Y_by_idx <- split(1:length(x$Y), x$Y)
  resam <- unlist(lapply(Y_by_idx, function(i) sample(i, replace=TRUE)))
  list(Y=x$Y[resam], X=x$X[resam,])
}

#' @rdname bootstrap
#' @examples 
#' 
#' ## mubada bootstrap analysis
#' Xl <- lapply(1:5, function(i) matrix(rnorm(100*20), 100, 20))
#' Yl <- lapply(1:5, function(i) factor(rep(letters[1:5], length.out=100)))
#' 
#' mb <- mubada(Yl, Xl)
#' @export
bootstrap.mubada <- function(x, nboot=100, ncomp=x$ncomp,type=c("projection", "rotated", "unrotated")) {
  type <- match.arg(type)
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  
  do_boot_projection <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample(x)
    
    ## calculate barycenters
    Xb <- do.call(cbind, lapply(1:length(resam$X), function(i) {
      reprocess(x, group_means(resam$Y[[i]], resam$X[[i]]), block_index=i)
    }))
    
    #browser()
    #browser()
    ## project rows and columns of barycenters
    br <- project(x, Xb)
    bc <- project_cols(x, Xb)
    
    list(boot4R=br,
         boot4C=bc)
  }
  
  do_boot_svd <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample(x)
    xboot <- mubada(resam$Y, resam$X, ncomp=x$ncomp, center=x$center, scale=x$scale, normalization=x$normalization)
    
    if (type == "rotated") {
      A <- proc_rot(scores(x), scores(xboot))
      boot_scores <- scores(xboot) %*% A
      boot_loadings <- loadings(xboot) %*% A
    } else {
      sign_flip <- sign(diag(t(scores(x)) %*% scores(xboot)))
      boot_scores <- t(t(scores(xboot)) * sign_flip)
      boot_loadings <- t(t(loadings(xboot)) * sign_flip)
    }
    
    list(boot4R=boot_scores,
         boot4C=boot_loadings)
  }
  
  ret <- if (type == "projection") {
    boots <- replicate(nboot, do_boot_projection(), simplify = FALSE)
    list(zboot_scores=boot_ratio(lapply(boots, "[[", 1)),
         zboot_loadings=boot_ratio(lapply(boots, "[[", 2)))
    
  } else {
    boots <- replicate(nboot, do_boot_svd(), simplify=FALSE)
    list(zboot_scores=boot_ratio(lapply(boots, "[[", 1)),
         zboot_loadings=boot_ratio(lapply(boots, "[[", 2)),
         nboot=nboot)
  }
  
  
  class(ret) <- c("bootstrap_result", "list")
  ret
  
}


#' @export
bootstrap.bada <- function(x, nboot=1000, ncomp=x$ncomp, 
                           type=c("projection", "permutation", "rotated", "unrotated"),
                           ci=c(.025, .975)) {
  
  type <- match.arg(type)
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  do_perm <- function() {
    split_idx <- split(1:length(x$Y), x$S)
    Yperm <- unlist(lapply(split_idx, function(idx) {
      sample(x$Y[idx])
    }))

   
    bperm <- bada(Yperm, x$X, S=x$S, preproc=x$preproc)
    
    list(boot4R=scores(bperm)[,1:ncomp],
         boot4C=loadings(bperm)[,1:ncomp])
  }
  
 
  do_boot_projection <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample(x)
    
    ## calculate barycenters of pre-processed data
    Xb <- group_means(resam$Y, reprocess.bada(x, resam$X))
    
    ## project rows and columns of barycenters
    bc <- project_cols(x, Xb)
    br <- project(x, Xb)
    
    list(boot4R=br,
         boot4C=bc)
  }
  
  do_boot_svd <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample.bada(x)
    xboot <- bada(resam$Y, resam$X, x$S, preproc=x$preproc)
    
    if (type == "rotated") {
      A <- proc_rot(scores(x), scores(xboot))
      boot_scores <- scores(xboot) %*% A
      boot_loadings <- loadings(xboot) %*% A
    } else {
      sign_flip <- sign(diag(t(scores(x)) %*% scores(xboot)))
      boot_scores <- t(t(scores(xboot)) * sign_flip)
      boot_loadings <- t(t(loadings(xboot)) * sign_flip)
    }
    
    list(boot4R=boot_scores,
         boot4C=boot_loadings)
  }
  
  ret <- if (type == "projection") {
    boots <- replicate(nboot, do_boot_projection(), simplify = FALSE)
    list(zboot_scores=boot_ratio(lapply(boots, "[[", 1)),
         zboot_loadings=boot_ratio(lapply(boots, "[[", 2)))
         #bootci_loadings=bootci(lapply(boots, "[[", 2)),
         #bootci_scores=bootci(lapply(boots, "[[", 1)))
         
    
  } else if (type == "permutation") {
    boots <- replicate(nboot, do_perm(), simplify = FALSE)
    list(bootci_scores=bootci(lapply(boots, "[[", 1)),
         bootci_loadings=bootci(lapply(boots, "[[", 2)))
    
  } else {
    boots <- replicate(nboot, do_boot_svd(), simplify=FALSE)
    list(zboot_scores=boot_ratio(lapply(boots, "[[", 1)),
         zboot_loadings=boot_ratio(lapply(boots, "[[", 2)),
         nboot=nboot)
  }

 
  class(ret) <- c("bootstrap_result", "list")
  ret
  
}




#' @export
#' @examples 
#' 
#' X <- matrix(rnorm(10*100), 10, 100)
#' x = pca(X, ncomp=9)
bootstrap.pca <- function(x, nboot=100, k=x$ncomp) {
  DUt <- t(scores(x))
  n <- dim(DUt)[2]
  
  gen <- function() {
    sidx <- sample(1:ncol(DUt), replace=TRUE)
    list(DUt=DUt[,sidx,drop=FALSE], idx=sidx)
  }
  
  
  res <- boot_svd(nboot=nboot, k=k, as.matrix(loadings(x)), gen)

  zboot <- do.call(cbind, lapply(1:k, function(ki) {
    res$EVs[[ki]]/res$sdVs[[ki]]
  }))
  
  zscores <- do.call(cbind, lapply(1:k, function(ki) {
    res$EScores[[ki]]/res$sdScores[[ki]]
  }))
  
  ret <- list(zboot_loadings=zboot, zboot_scores=zscores, nboot=nboot, k=k)
  class(ret) <- c("bootstrap_result", "list")
  ret
 
}




#' @keywords internal
svd_dutp <- function(DUtP,k) {
  n<-dim(DUtP)[2]
  svdDUtP <- svd(DUtP)
  sb <- svdDUtP$d
  
  sign_flip <- sign(diag(svdDUtP$u))
  
  sign_flip[sign_flip==0]<-1 
  sign_flip <- sign_flip[1:k] 
  
  
  Ab <- svdDUtP$u[1:min(dim(DUtP)), 1:k]
  Ub <- svdDUtP$v[1:n, 1:k] 
  
  Ab <- sweep(Ab,2,sign_flip, "*")
  Ub <- sweep(Ub,2,sign_flip, "*")
  
  list(d=sb, Ab=Ab, Ub=Ub)
}

#' @keywords internal
boot_sum <- function(res,k, v) {
  
  ## each of k elements has nboot rows and n columns
  AsByK <- lapply(1:k, function(ki) {
    do.call(rbind, lapply(res, function(a) {
      a$svdfit$Ab[,ki]
    }))
  })
  
  ScoresByK <- lapply(1:k, function(ki) {
    do.call(rbind, lapply(res, function(a) {
      u <- a$svdfit$Ub[,ki] * a$svdfit$d[ki]
      u2 <- rep(NA, length(u))
      u2[a$idx] <- u
      u2
    }))
  })
  
  ## mean scores
  EScores <- lapply(ScoresByK, function(s) {
    apply(s, 2, mean, na.rm=TRUE)
  })
  
  ## sd of scores
  sdScores <- lapply(ScoresByK, function(s) {
    apply(s, 2, sd, na.rm=TRUE)
  })
  
  
  ## produces 5 mean vectors, 1 per component
  EAs <- lapply(AsByK, colMeans) 
  
  ## 5 loadings components
  EVs <- lapply(EAs, function(EA) v %*% matrix(EA,ncol=1))
  
  
  ## k nXn covariance matrices
  varAs <- lapply(AsByK,var) #indexed by k
  
  varVs <- lapply(1:length(AsByK), function(ki) {
    rowSums((v %*% varAs[[ki]]) * v)
  })

  sdVs <- lapply(varVs,sqrt)
  list(res=res, EAs=EAs, EVs=EVs, varAs=varAs, sdVs=sdVs, EScores=EScores, sdScores=sdScores)
  
}

#' @keywords internal
boot_svd <- function(nboot, k, v, gen_DUtP) {
  
  ## Generate nboot resamples of scores
  res <- lapply(1:nboot, function(i) {
    sam <- gen_DUtP()
    DUtP <- sam$DUt
    #DUtP <- if(x$center) t(scale(t(DUt[,sidx]),center=TRUE,scale=FALSE)) else DUt[,sidx]
    list(svdfit=svd_dutp(DUtP,k), idx=sam$idx)
  })
  
  
  boot_sum(res,k, v)
  
}


# bootstrap2.bada <- function(x, nboot=1000, k=x$ncomp) {
#   ## get the latent projections
#   proj <- project(x, comp=1:k)
#   
#   ## split all indices by category
#   Y_by_idx <- split(1:length(x$Y), x$Y)
#   
#   gen <- function() {
#     ## sample with replacement within each category
#     sam_idx <- sort(unlist(lapply(Y_by_idx, function(i) sample(i, replace=TRUE))))
#     
#     ## recompute category means
#     t(group_means(x$Y[sam_idx], proj[sam_idx,]))
#   }
#   
#   res <- boot_svd(nboot=nboot, k=k, loadings(x), gen)
#   
#   zboot <- do.call(cbind, lapply(1:k, function(ki) {
#     res$EVs[[ki]]/res$sdVs[[ki]]
#   }))
#   
# }
# 

