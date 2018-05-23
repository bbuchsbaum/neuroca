
#' procrustes rotation
proc_rot <- function(X, Y) {
  XY <- crossprod(X, Y)
  sol <- svd(XY)
  A <- sol$v %*% t(sol$u)
  A
}

#' compute bootstrap ratios for a list of matrices
boot_ratio <- function(bootlist) {
  boot.mean <- Reduce("+", bootlist)/length(bootlist)
  boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
  boot.mean/boot.sd	
}


resampleXY <- function(X, Y) {
  ysplit <- split(1:length(Y), Y)
  yidx <- sort(unlist(lapply(ysplit, function(ids) {
      sample(ids, replace=TRUE)
  })))
  
  list(X=X[yidx,], Y=Y[yidx])
  
}

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

resample.bada <- function(x) {
  ## split all indices by category
  Y_by_idx <- split(1:length(x$Y), x$Y)
  resam <- unlist(lapply(Y_by_idx, function(i) sample(i, replace=TRUE)))
  list(Y=x$Y[resam], X=x$X[resam,])
}


#' @export
bootstrap.bada <- function(x, nboot=1000, ncomp=x$ncomp, type=c("projection", "rotated", "unrotated", "split_half")) {
  
  type <- match.arg(type)
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  project.rows <- function(xr) {
    xr %*% x$fit$v[,1:ncomp,drop=FALSE] 
  }
  
  project.cols <- function(xc) {
    t(xc) %*% (x$fit$u %*% diag(1/x$fit$d, nrow=ncomp, ncol=ncomp))
  }
  
  do_split <- function() {
      s1 <- sample(levels(x$S), length(levels(S))/2)
      s2 <- levels(x$S)[! (levels(x$S) %in% s1)]
      s1_idx <- which(x$S %in% s1)
      s2_idx <- which(x$S %in% s2)
      
      X1 <- x$X[s1_idx,,drop=FALSE]
      X2 <- x$X[s2_idx,,drop=FALSE]
      Y1 <- x$Y[s1_idx]
      Y2 <- x$Y[s2_idx]
      
      xboot <- bada(Y1, X1, S=x$S[s1_idx], center=x$center, scale=x$scale)
      A <- proc_rot(scores(x), scores(xboot))
      xboot <- rotate(xboot, A)
      X2b <- group_means(Y2, X2)
      boot_scores <- project(xboot, X2b)
      boot_loadings <- project_cols(xboot, X2b)
      
      list(boot4R=boot_scores,
           boot4C=boot_loadings)
  }
  

  do_jack <- function() {
    res <- lapply(levels(x$S), function(lev) {
      print(lev)
      jack_idx <- which(x$S == lev)
      Xrest <- x$X[-jack_idx,]
      Yrest <- x$Y[-jack_idx]
      
      Xjack <- x$X[jack_idx,]
      Yjack <- x$Y[jack_idx]
      
      Xjack_bc <- group_means(Yjack, Xjack)
      
      xboot <- bada(Yrest, Xrest, S=x$S[-jack_idx], center=x$center, scale=x$scale)
      sign_flip <- sign(diag(t(scores(x)) %*% scores(xboot)))
      
      if (any(sign_flip != 1)) {
        print(sign_flip)
        xboot <- rotate(xboot, diag(sign_flip))
      }
      
      boot_scores <- project(xboot, Xjack_bc)
      boot_loadings <- t(reprocess(xboot, Xjack_bc)) %*% (xboot$fit$u %*% diag(1/xboot$fit$d, nrow=ncomp, ncol=ncomp))
      
      list(boot_scores=boot_scores, boot_loadings=boot_loadings)
    })
    
  }
  
 
  do_boot_projection <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample.bada(x)
    
    ## calculate barycenters of pre-processed data
    Xb <- group_means(resam$Y, reprocess(x, resam$X))
    
    ## project rows and columns of barycenters
    bc <- project_cols(x, Xb)
    br <- project(x, Xb)
    
    list(boot4R=br,
         boot4C=bc)
  }
  
  do_boot_svd <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample.bada(x)
    xboot <- bada(resam$Y, resam$X, x$S, center=x$center, scale=x$scale)
    
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
  
  boot.res <- 
    lapply(1:nboot, function(i) {
      if (type == "projection") {
        do_boot_projection()
      } else if (type == "split_half") {
        do_split()
      } else {
        do_boot_svd()
      }
    })
  
  ret <- list(zboot_scores=boot_ratio(lapply(boot.res, "[[", 1)),
              zboot_loadings=boot_ratio(lapply(boot.res, "[[", 2)))
  
  class(ret) <- c("list", "bootstrap_result")
  ret
  
}




#' @export
bootstrap.mubada <- function(x, niter, ncomp=x$ncomp) {
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  project.rows <- function(xr) {
     xr %*% x$fit$v[,1:ncomp,drop=FALSE] 
   }
   
  project.cols <- function(xc) {
    t(xc) %*% (x$fit$u[,1:ncomp,drop=FALSE]) %*% diag(1/x$fit$d[1:ncomp], nrow=ncomp, ncol=ncomp)
  }
  
  do_boot <- function() {
    ## resample original data with replacement, generating a new dataset
    resam <- resample.mubada(x)
    
    ## calculate barycenters of pre-processed data
    XBary <- do.call(cbind, lapply(1:length(resam$X), function(i) {
      group_means(resam$Y[[i]], reprocess(x, resam$X[[i]], i))
    }))
    
    ## project rows and columns of barycenters
    br <- project.rows(XBary)
    bc <- project.cols(XBary)
    
    list(boot4R=br,
         boot4C=bc)
  }

  
  boot_ratio <- function(bootlist) {
    boot.mean <- Reduce("+", bootlist)/length(bootlist)
    boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
    boot.mean/boot.sd	
  }
  
  boot.res <- 
    lapply(1:niter, function(i) {
      message("bootstrap iteration: ", i)
        do_boot()
    })
  
  
  ret <- list(boot_ratio_rows=boot_ratio(lapply(boot.res, "[[", 1)),
              boot_ratio_cols=boot_ratio(lapply(boot.res, "[[", 2)),
              boot_raw_rows=lapply(boot.res, "[[", 1),
              boot_raw_cols=lapply(boot.res, "[[", 2))
  
  class(ret) <- c("list", "bootstrap_result")
  ret
  
}





#' @export
bootstrap.pca <- function(x, nboot=100, k=x$ncomp) {
  DUt <- t(scores(x))
  n <- dim(DUt)[2]
  
  gen <- function() {
    sidx <- sample(1:ncol(DUt), replace=TRUE)
    DUt[,sidx,drop=FALSE]
  }
  
  res <- boot_svd(nboot=nboot, k=k, loadings(x), gen)
  
  zboot <- do.call(cbind, lapply(1:k, function(ki) {
    res$EVs[[ki]]/res$sdVs[[ki]]
  }))
  
  ret <- list(boot_ratio_vars=zboot, nboot=nboot, k=k)
  class(ret) <- "bootstrap_pca"
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
  
  AsByK <- lapply(1:k, function(ki) {
    do.call(rbind, lapply(res, function(a) {
      a$Ab[,ki]
    }))
  })
  
  
  EAs <- lapply(AsByK, colMeans) #EAs is indexed by k
  EVs <- lapply(EAs, function(EA) v %*% matrix(EA,ncol=1))
  
  varAs <- lapply(AsByK,var) #indexed by k
  
  varVs <- lapply(1:length(AsByK), function(ki) {
    rowSums((v %*% varAs[[ki]]) * v)
  })
  
  sdVs <- lapply(varVs,sqrt)
  list(res=res, EAs=EAs, EVs=EVs, varAs=varAs, sdVs=sdVs)
  
}

#' @keywords internal
boot_svd <- function(nboot, k, v, gen_DUtP) {
  res <- lapply(1:nboot, function(i) {
    DUtP <- gen_DUtP()
    #DUtP <- if(x$center) t(scale(t(DUt[,sidx]),center=TRUE,scale=FALSE)) else DUt[,sidx]
    svd_dutp(DUtP,k)
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

