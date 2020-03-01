
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
bootstrap.mubada <- function(x, niter, nboot=100, ncomp=x$ncomp,type=c("projection", "rotated", "unrotated")) {
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
bootstrap.bada <- function(x, nboot=1000, ncomp=x$ncomp, type=c("projection", "rotated", "unrotated")) {
  
  type <- match.arg(type)
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
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
      boot_loadings <- t(reprocess.bada(xboot, Xjack_bc)) %*% (xboot$fit$u %*% diag(1/xboot$fit$d, nrow=ncomp, ncol=ncomp))
      
      list(boot_scores=boot_scores, boot_loadings=boot_loadings)
    })
    
  }
  
  do_perm <- function() {
    split_idx <- split(1:length(x$Y), x$S)
    Yperm <- unlist(lapply(split_idx, function(idx) {
      sample(x$Y[idx])
    }))

   
    bperm <- bada(Yperm, x$X, S=x$S, center=x$center, scale=x$scale)
    
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
bootstrap.pca <- function(x, nboot=100, k=x$ncomp) {
  DUt <- t(scores(x))
  n <- dim(DUt)[2]
  
  gen <- function() {
    sidx <- sample(1:ncol(DUt), replace=TRUE)
    list(DUt=DUt[,sidx,drop=FALSE], idx=sidx)
  }
  
  res <- boot_svd(nboot=nboot, k=k, loadings(x), gen)

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
  
  EScores <- lapply(ScoresByK, function(s) {
    apply(s, 2, mean, na.rm=TRUE)
  })
  
  sdScores <- lapply(ScoresByK, function(s) {
    apply(s, 2, sd, na.rm=TRUE)
  })
  
  
  
  EAs <- lapply(AsByK, colMeans) #EAs is indexed by k
  EVs <- lapply(EAs, function(EA) v %*% matrix(EA,ncol=1))
  
  varAs <- lapply(AsByK,var) #indexed by k
  
  varVs <- lapply(1:length(AsByK), function(ki) {
    rowSums((v %*% varAs[[ki]]) * v)
  })

  sdVs <- lapply(varVs,sqrt)
  list(res=res, EAs=EAs, EVs=EVs, varAs=varAs, sdVs=sdVs, EScores=EScores, sdScores=sdScores)
  
}

#' @keywords internal
boot_svd <- function(nboot, k, v, gen_DUtP) {
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

