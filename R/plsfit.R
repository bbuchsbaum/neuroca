
rbindN <- function(vec, n) {
  do.call("rbind", replicate(n, vec, simplify=FALSE)) 
}

mean_center <- function(Y, X) {
  G <- model.matrix(~ Y - 1)
  colnames(G) <- levels(Y)

  #GW <- G /(ncol(G) * rbindN(colSums(G),nrow(G))            
  #weights <- as.matrix(rep(1, ncol(X))*(1/ncol(X)))
  #rowMasses <- rowSums(GW)
  #colMasses <- colSums(GW)
  GW <- G/colSums(G)
  R <- crossprod(X, GW)
  centroid <-  rowMeans(R)
  Rcent <- sweep(R, 1, centroid)
  
  list(Rcent=Rcent, Ymat=G, centroid=centroid)
}

svd.pls <- function(XC, ncomp=min(dim(XC)), method=c("base", "fast", "irbla")) {
  res <- switch(method[1],
                base=svd(XC),
                fast=corpcor:::fast.svd(XC),
                irlba=irlba:::irlba(XC, nu=min(ncomp, min(dim(XC)) -3), nv=min(ncomp, min(dim(XC)) -3)))
  
  k <- which(res$d > .Machine$double.eps)
  k <- min(ncomp, length(k))
  res$d <- res$d[1:k]
  res$u <- res$u[,1:k]
  res$v <- res$v[,1:k]
  res$ncomp <- k
  res
}

project.rows <- function(xr, svd.fit, centroid=NULL) {
  if (!is.null(centroid)) {
    sweep(xr, 2, centroid) %*% svd.fit$u  
  } else {
    xr %*% svd.fit$u
  }
}

project.cols <- function(xc, svd.fit) {
  xc %*% svd.fit$v
}


split.sample <- function(y) {
  idx <- 1:length(y)
  f <- split(idx, y)
  
  res <- lapply(f, function(x) {
    xs <- sample(x)
    M <- as.integer(length(xs)/2)
    list(a=xs[1:M], b=xs[(M+1):length(xs)])
  })
  
  s1 <- unlist(lapply(res, "[[", 1))
  s2 <- unlist(lapply(res, "[[", 2))
  list(s1=s1, s2=s2)
}

mcen.pls.boot <- function(X, Y, svd.fit, ncomp=svd.fit$ncomp, boot.iter=100, strata=NULL) {
  do_boot <- function(indices) {
    XBoot <- X[indices,]
    YBoot <- Y[indices]
    bary <- mean_center(YBoot, XBoot)
    
    ### these are in factor space
    br <- project.rows(t(bary$Rcent), svd.fit)
    bc <- project.cols(bary$Rcent, svd.fit)
    
    list(boot4R=br,
         boot4C=bc)
  }
  
  boot.ratio <- function(bootlist) {
    boot.mean <- Reduce("+", bootlist)/length(bootlist)
    boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
    boot.mean/boot.sd	
  }
  
  .resample <- function(x, ...)  { x[sample.int(length(x), ...)] }
  
  boot.res <- if (is.null(strata)) {
    lapply(1:boot.iter, function(i) {
      indices <- lapply(levels(Y), function(lev) {
        sample(which(model$Y == lev), replace=TRUE)
      })
      
      do_boot(sort(unlist(indices)))
    })
  } else {
    
    lapply(1:boot.iter, function(i) {
      boot.strat <- sample(levels(strata), replace=TRUE)
      indices <- lapply(boot.strat, function(stratum) {
        unlist(lapply(levels(Y), function(lev) {
          .resample(which(Y == lev & strata==stratum), replace=TRUE)					
        }))
      })
      
      do_boot(sort(unlist(indices)))
    })		
  }
  
  ret <- list(boot.ratio.R=boot.ratio(lapply(boot.res, "[[", 1)),
                   boot.ratio.C=boot.ratio(lapply(boot.res, "[[", 2)),
                   boot.raw.R=lapply(boot.res, "[[", 1),
                   boot.raw.C=lapply(boot.res, "[[", 2))
    
  
}



mcen.pls.cv <- function(X, Y, svd.fit, cv.iter=200) {
  cv.res <- do.call(rbind, mclapply(1:cv.iter, function(i) {
    print(i)
    halves <- split.sample(Y)
    fit1 <- svd.pls(mean_center(Y[halves$s1],X[halves$s1,])$Rcent, method="irlba", ncomp=ncomp)
    fit2 <- svd.pls(mean_center(Y[halves$s2],X[halves$s2,])$Rcent, method="irlba", ncomp=ncomp)
    f1 <- fit1$v %*% diag(fit1$d)
    f2 <- fit2$v %*% diag(fit2$d)
    sapply(1:ncomp, function(k) {
      coeffRV(f1[,1:k, drop=FALSE], f2[,1:k,drop=FALSE])$rvstd
    })
  }))
  
  ### is first component significant?
  p1 <- wilcox.test(res[,1])$p.value
  nsig <- if (p1 > .05) {
    ## it's not
    0
  } else {
    ## rank RV coefficients for each iteration
    cv.ranks <- apply(res,1,rank)
    
    ## compute sign rank test for successive compenents
    pvals <- sapply(1:(ncomp-1), function(i) {
      wilcox.test(cv.ranks[i,], cv.ranks[i+1,])$p.value     
    })
    
    if (all(pvals > .05)) {
      1
    } else {
      ## find first sequential comparison that is not significant
      wp <- which(pvals  > .05)
      wp[1] 
    }
  }
  
  list(cv=res, nsig=nsig) 
}
                
mcen.pls <- function(Y, X, ncomp=5, strata, cv=TRUE, cv.iter=200, boot=TRUE, boot.iter=200, svd.method="fast") {
  if (!is.factor(Y)) {
    warning("converting Y to factor")
    Y <- as.factor(Y)
  }
  
  bary <- mean_center(Y,X)
  svd.fit <- svd.pls(bary$Rcent, ncomp, mehod=method)
  
  brainScores <- project.rows(X, bary$centroid, svd.fit)
  designScores <- svd.fit$v
  
  brainSaliences <- svd.fit$u
  designSaliences <- svd.fit$v
  
  if (cv && ncomp > 1) {
    cv.res <- mcen.pls.cv(X, Y, svd.fit, cv.iter=cv.iter)          
  }
  if (boot && ncomp >= 1) {
    boot.res <- mcen.pls.boot(X, Y, svd.fit, ncomp=svd.fit$ncomp, boot.iter=boot.iter, strata=strata)    
  }
 
}

setMethod("brainSaliences", "plsfit", function(x) {
  x@brainSaliences
})

setMethod("designSaliences", "plsfit", function(x) {
  x@designSaliences
})

setMethod("brainScores", "plsfit", function(x) {
  x@brainScores
})

setMethod("designScores", "plsfit", function(x) {
  x@designScores
})







