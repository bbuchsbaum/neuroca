
rbindN <- function(vec, n) {
  do.call("rbind", replicate(n, vec, simplify=FALSE)) 
}

mean_center <- function(Y, X) {
  G <- model.matrix(~ Y - 1)
  colnames(G) <- levels(Y)

  GW <- G/colSums(G)
  R <- crossprod(X, GW)
  centroid <-  rowMeans(R)
  Rcent <- sweep(R, 1, centroid)
  
  list(Rcent=Rcent, Ymat=G, centroid=centroid)
}

svd.pls <- function(XC, ncomp=min(dim(XC)), method=c("base", "fast", "irbla", "SPC", "PMD"), sumabsu=sqrt(nrow(XC)/4), sumabsv=sqrt(ncol(XC)/4)) {
  res <- switch(method[1],
                base=svd(XC),
                fast=corpcor:::fast.svd(XC),
                irlba=irlba:::irlba(XC, nu=min(ncomp, min(dim(XC)) -3), nv=min(ncomp, min(dim(XC)) -3)),
                SPC=as.list(unclass(SPC(XC, sumabsv, K=ncomp))),
                PMD=as.list(unclass(PMD(XC, sumabsv=sumabsv, sumabsu=sumabsu, K=ncomp))))
  
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

#mcen.pls.jack

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
  p1 <- wilcox.test(cv.res[,1])$p.value
  nsig <- if (p1 > .05) {
    ## it's not
    0
  } else {
    
    ## rank RV coefficients for each iteration
    cv.ranks <- apply(cv.res,1,rank)
    
    ## compute sign rank test for successive compenents
    pvals <- sapply(1:(ncomp-1), function(i) {
      D <- cv.ranks[i+1,] - cv.ranks[i,]
      if (sum(D > 0)/length(D) < .5) {
        1
      } else {
        wilcox.test(cv.ranks[i,], cv.ranks[i+1,])$p.value  
      }
    })
    
    if (all(pvals > .05)) {
      1
    } else {
      ## find first sequential comparison that is not significant
      wp <- which(pvals  > .05)
      wp[1] 
    }
  }
  
  list(cv=cv.res, nsig=nsig) 
}
                
pls.meancen <- function(Y, X, ncomp=5, strata=NULL, cv=TRUE, cv.iter=200, boot=TRUE, boot.iter=200, svd.method="fast") {
  if (!is.factor(Y)) {
    warning("converting Y to factor")
    Y <- as.factor(Y)
  }
  
  if (length(Y) != nrow(X)) {
    stop(paste("length of Y: ", length(Y), " must equal number of rows in X: ", nrow(X)))
  }
  
  if (is.null(strata)) {
    strata <- factor(rep(1, length(Y)))
  }
  
  
  
  bary <- mean_center(Y,X)
  svd.fit <- svd.pls(bary$Rcent, ncomp, method=svd.method)
  
  brainScores <- project.rows(X, svd.fit, bary$centroid)
  designScores <- svd.fit$v
  
  brainSaliences <- svd.fit$u
  designSaliences <- svd.fit$v
  
  nsigcomp <- if (cv && ncomp > 1) {
    cv.res <- mcen.pls.cv(X, Y, svd.fit, cv.iter=cv.iter)      
    cv.res$nsig
  } else {
    0
  }
  
  
  ### make crossvalidation and bootstapping separate.
  boot.res <- if (boot && ncomp >= 1) {
    mcen.pls.boot(X, Y, svd.fit, ncomp=svd.fit$ncomp, boot.iter=boot.iter, strata=strata)    
  } else {
    list(boot.ratio.R=matrix(),
         boot.ratio.C=matrix(),
         boot.raw.R=matrix(),
         boot.raw.C=matrix())
  }

  
  #representation(d="numeric",brainSaliences="matrix",designSaliences="matrix",brainScores="matrix", designScores="matrix", nsigcomp="integer",
  #               ZBootDesign="matrix", ZBootBrain="matrix", strata="factor"),
  
  new("plsfit",
      Y=Y,
      X=X,
      d=svd.fit$d,
      brainSaliences=svd.fit$u,
      designSaliences=svd.fit$v,
      brainScores=brainScores,
      designScores=svd.fit$v,
      nsigcomp=nsigcomp,
      Centroids=bary$Rcent,
      GlobalCentroid=bary$centroid,
      svd=svd.fit,
      ZBootDesign=boot.res$boot.ratio.R,
      ZBootBrain=boot.res$boot.ratio.C,
      strata=strata)
}

setMethod("brainSaliences", "plsfit", function(object) {
  object@brainSaliences
})

setMethod("designSaliences", "plsfit", function(object) {
  object@designSaliences
})

setMethod("brainScores", "plsfit", function(object) {
  object@brainScores
})

setMethod("designScores", "plsfit", function(object) {
  object@designScores
})

setMethod("contributions", "plsfit", function(object, by=NULL) {
  if (!is.null(by)) {
    res <- apply(object@brainScores, 2, function(vals) {
      aggregate(vals ~ by, FUN=function(vals) sum(vals^2))$vals
    })
    row.names(res) <- levels(by)
    apply(res^2, 2, function(vals) vals/sum(vals))
  } else {
    apply(object@brainScores^2, 2, function(vals) vals/sum(vals))
  }
})


setMethod("cv", signature(object="plsfit", nfolds="numeric", ncomp="numeric"),
          function(object, nfolds, ncomp, svd.method="fast", featSel=NULL, metric="accuracy") {
            foldList <- createFolds(object@Y, nfolds)
            res <- lapply(foldList, function(idx) {
              Xtrain <- object@X[-idx,]
              Xtest <- object@X[idx,]
              Ytrain <- object@Y[-idx]
              Ytest <- object@Y[idx]
              if (!is.null(featSel)) {
                
                keep.idx <- which(featSel(Ytrain, Xtrain))
                print(length(keep.idx))
              } else {
                keep.idx <- 1:ncol(Xtrain)
              }
              
              pfit <- pls.meancen(Ytrain, Xtrain[,keep.idx], ncomp=ncomp, cv=FALSE, boot=FALSE)            
              ypred <- predict(pfit, newdata=Xtest[,keep.idx], ncomp=ncomp)
              
              #c10 <- ypred[[1]]
              #proj <- c10$projection
              #ylevs <- levels(object@Y)
              #Fscores <- pfit@svd$v %*% diag(pfit@d)
              #SSTot <- apply(proj, 1, function(vals) sum(vals ^2))
              #SSWithin <- sapply(1:nrow(proj), function(i) {
              #  sum((proj[i,] - Fscores[which(Ytest[i] == ylevs),])^2)
              #})
                     
            })
            
            if (metric == "accuracy") {
              R <- do.call(rbind, lapply(res, function(M) {
                do.call(cbind, lapply(M, function(m) m$class))
              }))                     
              R[order(unlist(foldList)),]
              
            } else if (metric == "distanceRank") {
              ylevs <- levels(object@Y)
             
              rankD <- do.call(rbind, lapply(1:length(foldList), function(i) {
                idx <- foldList[[i]]
                yidx <- sapply(object@Y[idx], function(y) which(y == ylevs))
                M <- res[[i]]
                D1 <- do.call(cbind, lapply(M, function(m) {
                  m$distanceRank[cbind(yidx,1:ncol(m$distanceRank))]
                }))
              }))                   
              rankD <- rankD[order(unlist(foldList)),]         
            } 
            
          })

setMethod("predict", signature(object="plsfit", newdata="matrix", ncomp="numeric"),
          function(object, newdata, ncomp) {
            if (ncol(newdata) != nrow(object@Centroids)) {
              stop(paste("newdata must have ", nrow(object@Centroids), "columns"))
            }
            
            ## remove global mean
            R <- sweep(newdata, 2, object@GlobalCentroid)
            
            Fscores <- project.rows(R, object@svd)
            Forig <- object@svd$v %*% diag(object@d)
            
            res <- lapply(ncomp, function(nc) {
              
              D <- rdist(Fscores[,1:nc], Forig[,1:nc])
              D2 <- D^2
              min.d <- apply(D2, 1, which.min)
              rank.d <- apply(D2, 1, rank)
              classPred <- levels(object@Y)[min.d]
              
              list(class=classPred, dsquared=D2, projection=Fscores, distanceRank=rank.d)
            })
                        
            names(res) <- paste0("Components", ncomp)
            res           
          })






