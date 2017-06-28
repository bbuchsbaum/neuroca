
resampleXY <- function(X, Y) {
  ysplit <- split(1:length(Y), Y)
  yidx <- sort(unlist(lapply(ysplit, function(ids) {
      sample(ids, replace=TRUE)
  })))
  
  list(X=X[yidx,], Y=Y[yidx])
  
}

resample.musu_bada <- function(x) {
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
bootstrap.musu_bada <- function(x, niter, ncomp=x$ncomp) {
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  project.rows <- function(xr) {
     xr %*% x$pca_fit$v 
   }
   
  project.cols <- function(xc) {
    t(xc) %*% (x$pca_fit$u[,1:ncomp,drop=FALSE]) %*% diag(1/x$pca_fit$d[1:ncomp])
  }
  
  do_boot <- function() {
    resam <- resample.musu_bada(x)
    XBary <- do.call(cbind, lapply(1:length(resam$X), function(i) {
      group_means(resam$Y[[i]], x$reprocess(resam$X[[i]], i))
    }))
    
    
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


# 
# bootstrap.musubada_result <- function(x, niter) {
#   
# }
# 
# 
# boot_ratio <- function(bootlist) {
#   boot.mean <- Reduce("+", bootlist)/length(bootlist)
#   boot.sd <- sqrt(Reduce("+", lapply(bootlist, function(mat) (mat - boot.mean)^2))/length(bootlist))
#   boot.mean/boot.sd	
# }
# 
# .resample <- function(x, ...)  { x[sample.int(length(x), ...)] }
# 
# sample_indices <- function(Y) {
#   sort(unlist(lapply(levels(x$Y), function(lev) {
#     sample(which(x$Y == lev), replace=TRUE)
#   })))
# }
# 
# boot.res <- 
#   lapply(1:niter, function(i) {
#     message("bootstrap iteration: ", i)
#     if (is.null(x$strata)) {
#       row_indices <- sample_indices(x$Y)
#       do_boot(sort(unlist(row_indices)))
#     } else {
#       levs <- levels(x$strata)
#       slevs <- sample(levs, length(levs), replace=TRUE)
#       ind <- sort(sapply(slevs, function(lev) which(lev == levs)))
#       ret <- do_stratified_boot(ind)
#     }
#   })
# 
# 
# do_boot <- function(indices) {
#   XBoot <- x$X[indices,]
#   YBoot <- x$Y[indices]
#   
#   #Xc <- x$pre_process(XBoot)
#   
#   ## barycenters of bootstrap samples
#   #XB <- x$reduce(XBoot, YBoot)
#   XB <- group_means(YBoot, XBoot)
#   
#   ## apply pre-processing from full model
#   XBc <- x$pre_process(XB) 
#   
#   br <- project.rows(XBc, x$u)
#   bc <- project.cols(t(XBc), x$v)
#   
#   list(boot4R=br,
#        boot4C=bc)
# }