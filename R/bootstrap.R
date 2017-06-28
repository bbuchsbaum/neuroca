
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
    ## resample original data with replacement, generating a new dataset
    resam <- resample.musu_bada(x)
    
    ## calculate barycenters of pre-processed data
    XBary <- do.call(cbind, lapply(1:length(resam$X), function(i) {
      group_means(resam$Y[[i]], x$reprocess(resam$X[[i]], i))
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

