
resampleXY <- function(X, Y) {
  ysplit <- split(1:length(Y), Y)
  yidx <- sort(unlist(lapply(ysplit, function(ids) {
      sample(ids, replace=TRUE)
  })))
  
  list(X=X[yidx,], Y=Y[yidx])
  
}

resample.musubada <- function(x, resample_blocks=TRUE, resample_observations=TRUE) {
  ## check for multiple repetitions per Y levels per block

  if (resample_blocks) {
    table_idx <- sort(sample(1:x$ntables, replace=TRUE))
    Xresam <- x$Xlist[table_idx]
    Yresam <- x$Y[table_idx]
  } else {
    Xresam <- x$Xlist
    Yresam <- x$Y
  }
  
  ret <- lapply(1:length(Yresam), function(i) {
    resampleXY(Yresam[[i]], Xresam[[i]])
  })
  
  Xret <- lapply(ret, "[[", "X")
  Yret <- lapply(ret, "[[", "Y")
  list(X=Xret, Y=Yret)
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