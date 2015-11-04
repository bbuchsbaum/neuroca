#' @importFrom assertthat assert_that 

#' @export
musubada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast.svd", normalization="MFA") {
  assert_that(all(sapply(Xlist, is.matrix)))
  assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
  
  assert_that(is.factor(Y))
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  ## compute barycenters for each table
  XB <- lapply(1:length(Xlist), function(i) group_means(Y[YIndices==i], Xlist[[i]]))
  # center/scale barcycenters of each table
  XBc <- lapply(XB, function(X) scale(X, scale=scale, center=center))
  
  ## create block variable
  blockVar <- factor(unlist(lapply(1:length(Xlist), function(i) rep(i, ncol(Xlist[[i]])))))
  
  ## get column indices associated with each block
  blockIndices <- lapply(1:length(Xlist), function(i) which(blockVar == i))
  
  ## compute scaling factor based on reciprocal of first squared singular value
  alpha <- lapply(XBc, function(X) 1/svd.wrapper(X, ncomp=ncomp, method="irlba")$d[1]^2)
  
  ## scale each table by 'alpha'
  normXBc <- lapply(1:length(XBc), function(i) XBc[[i]] * alpha[[i]])
  
  YB <- factor(row.names(normXBc[[1]]))
  pls_res <- plscorr_da(YB, do.call(cbind, normXBc), ncomp=ncomp, center=FALSE, scale=FALSE)
  
  partial_fscores <- lapply(1:length(Xlist), function(i) {
    ind <- blockIndices[[i]]
    normXBc[[i]] %*% pls_res$u[ind,]
  })
  
 table_contr <- do.call(cbind, lapply(1:ncomp, function(i) {
   sapply(1:length(blockIndices), function(j) {
     sum(pls_res$u[blockIndices[[j]],i]^2)
   })
 }))
     
  
  
  
}