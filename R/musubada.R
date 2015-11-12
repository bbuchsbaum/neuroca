#' @importFrom assertthat assert_that 



### implement soft-thresholding that spans datasets...? similar to spls?
#' @importFrom assertthat assert_that 
#' @export
musubada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast.svd", normalization="MFA") {
  assert_that(all(sapply(Xlist, is.matrix)))
  #assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
  
  assert_that(is.factor(Y))
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  ## create block variable
  blockVar <- factor(unlist(lapply(1:length(Xlist), function(i) rep(i, ncol(Xlist[[i]])))))
  
  ## get column indices associated with each block
  blockIndices <- lapply(1:length(Xlist), function(i) which(blockVar == i))
  
  pre_process <- function(X) {
    Xc <- scale(X, center=center, scale=scale)
    
    if (normalization == "MFA") {
      alpha <- 1/svd.wrapper(X, ncomp=1, method="irlba")$d[1]
      Xc <- Xc * alpha
    } else {
      alpha <- 1
    }
    
    attr(Xc, "alpha") <- alpha
    Xc
  }
  
  ## compute barycenters for each table
  XB <- lapply(1:length(Xlist), function(i) group_means(Y[YIndices==i], Xlist[[i]]))
  
  # center/scale barcycenters of each table
  normXBc <- lapply(XB, pre_process)
  
  
  
  ## compute scaling factor based on reciprocal of first squared singular value
  #alpha <- lapply(XBc, function(X) 1/svd.wrapper(X, ncomp=1, method="irlba")$d[1])
  
  ## scale each table by 'alpha'
  #normXBc <- lapply(1:length(XBc), function(i) XBc[[i]] * alpha[[i]])
  
  YB <- factor(row.names(normXBc[[1]]))
  pls_res <- plscorr_da(YB, do.call(cbind, normXBc), ncomp=ncomp, center=FALSE, scale=FALSE)
  
  
  result <- list(
    scores=pls_res$scores,
    partial_fscores = 
      lapply(1:length(Xlist), function(i) {
        ind <- blockIndices[[i]]
        length(Xlist) * normXBc[[i]] %*% pls_res$u[ind,]
      }),
    
    table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
      sapply(1:length(blockIndices), function(j) {
        sum(pls_res$u[blockIndices[[j]],i]^2)
      })
    })),
    ntables=length(Xlist),
    ncond=nrow(XB[[1]]),
    plsfit=pls_res,
    center=center,
    scale=scale,
    ncomp=ncomp,
    alpha=sapply(normXBc, function(x) attr(x, "alpha")),
    pre_process=pre_process
  )
  
  class(result) <- "musubada_result"
  result
}

#' @export
project_table <- function(x, newdata, ncomp, ...) UseMethod("project_table")


#' @importFrom assertthat assert_that 
#' @export
project_table.musubada_result <- function(x, suptab, ncomp=x$ncomp) {
  assert_that(is.matrix(suptab))
  assert_that(nrow(suptab) == x$ncond )
  assert_that(ncomp >= 1)
  
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$plsfit$v[,1:ncomp,drop=FALSE] %*% diag(1/x$plsfit$d[1:ncomp]))
  Fsup <- x$ntables * suptab %*% Qsup 
}

#' @export
predict.musubada_result <- function(x, newdata, type=c("class", "prob", "scores", "crossprod"), ncomp=x$ncomp) {
  #if (is.vector(newdata) && (length(newdata) == ncol(x$condMeans))) {
  #  x <- matrix(newdata, nrow=1, ncol=ncol(x$condMeans))
  #
  #}
  
  assert_that(is.matrix(newdata))
  assert_that(nrow(newdata) == x$ncond)
  assert_that(type[1] %in% c("class", "prob", "scores", "crossprod"))
  
  fscores <- project_table(x, newdata, ncomp)
  
  if (type == "scores") {
    fscores
  } else if (type == "crossprod") {
    tcrossprod(fscores, x$scores[,1:ncomp,drop=FALSE])
  }else if (type =="class") {
    D <- rdist(fscores, x$scores[,1:ncomp,drop=FALSE])
    D2 <- D^2
    min.d <- apply(D2, 1, which.min)
    levels(x$plsfit$Y)[min.d]
  } else {
    ## type is 'prob'
    probs <- tcrossprod(fscores, x$scores[,1:ncomp,drop=FALSE])
    maxid <- apply(probs,1,which.max)
    maxp <- probs[cbind(1:nrow(probs), maxid)]
    probs <- exp(probs - maxp)
    probs <- zapsmall(probs/rowSums(probs))
  }
  
}
