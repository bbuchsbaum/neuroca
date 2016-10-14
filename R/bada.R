#' @export
bada <- function(Y, X, ncomp=length(levels(Y)), center=TRUE, scale=FALSE, svd.method="base", strata=NULL) {
  assert_that(is.factor(Y))
  assert_that(length(Y) == nrow(X)) 
  
  if (!is.null(strata)) {
    assert_that(length(strata) == length(Y))
    strata <- as.factor(strata)
  }
  
  
  XBc <- if (any(table(Y) > 1) && is.null(strata)) {
    ## compute barycenters, no strata
    XB <- group_means(Y, X)
    scale(XB, center=center, scale=scale)
  } else if (!is.null(strata)) {
    ## compute barycenters by strata
    XBlock <- lapply(levels(strata), function(i) {
      idx <- which(strata==i)
      scale(group_means(Y[idx], X[idx,]), center=center, scale=scale)
    })
   
    Reduce("+", XBlock)/length(XBlock)
  } else {
    scale(X, center=center, scale=scale)
  }
  
  pcres <- pca_core(t(XBc), ncomp, center=FALSE, scale=FALSE, svd.method)
  scores <- pcres$scores
  row.names(scores) <- levels(Y)
  
  refit <- function(.Y, .X, .ncomp) { bada(.Y, .X, .ncomp, center, scale, svd.method, strata) }
  
  permute_refit <- function(.ncomp=ncomp) {
    if (is.null(strata)) {
      idx <- sample(seq_along(Y))
      refit(Y[idx], X, .ncomp)
    } else {
      Xperm <- do.call(rbind, lapply(levels(strata), function(lev) {
        idx <- which(lev == strata)
        X[sample(idx),]
      }))
      refit(Y, Xperm, .ncomp)
    }
  }
  
  ret <- list(Y=Y,X=X,ncomp=ncomp, condMeans=XBc, center=center, scale=scale, pre_process=apply_scaling(XBc), 
              svd.method=svd.method, scores=scores, v=pcres$v, u=pcres$u, d=pcres$d, refit=refit, permute_refit=permute_refit, 
              strata=strata, XBlock=XBlock)
 
  class(ret) <- c("pca_result", "bada_result")
  ret
  
}


singular_values.bada_result <- function(x) {
  x$d
}

partial_scores.bada_result <- function(x) {
  if (is.null(x$strata)) {
    x$scores
  } else {
    levs <- levels(x$strata)
    lapply(1:length(levs), function(i) {
      x$XBlock[[i]] %*% x$u
    })
  }
}

