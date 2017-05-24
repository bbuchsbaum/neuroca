#' @export
supplementary_loadings <- function(x,...) UseMethod("supplementary_loadings")


colCors = function(x, y) { 
  sqr = function(x) x*x
  if (!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y))) {
    stop("Please supply two matrices of equal size.")
  }
  
  x   <- sweep(x, 2, colMeans(x))
  y   <-  sweep(y, 2, colMeans(y))
  cor <- colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
  return(cor)
}

corMat <- function(Xl) {
  pairs <- combn(length(Xl),2)
  M <- matrix(0, length(Xl), length(Xl))
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    rmean <- mean(colCors(t(Xl[[p1]]), t(Xl[[p2]])))
    M[p1,p2] <- rmean
    M[p2,p1] <- rmean
  }
  M
  
}

compute_sim_mat <- function(block_mat, FUN, ...) {
  pairs <- combn(nblocks(block_mat),2)
  M <- matrix(0, nblocks(block_mat), nblocks(block_mat))
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    sim <- FUN(get_block(block_mat, p1), get_block(block_mat, p2), ...)
    M[p1,p2] <- sim
    M[p2,p1] <- sim
  }
  
  M
}


coefficientRV <- function(X, Y, center=TRUE, scale=FALSE) {
  if (dim(X)[[1]] != dim(Y)[[1]]) 
    stop("no the same dimension for X and Y")
  if (dim(X)[[1]] == 1) {
    print("1 configuration RV is  NA")
    rv = NA
  }
  else {
    if (center || scale) {
      Y <- scale(Y, center=center, scale=scale)
      X <- scale(X, center=center, scale=scale)
    }
    W1 <- X %*% t(X)
    W2 <- Y %*% t(Y)
    rv <- sum(diag(W1 %*% W2))/(sum(diag(W1 %*% W1)) * 
                                  sum(diag(W2 %*% W2)))^0.5
  }
  return(rv)
}

normalization_factors <- function(block_mat, type=c("MFA", "RV", "DCor", "None", "RV_MFA")) {
  type <- match.arg(type)
  
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(as.list(block_mat), function(X) 1/(svd_wrapper(X, ncomp=1, method="svds")$d[1]^2)))
  } else if (type == "RV" && nblocks(block_mat) > 2) {
    message("rv normalization")
    smat <- compute_sim_mat(block_mat, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    wts <- abs(svd_wrapper(smat, ncomp=1, method="propack")$u[,1])
  } else if (type == "DCor" && nblocks(block_mat) > 2) {
    message("dcor normalization")
    smat <- compute_sim_mat(block_mat, function(x1,x2) energy::dcor.ttest(x1,x2)$estimate)
    diag(smat) <- 1
    wts <- abs(svd_wrapper(smat, ncomp=1, method="propack")$u[,1])
  } else {
    rep(1, nblocks(block_mat))
  }
}


## take a new design and reduce original data
#' @importFrom assertthat assert_that
reduce_rows <- function(Xlist, Ylist, center=TRUE, scale=FALSE) {
  
  assert_that(length(Xlist) == length(Ylist))
  #normalization <- match.arg(normalization)
  
  ## compute barycenters for each table
  XB <- lapply(seq_along(Xlist), function(i) group_means(Ylist[[i]], Xlist[[i]]))
  
  # center/scale barcycenters of each table
  XBc <- lapply(XB, pre_process, center, scale)
  
  centroids <- lapply(XBc, attr, "scaled:center")
  
  scales <- lapply(XBc, attr, "scaled:scale")
  
  #message("normalization is: ", normalization)
  
  #if (normalization != "Pre") {
    #message("computing normalization factors")
    # compute normalization factor
    
    #alpha <- sqrt(normalization_factors(XBc, normalization))
    #print(alpha)
  #}
  
  # normalize
  #normXBc <- lapply(1:length(alpha), function(i) XBc[[i]] * alpha[i])
  

  # construct matrix of barycenters
  #Xr <- do.call(cbind, normXBc)
  
  Xr <- block_matrix(XBc)
  
  list(Xr=Xr, centroids=centroids, scales=scales)
}



### implement soft-thresholding that spans datasets...? similar to spls?
### musubada is really a bada with block structure -- bada can be engine

#' musubada
#' 
#' @importFrom assertthat assert_that 
#' @importFrom energy dcor.ttest
#' @param Y dependent \code{factor} variable. If All X matrices have same number of rows, Y can be a single factor.
#'        If there are a different nume rof rows (e..g different numbers of replications per subject), Y can be a list of factors.
#' @param Xlist a list of X matrices, one per subject. 
#' @param ncomp
#' @param center
#' @param scale
#' @param svd.method
#' @param normalization
#' @param rank_k use reduce data to k components per block
#' @importFrom assertthat assert_that 
#' @export
musubada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE, svd.method="fast", 
                     normalization=c("MFA", "RV", "None","DCor"), rank_k=NULL) {
  normalization <- normalization[1]

  assert_that(all(sapply(Xlist, is.matrix)))
  
  if (is.null(names(Xlist))) {
    names(Xlist) <- paste0(1:length(Xlist))
  }
   
  table_names <- names(Xlist)
 
  if (is.factor(Y)) {
    ## Y is a single factor, therefore all matrices must have same number of rows.
    assert_that(all(sapply(Xlist, nrow) == nrow(Xlist[[1]])))
    assert_that(length(Y) == nrow(Xlist[[1]]))
    Yl <- replicate(length(Xlist), Y, simplify=FALSE)
  } else if (is.list(Y)) {
    ## Y is a list of factors
    assert_that(all(sapply(Y, is.factor)))
    Yl <- Y
    ## check that all Ys have same levels
    assert_that(length(unlist(Y)) == sum(sapply(Xlist, nrow)))
  }
  
  Y_reps <- 
    lapply(Yl, function(f) {
      m <- outer(f, unique(f), "==")
      apply(m* apply(m,2,cumsum), 1, sum)
    })
  
  YIndices <- rep(1:length(Xlist), sapply(Xlist, nrow))
  
  ## create block variable
  blockInd <- blockIndices(Xlist)
  
  refit <- function(.Y, .Xlist, .ncomp) { 
    musubada(.Y, .Xlist, .ncomp, center, scale, svd.method, normalization) 
  }
  
  permute_refit <- function(.ncomp=ncomp) {
    nreps <- sapply(lapply(Yl, table), min)
    if (any(nreps == 1)) {
      ## permute data
      .Xlist <- lapply(Xlist, function(x) x[sample(1:nrow(x)),])
      refit(Y, .Xlist, .ncomp)
    } else {
      YPerm <- lapply(Yl, function(y) y[sample(1:length(y))])
      refit(YPerm, Xlist, .ncomp)
    }
    
  }
    
  ## average categories within block
  Xreduced <- reduce_rows(Xlist, Yl, center, scale)
  
  Xr <- if (!is.null(rank_k)) {
    is_reduced <- TRUE
    reducer <- reduce_rank(Xreduced$Xr, rank_k)
    reducer$x
  } else {
    is_reduced=FALSE
    reducer <- NULL
    Xreduced$Xr
  }
  
  alpha <- normalization_factors(Xr, type=normalization)
  
  Xr <- block_apply(Xr, function(x,i) x * alpha[i])
  bind <- attr(Xr, "block_ind")
  
  reprocess <- function(newdat, table_index) {
    ## given a new observation(s), pre-process it in the same way the original observations were processed
    if (!is.null(Xreduced$centroids)) {
      ## recenter
      newdat <- sweep(newdat, 2, Xreduced$centroids[[table_index]], "-")
    }
    
    if (!is.null(Xreduced$scales[[table_index]])) {
      ## rescale
      newdat <- sweep(newdat, 2, Xreduced$scales[[table_index]], "/")
    }
    
    if (is_reduced) {
      newdat <- project(reducer, newdat, table_index)
    }
    
    newdat * alpha[table_index]
    
  }
  
  
  YB <- factor(levels(Yl[[1]]), levels=levels(Yl[[1]]))
  
  pca_fit <- pca_core(as.matrix(Xr), ncomp=ncomp, 
                      center=FALSE, 
                      scale=FALSE, 
                      svd.method=svd.method)
  
  ncomp <- length(pca_fit$d)
  
 
  partial_fscores = 
    lapply(1:length(Xlist), function(i) {
      ind <- attr(Xr, "block_indices")[i,]
      length(Xlist) * get_block(Xr, i) %*% pca_fit$v[ind[1]:ind[2],]
    })
  
  sc <- pca_fit$scores
  row.names(sc) <- levels(Yl[[1]])
  
  
  
  result <- list(
    Xlist=Xlist,
    Y=Yl,
    Y_reps=Y_reps,
    conditions=YB,
    XB=Xr,
    scores=sc,
    partial_scores=partial_fscores,
    
    table_contr = do.call(cbind, lapply(1:ncomp, function(i) {
      sapply(1:length(Xlist), function(j) {
        ind <- bind[j,1]:bind[j,2]
        sum(pca_fit$v[ind,i]^2 * alpha[j])
      })
    })),
    
    ntables=length(Xlist),
    ncond=nrow(Xr),
    pca_fit=pca_fit,
    center=center,
    scale=scale,
    ncomp=ncomp,
    blockIndices=bind,
    alpha=alpha,
    normalization=normalization,
    refit=refit,
    table_names=table_names,
    reprocess=reprocess,
    rank_k=rank_k,
    permute_refit=permute_refit
  )
  
  class(result) <- c("musubada")
  result
}

contributions.musubada <- function(x, component=1) {
  out <- lapply(1:x$ntables, function(j) {
    ind <- x$blockInd[j,1]:x$blockInd[j,2]
    contr <- pca_fit$v[ind,component]^2 * x$alpha[j]
    sum(contr)
  })
  
  do.call(cbind, out)
  
}

permute_refit.musubada <- function(x) {
  x$permute_refit()
}

singular_values.musubada <- function(x) {
  x$pca_core$d
}


#' @export
project_table <- function(x, supY, supX, ncomp, ...) UseMethod("project_table")

project_copy <- function(x, ...) UseMethod("project_copy")

project_copy.musubada <- function(x, Ylist, Xlist) {
  ## reduce but no centering or scaling
  XB <- lapply(1:length(Xlist), function(i) group_means(Ylist[[i]], Xlist[[i]]))
  
  XBlocks <- lapply(1:length(XB), function(i) x$reprocess(XB[[i]], i))
  browser()
  partial_scores = 
    lapply(1:length(XBlocks), function(i) {
      ind <- x$blockInd[i,1]:x$blockInd[i,2]
      length(XBlocks) * XBlocks[[i]] %*% x$pca_fit$v[ind,]
    })
  
  lds <- do.call(cbind, XBlocks %*% x$pca_fit$u)
  
  scores <- Reduce("+", partial_fscores)/length(partial_fscores)
  list(scores=scores, partial_scores=partial_scores, loadings=lds)
  
}

#' @importFrom assertthat assert_that 
#' @export
project_table.musubada <- function(x, supY, supX, ncomp=x$ncomp) { #table_index=NULL, table_name = NULL) {
  assert_that(ncomp >= 1)
  assert_that(is.factor(supY))
  assert_that(length(levels(supY)) == x$ncond)
  
  ## pre-process new table
  suptab <- block_reduce(list(supX), list(supY), normalization=x$normalization, scale=x$scale, center=x$center)
  
  ## compute supplementary Q
  Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] %*% diag(1/x$pca_fit$d[1:ncomp]))
  
  ##Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] )
  
  ## compute supplementary factor scores
  Fsup <- x$ntables * (suptab$Xr %*% Qsup) 
  list(scores=Fsup, loadings=Qsup)
}

#' @export
supplementary_loadings.musubada <- function(x, suptab, ncomp=x$ncomp) {
  suptab <- x$pre_process(suptab)
  Qsup <- t(suptab) %*% (x$pca_fit$u[,1:ncomp,drop=FALSE]) %*% diag(1/x$pca_fit$d[1:ncomp])
}


#' @export
loadings.musubada <- function(x, table_index=1:nrow(x$blockIndices), comp=1:ncol(x$pca_fit$v)) {
  do.call(cbind, lapply(table_index, function(i) {
    ind <- x$blockIndices[i,1]:x$blockIndices[i,2]
    1/x$alpha[i] * x$pca_fit$v[ind, comp,drop=FALSE]
  }))
}

#' @export
correlations.musubada <- function(x, table_index=1:nrow(x$blockIndices), comp=1:ncol(x$pca_fit$u)) {
  do.call(cbind, lapply(table_index, function(i) {
    ind <- x$blockIndices[i,1]:x$blockIndices[i,2]
    cor(x$XB[,ind,drop=FALSE], x$pca_fit$u[, comp])
  }))
}


# rename "embed_table"?

#' @export
#' @importFrom assertthat assert_that
supplementary_predictor.musubada <- function(x, supX, supY, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
                                             ncomp=x$ncomp) {
  assert_that(length(supY) == nrow(supX))
  type <- match.arg(type)
  ## pre-process new table
  suptab <- block_reduce(list(supX), list(supY), normalization=x$normalization, scale=x$scale, center=x$center)
  
  ## compute supplementary Q
  Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE] %*% diag(1/x$pca_fit$d[1:ncomp]))
  
  ##Qsup <- t(suptab$Xr) %*% (x$pca_fit$u[, 1:ncomp, drop=FALSE])
  
  predfun <- function(newdat) {
    if (!is.null(suptab$centroids)) {
      newdat <- sweep(newdat, 2, suptab$centroids[[1]], "-")
    }
    if (!is.null(suptab$scales[[1]])) {
      sweep(newdat, 2, suptab$scales[[1]], "/")
    }
    
    newdat <- newdat * suptab$alpha
    fscores <- x$ntables * newdat %*% Qsup
    scorepred(fscores, x$scores, type, ncomp)
  }
}

 
#' @importFrom abind abind
project.musubada <- function(x, newX=NULL, ncomp=x$ncomp, table_index=1:x$ntables) {
  if (length(table_index) > 1) {
    assert_that(is.list(newX) && length(newX) == length(table_index))
  }
  
 
  # if (is.null(newX)) {
  #   ## project each replication onto the subject's components
  #   lapply(table_index, function(i) {
  #     yr <- x$Y_reps[[i]]
  #     ysplit <- split(1:length(yr), yr)
  #     newdata <- x$reprocess(x$Xlist[[i]],i)
  #     #browser()
  #     lv <- lapply(ysplit, function(idx) {
  #       xnewdat <- newdata[idx,]
  #       #t(t(xnewdat) %*% x$pca_fit$u[, 1:ncomp])
  #       t(t(xnewdat) %*% (x$pca_fit$u[, 1:ncomp] %*% diag(1/x$pca_fit$d[1:ncomp])))
  #     })
  #     
  #     abind::abind(lv, along=3)
  #     
  #   })
  #   
  # } else {
  
  
  
  if (is.vector(newX)) {
    assert_that(length(table_index) == 1)
    newX <- list(matrix(newX, 1, length(newX)))
  }
   
  if (is.matrix(newX)) {
    assert_that(length(table_index) == 1)
    newX <- list(newX)
  }
    
  ## project new data-point 
  res <- lapply(1:length(table_index), function(i) {
    
    tbind <- table_index[i]
    print(tbind)
    xnewdat <- x$reprocess(newX[[i]], tbind)
    ind <- x$blockIndices[tbind,]
    xnewdat %*% x$pca_fit$v[ind[1]:ind[2], 1:ncomp] * x$ntables
  })
  
  names(res) <- paste0("table_", table_index)
  res
}
  


## project from existing table
#' @export
predict.musubada <- function(x, newdata, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
                                    ncomp=x$ncomp, table_index=1:x$ntables) {
  type <- match.arg(type)
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  fscores <- if (length(table_index) == x$ntables) {
    assert_that(ncol(newdata) == sum(sapply(x$Xlist, ncol)))
    Reduce("+", lapply(table_index, function(i) {
      ind <- x$blockIndices[i,]
      Xp <- x$reprocess(newdata[, ind[1]:ind[2]], i)
      fscores <- Xp %*% x$pca_fit$v[ind[1]:ind[2],,drop=FALSE]
    }))
  } else if (length(table_index) == 1) {
    ind <- x$blockIndices[table_index,]
    Xp <- x$reprocess(newdata, table_index)
    Xp %*% x$pca_fit$v[ind[1]:ind[2], 1:ncomp, drop=FALSE] * x$ntables
  }
  
  
  scorepred(fscores, x$scores, type=type, ncomp=ncomp)
  
}
