



#' @export
#' @param Y \code{factor} variable defining the groups
#' @param X \code{matrix} defining the matrix data to be group-wise averaged
group_means <- function(Y, X) {
  assertthat::assert_that(is.character(Y) || is.factor(Y))
  if (all(table(Y) == 1)) {
    row.names(X) <- as.character(Y)
    X
  } else {
    if (any(is.na(X))) {
      xspl <- split_matrix(X, Y)
      ret <- do.call(rbind, lapply(xspl, function(x) matrixStats::colMeans2(x, na.rm=TRUE)))
      row.names(ret) <- names(xspl)
      ret
    } else { 
      Rs <- rowsum(X,Y,na.rm=TRUE)
      yt <- table(Y)
      ret <- sweep(Rs, 1, yt, "/")
      row.names(ret) <- names(yt)
      ret
    }
  }
}


#' @export
split_matrix <- function(X, fac) {
  idx <- split(1:nrow(X), fac)
  lapply(idx, function(i) X[i,])
}



#' @keywords internal
create_folds <- function (y, k = 10, list = TRUE, returnTrain = FALSE) {
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2) 
      cuts <- 2
    if (cuts > 5) 
      cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0) 
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                                                               size = numInClass[i])
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                        sep = "")
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
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

colCors <- function(x, y) { 
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


code_replications <- function(f) {
  m <- outer(f, unique(f), "==")
  ordered(apply(m* apply(m,2,cumsum), 1, sum))
}


normalize <- function(x) {
  x/norm(x, "F")
}


#' @export
apply_scaling <- function(Xc) {
  cen <- if (is.null(attr(Xc, "scaled:center"))) FALSE else attr(Xc, "scaled:center")
  stdev <- if (is.null(attr(Xc, "scaled:scale"))) FALSE else attr(Xc, "scaled:scale")
  
  function(M) {
    if (!cen && !stdev) {
      M
    } else if (cen && !stdev) {
      sweep(M, 2, cen, "-")
    } else if (!cen && stdev) {
      sweep(M, 2, stdev, "/")
    } else if (cen && stdev) {
      M <- sweep(M, 2, cen)
      sweep(M, 2, stdev, "/")
    } else {
      stop()
    }
  }
}


combinedACC <- function(Pred, Obs) {
  levs <- levels(as.factor(Obs))
  maxind <- apply(Pred, 1, which.max)
  pclass <- levs[maxind]
  sum(pclass == Obs)/length(pclass)
  
}

combinedAUC <- function(Pred, Obs) {
  mean(sapply(1:ncol(Pred), function(i) {
    lev <- levels(Obs)[i]
    pos <- Obs == lev
    pclass <- Pred[,i]
    pother <- rowMeans(Pred[,-i,drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  }))
}





