


#' group_means
#' 
#' @param Y \code{factor} variable defining the groups
#' @param X \code{matrix} defining the matrix data to be group-wise averaged
#' @export
group_means <- function(Y, X) {
  
  if (all(table(Y) == 1)) {
    row.names(X) <- names(table(Y))
    X
  } else {
    Rs <- rowsum(X,Y)
    yt <- table(Y)
    ret <- sweep(Rs, 1, yt, "/")
    row.names(ret) <- names(yt)
    ret
  }
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


pre_process <- function(X, center=TRUE, scale=FALSE) {
  Xc <- scale(X, center=center, scale=scale)
  applyFun <- apply_scaling(Xc)
  attr(Xc, "applyFun") <- applyFun
  Xc
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



#group_means1 <- function(Y, X) {
# G <- model.matrix(~ Y - 1)
# colnames(G) <- levels(Y)
# GW <- G/colSums(G)
# R <- t(crossprod(X, GW))
#centroid <-  rowMeans(R)
#Rcent <- sweep(R, 1, centroid)
#list(Rcent=Rcent, Ymat=G, centroid=centroid)
#}\\



