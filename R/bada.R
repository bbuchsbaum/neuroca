#' Barycentric Discriminant Analysis
#' 
#' A component technique that maximizes the between group variance over a set of variables. 
#' 
#' @importFrom assertthat assert_that 
#' @param Y dependent \code{factor} variable. The categories of the observations.
#' @param X the data matrix with \code{nrow} equal to \code{length} of \code{Y}
#' @param S an integer \code{vector} or \code{factor} of subject ids which is the same length as \code{Y}.
#' @param ncomp number of components to estimate
#' @param preproc pre-processing function, defaults to \code{center}
#' @param A the column constraints
#' @param M the row constraints
#' @param ... arguments to pass through
#' 
#' @details 
#' 
#' The \code{S} argument can be used to model multi-level structure. If \code{S} is included, then pre-processing is applied 
#' separately to each unique value of S. This has the effect, for example when preproc = \code{center} of removing subject-specific 
#' means before computing the barycetners.
#' 
#' @references
#' Abdi, H., Williams, L. J., & Bera, M. (2017). Barycentric discriminant analysis. \emph{Encyclopedia of Social Network Analysis and Mining}, 1-20.
#' 
#' @examples 
#' 
#' X <- matrix(rnorm(100*1000), 100, 1000)
#' Y <- factor(rep(letters[1:4], length.out=100))
#' S <- factor(rep(1:10, each=10))
#' 
#' bres <- bada(Y, X, S, ncomp=3)
#' 
#' project(bres, X[S==1,])
#' 
#' ## no strata
#' bres <- bada(Y, X, ncomp=3)
#' 
#' 
#' 
#' xbar <- matrix(rnorm(4*1000), 4, 1000)
#' Y <- factor(rep(letters[1:4], length.out=100))
#' X <- do.call(rbind, lapply(as.character(Y), function(y) { switch(y,
#'                                "a" = xbar[1,] + rnorm(1000)*5,
#'                                "b" = xbar[2,] + rnorm(1000)*5,
#'                                "c" = xbar[3,] + rnorm(1000)*5,
#'                                "d" = xbar[4,] + rnorm(1000)*5)
#'                                }
#'                                ))
#'                                
#' bres2 <- bada(Y, X, S, ncomp=3)
#' 
#' 
#' @export
bada <- function(Y, X, S=rep(1, nrow(X)), ncomp=length(levels(as.factor(Y)))-1, preproc=center(), A=NULL, M=NULL, ...) {
  assert_that(is.factor(Y))
  assert_that(length(Y) == nrow(X)) 
  assert_that(length(S) == nrow(X))
  
  if (!is.null(A)) {
    assert_that( (ncol(A) == nrow(A)) && (ncol(A) == ncol(X)))
  }
  
  if (!is.null(M)) {
    assert_that( (ncol(M) == nrow(M)) && (nrow(M) == nrow(X)))
  }
  
  S <- factor(S)
  
  procres <- if (length(levels(S)) > 1) {
    Xs <- lapply(levels(S), function(lev) {
      Xs <- X[S == lev,,drop=FALSE]
    })
    
    procres <- lapply(levels(S), function(lev) {
      p <- fresh(preproc)
      prep(p)
    })
    names(procres) <- levels(S)
    Xp <- do.call(rbind, lapply(1:length(procres), function(i) procres[[i]]$init(Xs[[i]])))
    procres
  } else {
    p <- prep(preproc)
    Xp <- p$init(X)
    p
  }
  
  
  ncomp <- min(ncomp, length(levels(Y)))
  
  Xr <- group_means(Y, Xp)
  
  fit <- if (!is.null(A) || !is.null(M) ) {
    genpca(Xr, ncomp=ncomp, preproc=pass(), A=A, M=M)
  } else {
    pca(Xr, ncomp=ncomp, preproc=pass())
  }
  
  ret <- list(
    preproc=preproc,
    procres=procres,
    ncomp=fit$ncomp,
    fit=fit,
    Y=Y,
    X=X,
    Xr=Xr,
    S=S)
  
  class(ret) <- c("bada")
  ret
}

#' @export
rotate.bada <- function(x, rot) {
  rfit <- rotate(x$fit, rot)
  ret <- list(
    preproc=x$preproc,
    procres=x$procres,
    fit=rfit,
    Y=x$Y,
    X=x$X,
    Xr=x$Xr,
    S=x$S)
  class(ret) <- c("bada")
  ret
}

#' @export
project.bada <- function(x, newdata=NULL, comp=1:x$fit$ncomp, colind=NULL, stratum=NULL) {
  if (is.null(newdata)) {
    project(x$fit, newdata=x$Xr, comp=comp, colind=colind)
  } else {
    ## stratum-specific pre-processing
    Xp <- reprocess(x, newdata, colind=colind, stratum=stratum)
    project(x$fit, Xp, comp=comp, colind=colind)
  }
}

## project from existing table
#' @export
predict.bada <- function(x, newdata, ncomp=x$ncomp, colind=NULL, stratum=NULL, 
                         type=c("class", "prob", "scores", "crossprod", "distance", "cosine")) {
  
  
  type <- match.arg(type)
  Xp <- reprocess(x, newdata, colind=colind, stratum=stratum)
  
  fscores <- project(x$fit, Xp, comp=1:ncomp, colind=colind)
  
  if (type == "scores") {
    fscores
  } else {
    scorepred(as.matrix(fscores), scores(x), type=type, ncomp=ncomp)
  }
  
}


#' @export
project_cols.bada <- function(x, newdata=NULL, comp=1:x$ncomp) {
  project_cols(x$fit,newdata,comp=comp)
}

#' @export
singular_values.bada <- function(x) x$fit$d

#' @export
loadings.bada <- function(x) loadings(x$fit)

#' @export
scores.bada <- function(x) {
  sc <- scores(x$fit)
  row.names(sc) <- levels(x$Y)
  sc
}


#' @export
reprocess.bada <- function(x,
                           newdata,
                           colind = NULL,
                           stratum = NULL) {
  

  if (is.null(colind)) {
    assert_that(ncol(newdata) == nrow(loadings(x)))
  }
  
  avg_preproc <- function(newdata, colind = NULL) {
    Reduce("+", lapply(1:length(levels(x$S)), function(i) {
      p <- x$procres[[i]]
      p$transform(newdata, colind)
    })) / length(levels(x$S))
  }
  
  stratum_preproc <- function(newdata, stratum, colind = NULL) {
    assert_that(stratum %in% levels(x$S))
    p <- x$procres[[stratum]]
    p$transform(newdata, colind)
  }
  
  
  ns <- length(unique(x$S))
  
  
  
  Xp <- if (ns > 1 && is.null(stratum)) {
    ## we have multiple strata
    avg_preproc(newdata, colind)
  } else if (ns > 1 && !is.null(stratum)) {
    assert_that(stratum %in% levels(x$S))
    stratum_preproc(newdata, stratum, colind)
  } else {
    ## no strata
    x$procres$transform(newdata, colind)
  }
}


#' @export
permute.bada <- function(x) {
  split_idx <- split(1:length(x$Y), x$S)
  Yperm <- unlist(lapply(split_idx, function(idx) {
    sample(x$Y[idx])
  }))
  
  Yperm <- Yperm[unlist(split_idx)]
}

#' @export
reconstruct.bada <- function(x, ncomp=x$ncomp) {
  recon <- reconstruct(x$fit,comp=1:ncomp)
  row.names(recon) <- levels(x$Y)
  do.call(rbind, lapply(levels(x$S), function(s) {
    p <- x$procres[[s]]
    out <- p$reverse_transform(recon)
    y <- x$Y[x$S == s]
    out <- recon[match(y, row.names(recon)),]
    row.names(out) <- paste0("S_", s, "_", y)
    out
  }))
  
}

#' @export
residuals.bada <- function(x, ncomp=x$ncomp) {
  recon <- reconstruct(x,ncomp)
  recon - x$X
}


#' @export
refit.bada <- function(x, Y, Xlist, S = x$S, ncomp=x$ncomp,...) { 
  ##bada(Y, Xlist, ncomp=ncomp, x$preproc, normalization=x$normalization, A=x$fit$A, M=x$fit$M) 
  bada(Y, Xlist, ncomp=ncomp, x$preproc, ...) 
}




#' @export
print.bada <- function(object) {
  showk <- 1:min(object$ncomp, 5)

  cat("barycentric discriminant analysis (BADA) decomposition", "\n")
  cat("  number of components: ", object$ncomp, "\n")
  cat("  number of variables: ", nrow(object$fit$v), "\n")
  cat("  number of categories: ", nrow(object$fit$u), "\n")
  
  if (length(levels(object$Y)) > 8) {
    cat("  category levels: ", paste(levels(object$Y)[1:8], collapse = " "), " ... \n")
  } else {
    cat("  category levels: ", paste(levels(object$Y), collapse = " "), "\n")
  }
  cat("  number of original observations: ", nrow(object$X), "\n")
  cat("  center: ", object$fit$preproc$center, "\n")
  cat("  scale: ", object$fit$preproc$scale, "\n")
  cat("  % variance explained (top ", showk, "): ", ((object$fit$d[showk]^2)/sum(object$fit$d^2)) * 100, "\n")
}




