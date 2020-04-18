



#' construct a `bi-projector` instance
#' @export
#' @inheritParams projector
#' @param preproc
#' @param ncomp
#' @param u
#' @param d
#' @param scores
bi_projector <- function(preproc, ncomp, v, u, d, scores, classes=NULL, ...) {
  out <- list(
    preproc=preproc,
    ncomp=ncomp,
    v=v,
    u=u,
    d=d,
    scores=scores,
    ...)
  
  class(out) <- c(classes, "bi_projector", "projector")
  out
}
    



#shrink.pca <- function(x, method=c("GSURE", "QUT", "SURE"), sigma, lambda=NULL, gamma=NULL) {
#  
#}

#' shrink_pca
#' 
#' adaptive shrinkage pca from the \code{denoiseR} package
#' 
#'   
#' @param X
#' @param center
#' @param scale
#' @importFrom denoiseR adashrink
#' @export
shrink_pca <- function(X, preproc=center(), method = c("GSURE", "QUT", "SURE"), ...) {
  assert_that(is.matrix(X) || inherits(X, "Matrix"))
  
  procres <- prep(preproc, X)
  Xp <- procres$init(X)
  
  res <- denoiseR::adashrink(Xp, method=method, center=FALSE, ...)
  
  keep <- res$singval > 1e-06
  
  if (sum(keep) == 0) {
    keep <- 1
  } 
  
  v=res$low.rank$v[,keep,drop=FALSE]
  u=res$low.rank$u[,keep,drop=FALSE]
  d=res$low.rank$d[keep]

  ret <- bi_projector(
              preproc=procres,
              ncomp=length(d),
              v=v, 
              u=u,
              d=d,
              scores=t(t(as.matrix(u)) * d),
              classes=c("shrink_pca", "pca"))

}
                       

#' pseudo_svd
#' 
#' @param u the row weights
#' @param v the column weights
#' @param d the singular values
#' @param rnames row names of observations
#' @export
pseudo_svd <- function(u, v, d, rnames=NULL) {
  
  scores <- t(t(as.matrix(u)) * d)
  
  if (!is.null(rnames)) {
    row.names(scores) <- rnames
  }
  
  assert_that(length(d) == ncol(u))
  assert_that(ncol(v) == ncol(u))
  

  ret <- bi_projector(
              preproc=pass(),
              ncomp=length(d),
              v=v, 
              u=u, 
              d=d, 
              scores=scores,
              classes="pseudo_svd")
  ret
}

#' principal components analysis
#' 
#' @param X
#' @param ncomp
#' @param preproc
#' @export
#' 
#' @examples 
#' 
#' X <- matrix(rnorm(10*20), 10, 20)
#' res <- pca(X, ncomp=10, preproc=center())
pca <- function(X, ncomp=min(dim(X)), preproc=center(), ...) {
  assert_that(is.matrix(X) || inherits(X, "Matrix"))
  
  procres <- prep(preproc)
  
  Xp <- procres$init(X)
  
  svdres <- svd_wrapper(Xp, ncomp, ...)
  
  scores <- t(t(as.matrix(svdres$u)) * svdres$d)
  
  if (!is.null(row.names(scores))) {
    row.names(scores) <- row.names(Xp)[seq_along(svdres$d)]
  }
  
  ret <- bi_projector(
              preproc=procres,
              length(svdres$d),
              v=svdres$v, 
              u=svdres$u, 
              d=svdres$d, 
              scores=scores,
              classes=c("pca"))
  ret
}

#' @export
refit.pca <- function(x, X) {
  pca(X, ncomp=x$ncomp, preproc=x$preproc$preproc)
}

#' @export
permutation.pca <- function(x, X, nperm=100) {
  Q <- min(dim(X)) - 1
  evals <- x$d^2
  Fa <- sapply(1:length(evals), function(i) evals[i]/sum(evals[i:Q]))
  
  F1_perm <- sapply(1:nperm, function(i) {
    Xperm <- apply(X, 2, function(x) sample(x))
    fit <- refit(x, Xperm)
    evals <- fit$d^2
    F1_perm <- evals[1]/sum(evals[1:length(evals)])
  })
  
  if (x$ncomp > 1) {
    Fq <- parallel::mclapply(2:x$ncomp, function(a) {
      sapply(1:nperm, function(j) {
        Ea <- residuals(x, a-1, X)
        Ea_perm <- apply(Ea, 2, function(x) sample(x))
      
        I <- diag(nrow(X))
      
        uu <- Reduce("+", lapply(1:(a-1), function(i) {
          x$u[,i,drop=FALSE] %*% t(x$u[,i,drop=FALSE]) 
        }))
      
        Ea_perm_proj <- (I - uu) %*% Ea_perm
  
        fit <- refit(x, Ea_perm_proj)
        evals <- fit$d^2
        Fq_perm <- evals[1]/sum(evals)
      })
    })
      
    Fq <- do.call(cbind, Fq)
    Fq <- cbind(F1_perm, Fq)
  } else {
    Fq <- as.matrix(F1_perm)
  }
  
  nperm=500
  
    
}
 
#' @export   
residuals.pca <- function(x, ncomp, xorig) {
  recon <- reconstruct(x, comp=1:ncomp)
  xorig - recon
}


#' @export
rotate.pca <- function(x, rot) {
  u_rot <- x$u %*% rot
  v_rot <- x$v %*% rot
  sc_rot <- scores(x) %*% rot
  
  ret <- bi_projector(
    preproc=x$preproc,
    ncomp=length(x$d),
    v=v_rot, 
    u=u_rot, 
    d=apply(sc_rot, 2, function(x) sqrt(sum(x^2))),
    scores=sc_rot, 
    classes=c("rotated_pca", "pca"))
  
  ret
  
}

genreprocess <- function(x, newdata, colind=NULL) {
  if (is.null(colind)) {
    assert_that(ncol(newdata) == nrow(loadings(x)))
    x$preproc$transform(newdata)
  } else {
    assert_that(length(colind) == ncol(newdata), 
                msg=paste("length of colind not equal to number of columns of newdata", length(colind), "!=", ncol(newdata)))
    x$preproc$transform(newdata, colind)
  }
}

#' @export
reprocess.projector <- function(x, newdata, colind=NULL) {
 genreprocess(x,newdata,colind)
  
}

#' TODO is this needed
#' @export
reprocess.bi_projector <- function(x, newdata, colind=NULL) {
  genreprocess(x,newdata,colind)
}

#' @export
singular_values.bi_projector <- function(x) {
  x$d
}

#' @export
loadings.bi_projector <- function(x) {
  x$v
}


#' @export
projection_fun.bi_projector <- function(x, comp=1:ncomp(x), pre_process=TRUE) {
  .comp <- comp
  .pre_process <- pre_process
  
  function(newdata, colind=NULL) {
    .colind <- colind
    project(x, newdata, comp=.comp, pre_process=.pre_process, colind=.colind)
  }
}

#' @export
project_cols.bi_projector <- function(x, newdata, comp=1:ncomp(x)) {

  ## if no newdata, then simply return the loadings
  if (missing(newdata)) {
    return(loadings(x)[,comp, drop=FALSE])
  } else {
    if (is.vector(newdata)) {
      newdata <- as.matrix(newdata, ncol=1)
    }
    assert_that(nrow(newdata) == nrow(x$u))
    t(newdata) %*% (x$u[,comp,drop=FALSE]) %*% diag(1/x$d[comp], nrow=length(comp), ncol=length(comp))
  }
}

#' @export
project.projector <- function(x, newdata, comp=1:ncomp(x), colind=NULL) {
  assert_that(max(comp) <= ncomp(x))
  ## if no newdata, then simply return the factor scores
  if (missing(newdata)) {
    ## TODO deal with colind
    return(scores(x)[,comp, drop=FALSE])
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(colind)) {
    ## pre_process new data and project
    reprocess(x, newdata) %*% loadings(x)[,comp, drop=FALSE]
  } else {
    ## colind must equal number of columns of newdata
    ## colind cannot have more columns that original dataset
    assertthat::assert_that(length(colind) == ncol(newdata))
    reprocess(x, newdata, colind=colind) %*% (loadings(x)[colind, comp, drop=FALSE])
  }
}

#' @export
residuals.bi_projector <- function(x, ncomp=1, xorig) {
  recon <- reconstruct(x,comp=1:ncomp)
  xorig - recon
}


genreconstruct <- function(x, newdata, comp, colind=NULL, 
                           rowind=NULL, reverse_pre_process=TRUE) {
  if (!is.null(newdata)) {
    assert_that(ncol(newdata) == length(comp) && nrow(newdata) == nrow(scores(x)))
  } else {
    newdata <- scores(x)[,comp, drop=FALSE]
  }
  
  if (is.null(rowind)) {
    rowind <- 1:nrow(newdata)
  } else {
    assert_that(all(rowind > 0) && max(rowind) < nrow(newdata))
  }
  
  if (is.null(colind)) {
    if (reverse_pre_process) {
      x$preproc$reverse_transform(newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE]))
    } else {
      newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE])
    }
  } else {
    if (reverse_pre_process) {
      x$preproc$reverse_transform(newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE])[,colind], 
                                  colind=colind)
    } else {
      newdata[rowind,,drop=FALSE] %*% t(loadings(x)[,comp,drop=FALSE])[,colind]
    }
  }
  
}

#' @export
reconstruct.bi_projector <- function(x, newdata=NULL, comp=1:x$ncomp, colind=NULL, rowind=NULL, reverse_pre_process=TRUE) {
  genreconstruct(x,newdata,comp,colind,rowind,reverse_pre_process)
}


#' @export
ncol.projector <- function(x) {
  nrow(x$v)
}

#' @export
nrow.bi_projector <- function(x) {
  nrow(x$u)
}

#' @export
ncomp.projector <- function(x) {
  x$ncomp
}

#' @export
scores.bi_projector <- function(x) {
  x$scores
}


#' @export
dim.projector <- function(x) {
  c(nrow(x$u), nrow(x$v))
}



#' @export
truncate.bi_projector <- function(obj, ncomp) {
  if (ncomp >= obj$ncomp) {
    warning("number of components to keep is greater than or equal to rank of pca fit, returning original model")
    ret <- obj
  } else {
    ret <- list(v=obj$v[,1:ncomp,drop=FALSE], 
                u=obj$u[,1:ncomp,drop=FALSE], 
                d=obj$d[1:ncomp], 
                scores=obj$scores[,1:ncomp,drop=FALSE], 
                ncomp=ncomp, 
                preproc=obj$preproc)
    
    class(ret) <- c("bi_projector", "list")
  }
  
  ret
}


#' @export
contributions.projector <- function(x) {
  loadings(x)^2
}


#' @export
contributions.bi_projector <- function(x, type=c("column", "row")) {
  type <- match.arg(type)
  if (type == "column") {
    loadings(x)^2
  } else {
    x$u^2
  }
}




#' @export
print.bi_projector <- function(object) {
  showk <- 1:min(object$ncomp, 5)
  cat(class(object)[1], "\n")
  cat("  number of components: ", object$ncomp, "\n")
  cat("  number of variables: ", nrow(object$v), "\n")
  cat("  number of observations: ", nrow(object$u), "\n")
  cat("  % variance explained (top 5): ", ((object$d[showk]^2)/sum(object$d^2)) * 100, "\n")
}


#' @importFrom explor explor
explor.pca <- function(obj) {
  if (!inherits(obj, "pca")) stop("obj must be of class pca")
  
  ## results preparation
  res <- prepare_results(obj)
  
  ## Settings
  settings <- list()
  settings$var_columns <- c("Variable", "Coord")
  settings$varsup_columns <- c("Variable", "Coord")
  settings$ind_columns <- c("Name", "Coord")
  settings$indsup_columns <- c("Name", "Coord")
  settings$scale_unit <- FALSE
  settings$obj_name <- deparse(substitute(obj))    
  
  settings$has_count <- FALSE
  settings$has_contrib <- FALSE
  settings$has_cos2 <- FALSE
  settings$has_var_eta2 <- FALSE
  settings$has_varsup_eta2 <- FALSE
  
  
  ## Launch interface
  explor:::explor_multi_pca(res, settings)
  
}
  

  
#' @importFrom explor prepare_results
prepare_results.pca <- function(obj) {
    
    if (!inherits(obj, "pca")) stop("obj must be of class pca")
    
    vars <- loadings(obj)
    colnames(vars) <- paste0("PC", 1:length(obj$d))
    vars <- data.frame(vars)
    ## Axes names and inertia
    axes <- seq_len(length(obj$d))
    percent <- round(obj$d^2 / sum(obj$d^2) *100, 2)
    names(axes) <- paste("Axis", axes, paste0("(", percent,"%)"))
    ## Eigenvalues
    eig <- data.frame(dim = 1:length(obj$d), percent = percent)
    
    ## Variables data coordinates
    vars$varname <- rownames(vars)
    vars$modname <- ""
    vars$Type <- "Active"
    vars$Class <- "Quantitative"
    
    vars <- vars %>% gather(Axis, Coord, starts_with("PC")) %>%
      mutate(Axis = gsub("PC", "", Axis, fixed = TRUE),
             Coord = round(Coord, 3))
    
    #vars <- vars %>% rename(Variable = varname, Level = modname)
    
    
    ## TODO have way of retrieivng original data?
    recon <- reconstruct(obj)
    #tmp <- data.frame(t(cor( scores(obj),obj$preproc$Xp)))
    tmp <- data.frame(t(cor( scores(obj),recon)))
    colnames(tmp) <- paste0("PC", 1:ncol(tmp))
    
    tmp$varname <- rownames(tmp)
    tmp$modname <- ""
    tmp$Type <- "Active"
    tmp$Class <- "Quantitative"  
   
    tmp <- tmp %>% gather(Axis, Cor, starts_with("PC")) %>%
      mutate(Axis = gsub("PC", "", Axis, fixed = TRUE),
             Cor = round(Cor, 3))
    
    vars <- vars %>% left_join(tmp, by = c("varname", "modname", "Type", "Class", "Axis"), na_matches = "na")
    
  
    
    vars$Contrib <- NA
    vars$Cos2 <- NA
    
    vars <- vars %>% rename(Variable = varname, Level = modname)
    
    
    ## Individuals coordinates
    ind <- data.frame(scores(obj))
    colnames(ind) <- paste0("PC", 1:ncol(ind))
    ind$Name <- rownames(ind)
    ind$Type <- "Active"
    
    ind <- ind %>% gather(Axis, Coord, starts_with("PC")) %>%
      mutate(Axis = gsub("PC", "", Axis, fixed = TRUE),
             Coord = round(Coord, 3))
    ind$Contrib <- NA
    ind$Cos2 <- NA
    
    ## cor var
    ## cor(pres2$preproc$Xp, scores(pres2))
    
    return(list(vars = vars, ind = ind, eig = eig, axes = axes))
    
}




