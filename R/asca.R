

gen_threeway_mat <- function(n1,n2,n3, nvox, index) {
  nk <- n1*n2*n3

  grid <- expand.grid(A=paste0("A", 1:n1), B=paste0("B", 1:n2), C=paste0("C", 1:n3))
  grid$sid <- index
  mat <- matrix(rnorm(nk*nvox), nk, nvox)
  list(mat=mat, design=grid)

}

Xs <- lapply(1:10, function(i) gen_threeway_mat(3,3,5, 10, i))
Xlist <- lapply(Xs, "[[", "mat")
design <- Xs[[1]]$design
formula = ~ A + B + C + A:B + A:C + B:C + A:B:C
strata <- "sid"

#' musu_asca
#' 
#' 
#' @param Xlist the list of data matrices
#' @param formula a formula specifying the ANOVA design
#' @param design a \code{data.frame} providing the variables provided in \code{formula} argument.
#' @param center
#' @param scale
#' @param svd.method the svd method to use.
#' @param strata name of variable in design matrix that indicates any blocking structure, such as "subject".
#' @export
musu_asca <- function(Xlist, formula, design, center=TRUE, scale=FALSE, svd.method="base", strata) {
  tform <- terms(formula)
  facs <- attr(tform, "factors")
  
  termorder <- apply(facs,2, sum)
  orders <- seq(1, max(termorder))
  
  if (!is.null(strata)) {
    strata_fac <- as.factor(design[[strata]])
    design[[strata]] <- as.factor(design[[strata]])
    assert_that(length(strata_fac) == nrow(design))
  }
  
  blockInd <- blockIndices(Xlist)
  
  X <- do.call(cbind, Xlist)
  
  get_lower_form <- function(ord) {
    tnames <- colnames(facs)[termorder < ord]
    
    meat <- if (!is.null(strata)) {
      paste0("(", paste(tnames, collapse= " + "), ")*", strata)
    } else {
      paste(tnames, collapse= " + ")
    }
    
    as.formula(paste(" ~ ", meat))
  }
  
  expand_scores <- function(scores, y) {
    mind <- match(as.character(y), row.names(scores))
    scores[mind,]
  }
  
  res <- lapply(1:length(termorder), function(i) {
    ord <- termorder[i]
    message("pca decomposition for factor ", colnames(facs)[i])
    if (ord == 1) {
      form <- as.formula(paste("~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form,data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))], levels=cnames)
      bres <- musubada(Y, Xlist, ncomp=length(levels(Y)), center=TRUE, svd.method=svd.method, normalization="None")
      list(G=G, Y=Y, Glower=NULL, form=form, lower_form = ~ 1, term=colnames(facs)[i], ex_scores=expand_scores(bres$scores, Y), scores=bres$scores, bada_result=bres)
    } else {
      lower_form <- get_lower_form(ord)
      Glower <- model.matrix(lower_form, data=design)
      Xresid <- resid(lsfit(Glower, X, intercept=FALSE))
      Xresidlist <- lapply(1:nrow(blockInd), function(i) Xresid[, blockInd[i,1]:blockInd[i,2]])
      
      form <- as.formula(paste("~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form, data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))])
      bres <- musubada(Y, Xresidlist, ncomp=length(levels(Y)), center=TRUE, svd.method=svd.method, normalization="None")
      list(G=G, Y=Y, Glower, form=form, lower_form = lower_form, term=colnames(facs)[i], ex_scores=expand_scores(bres$scores, Y), scores=bres$scores, bada_result=bres)
    }
  })
  
  permute <- function(obj) {
    idx <- sample(1:nrow(obj$X))
    list(X=obj$X, Y=obj$Y[idx,], idx=idx, group=group[idx])
  }
  
  expanded_scores <- lapply(res, "[[", "ex_scores")
  scores <- lapply(res, "[[", "scores")
  names(scores) <- colnames(facs)
  names(expanded_scores) <- colnames(facs)
  names(res) <- colnames(facs)
  
  ret <- list(
    results=res,
    scores=scores,
    expanded_scores=expanded_scores,
    formula=formula,
    strata=strata,
    design=design,
    terms=colnames(facs)
  )
  
  
  class(ret) <- c("projector", "musu_asca")
  ret
  
}





#' 
#' @param X the data matrix
#' @param formula a formula specifying the ANOVA design
#' @param design a \code{data.frame} providing the variables provided in \code{formula} argument.
#' @param center
#' @param scale
#' @param svd.method the svd method to use.
#' @param strata name of variable in design matrix that indicates any blocking structure, such as "subject".
asca <- function(X, formula, design, center=TRUE, scale=FALSE, svd.method="base", strata=NULL) {
  tform <- terms(formula)
  facs <- attr(tform, "factors")
  
  termorder <- apply(facs,2, sum)
  orders <- seq(1, max(termorder))
  
  if (!is.null(strata)) {
    strata_fac <- as.factor(design[[strata]])
    design[[strata]] <- as.factor(design[[strata]])
    assert_that(length(strata_fac) == nrow(design))
  }
  
  get_lower_form <- function(ord) {
    tnames <- colnames(facs)[termorder < ord]
    
    meat <- if (!is.null(strata)) {
      paste0("(", paste(tnames, collapse= " + "), ")*", strata)
    } else {
      paste(tnames, collapse= " + ")
    }
    
    as.formula(paste("X ~ ", meat))
  }
  
  expand_scores <- function(scores, y) {
    mind <- match(as.character(y), row.names(scores))
    scores[mind,]
  }
  
  res <- lapply(1:length(termorder), function(i) {
    ord <- termorder[i]
    message("pca decomposition for factor ", colnames(facs)[i])
    if (ord == 1) {
      form <- as.formula(paste("X ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form,data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))], levels=cnames)
      bres <- bada(Y, X, ncomp=length(levels(Y)), center=TRUE, svd.method=svd.method, strata=strata_fac)
      list(G=G, Y=Y, Glower=NULL, form=form, lower_form = ~ 1, term=colnames(facs)[i], ex_scores=expand_scores(bres$scores, Y), scores=bres$scores, bada_result=bres)
    } else {
      lower_form <- get_lower_form(ord)
      Glower <- model.matrix(lower_form, data=design)
      Xresid <- resid(lsfit(Glower, X, intercept=FALSE))
      form <- as.formula(paste("Xresid ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form, data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))])
      bres <- bada(Y, Xresid, ncomp=length(levels(Y)), center=TRUE, svd.method=svd.method, strata=strata_fac)
      list(G=G, Y=Y, Glower, form=form, lower_form = lower_form, term=colnames(facs)[i], ex_scores=expand_scores(bres$scores, Y), scores=bres$scores, bada_result=bres)
    }
  })
  
  permute <- function(obj) {
    idx <- sample(1:nrow(obj$X))
    list(X=obj$X, Y=obj$Y[idx,], idx=idx, group=group[idx])
  }
  
  expanded_scores <- lapply(res, "[[", "ex_scores")
  scores <- lapply(res, "[[", "scores")
  names(scores) <- colnames(facs)
  names(expanded_scores) <- colnames(facs)
  names(res) <- colnames(facs)
  
  ret <- list(
    results=res,
    scores=scores,
    expanded_scores=expanded_scores,
    formula=formula,
    strata=strata,
    design=design,
    terms=colnames(facs)
  )
  
  
  class(ret) <- c("projector", "asca")
  ret
  
}