#' plscorr_aov
#' 
#' @export
#' @param Xlist the data matrices, 1 per subject
#' @param formula a formula specifying the ANOVA design
#' @param design a \code{data.frame} providing the variables provided in \code{formula} argument.
#' @param svd.method the svd method to use.
musu_asca <- function(Xlist, formula, design, svd.method="base") {
  tform <- terms(formula)
  facs <- attr(tform, "factors")
  
  termorder <- apply(facs,2, sum)
  orders <- seq(1, max(termorder))
  
  X <- do.call(cbind, Xlist)
  u <- colMeans(X)
  
  blockIndices <- rep(1:length(Xlist), sapply(Xlist, ncol))
  blocks <- sapply(Xlist, ncol)

  get_lower_form <- function(ord) {
    tnames <- colnames(facs)[termorder < ord]
    
    meat <- if (!is.null(random)) {
      paste0("(", paste(tnames, collapse= " + "), ")*", random)
    } else {
      paste(tnames, collapse= " + ")
    }
    
    as.formula(paste("X ~ ", meat))
  }
  
  res <- lapply(1:length(termorder), function(i) {
    ord <- termorder[i]
    if (ord == 1) {
      form <- as.formula(paste("X ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form,data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))], levels=cnames)
      
      #X0 <- t(t(X) %*% G)
      #X0c <- scale(X0, center=TRUE, scale=FALSE)
    
      mbres <- musubada(Y, Xlist, ncomp=length(levels(Y)), center=TRUE, svd.method=svd.method, normalization="None")
      list(G=G, Y=Y, Glower=NULL, form=form, lower_form = ~ 1, plsres=plsres)
    } else {
      lower_form <- get_lower_form(ord)
      Glower <- model.matrix(lower_form, data=design)
      Xresid <- resid(lsfit(Glower, X, intercept=FALSE))
      form <- as.formula(paste("Xresid ~ ", colnames(facs)[i], "-1"))
      G <- model.matrix(form, data=design)
      cnames <- colnames(G)
      Y <- factor(cnames[apply(G, 1, function(x) which(x==1))])
      plsres <- plscorr_contrast(Xresid, G, strata=strata, center=center, scale=scale, ncomp=ncomp, svd.method=svd.method)
      list(G=G, Glower, form=form, lower_form = lower_form, plsres=plsres)
    }
  })
  
  permute <- function(obj) {
    idx <- sample(1:nrow(obj$X))
    list(X=obj$X, Y=obj$Y[idx,], idx=idx, group=group[idx])
  }
  
  names(res) <- colnames(facs)
  ret <- list(
    results=res,
    formula=formula,
    strata=strata,
    design=design,
    terms=colnames(facs)
  )
  
  
  class(ret) <- "plscorr_result_aov"
  ret
  
}