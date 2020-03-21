
gen_basis <- function(n,p) {
  matrix(rnorm(n*p),n,p)
}

#' @export
regress <- function(basis, X, preproc=pass(), method=c("linear", "ridge"), intercept=FALSE) {
  method <- match.arg(method)
  
  procres <- prep(preproc, X)
  Xp <- procres$init(X)
  
  # X ~ basis*b
  # X * b_inv = basis
  
  
  if (method == "linear") {
    lfit = lsfit(basis, X, intercept=intercept)
    
    if (intercept) {
      scores <- cbind(rep(1, nrow(basis)), basis)
    } else {
      scores <- basis
    }
    
    loadings <- t(coef(lfit))
    
  } else {
    stop("method 'ridge' not implemented yet.")
  }
  
  projector(procres, ncomp=ncol(scores), 
            basis=basis,
            v=loadings, 
            projection=corpcor::pseudoinverse(loadings),
            classes="regress")
  
}

#' @export
scores.regress <- function(x) {
  x$basis
}

#' @export
project.regress <- function(x, newdata, colind=NULL) {
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(colind)) {
    ## pre_process new data and project
    reprocess(x, newdata) %*% t(x$projection)
  } else {
    
    ## colind must equal number of columns of newdata
    ## colind cannot have more columns that original dataset
    assertthat::assert_that(length(colind) == ncol(newdata))
    #browser()
    reprocess(x, newdata, colind=colind) %*% t(x$projection[, colind,drop=FALSE])
  }

}

### TODO create 'genreconstruct' function

#' @export
reconstruct.regress <- function(x, newdata=NULL, colind=NULL, 
                                rowind=NULL, reverse_pre_process=TRUE) {
  if (is.null(newdata)) {
    newdata <- x$basis
  }
  
  rowind <- 1:nrow(newdata)
  genreconstruct(x,newdata,comp=1:nrow(loadings(x)),colind,rowind,reverse_pre_process)
  
}

residuals.regress <- function(x, ncomp) {
  Xr <- reconstruct(x)
}
