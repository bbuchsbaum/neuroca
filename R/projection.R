#' @export
identity_projector <- function(mat) {
  assert_that(is.matrix(mat))
  class(mat) <- c("identity_projector", "projector", "matrix")
  mat
}

#' @export
loadings.identity_projector <- function(x) {
  rep(1:ncol(x))/ncol(x)
}

#' @export
ncomp.projector <- function(x) {
  ncol(x)
}

#' @export
project.identity_projector <- function(x, newdata=NULL, comp=1:ncomp(x), pre_process=FALSE, subind=NULL) {
  if (is.null(newdata) && is.null(subind)) {
    return(x)
  } else if (is.null(newdata) && !is.null(subind)) {
    return(x[,subind,drop=FALSE])
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata,nrow=1)
  }
  
  if (is.null(subind)) {
    newdata
  } else {
    assertthat::assert_that(length(subind) == ncol(newdata))
    newdata[,subind, drop=FALSE]
  }
}


#' @export
scores.identity_projector <- function(x) {
  rep(1:nrow(x))/nrow(x)
}

#' @export
project.matrix <- function(x, newdata=NULL, subind=NULL) {
  if (is.null(newdata)) {
    if (is.null(subind)) {
      x
    } else {
      x[,subind,drop=FALSE]
    }
  } else {
    if (is.null(subind)) {
      as.matrix(newdata)
    } else {
      as.matrix(newdata)[,subind,drop=FALSE]
    }
  }
}
