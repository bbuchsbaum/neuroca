#' @export
ncomp.projector <- function(x) {
  ncol(x)
}

#' @export
nrow.projector <- function(x) {
  nrow(scores(x))
}

dim.projector <- function(x) {
  c(nrow(x), ncol(x))
}


#' @export
as.projector.matrix <- function(x,...) {
  identity_projector(x)
}

#' @export
projection_fun.identity_projector <- function(x) {
  function(x) x
}

#' @export
projection_fun.matrix <- function(x) {
  function(x) x
}

#' @export
identity_projector <- function(mat) {
  assert_that(is.matrix(mat))
  class(mat) <- c("identity_projector", "projector", "matrix")
  mat
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
  x
}


#' @export
#' 
ncomp.matrix <- function(x) ncol(x)


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

#' @export
project.block_matrix <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    as.matrix(x)
  } else {
    as.matrix(newdata)
  }
}


#' @export
ncol.identity_projector <- function(x) {
  ncol(x)
}

#' @export
ncomp.identity_projector <- function(x) {
  ncol(x)
}

#' @export
reprocess.identity_projector <- function(x) {
  x
}

#' @export
predict.projector <- function(x, newdata, ncomp=ncomp(x)) {
  project(x, newdata, comp=1:ncomp)
}

