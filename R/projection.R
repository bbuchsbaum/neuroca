#' @export
ncomp.projector <- function(x) {
  ncol(x)
}

#' @export
nrow.projector <- function(x) {
  nrow(scores(x))
}

#' @export
dim.projector <- function(x) {
  c(nrow(x), ncol(x))
}

#' @export
is_transformer.projector <- function(x) TRUE

#' @export
is_transformer.matrix <- function(x) FALSE

#' @export
is_transformer.Matrix <- function(x) FALSE





#' @export
predict.projector <- function(x, newdata, ncomp=ncomp(x)) {
  project(x, newdata, comp=1:ncomp)
}


#' @export
projection_fun.matrix <- function(x) {
  function(newdata, subind=NULL) {
    newdata <- as.matrix(newdata)
    if (is.null(subind)) {
      assert_that(ncol(newdata) == ncol(x))
      newdata
    } else {
      assert_that(max(subind) <= ncol(x))
      assert_that(ncol(newdata) == length(subind))
      newdata
      #i <- rep(1:nrow(x), length(subind))
      #j <- rep(subind, each=nrow(x))
      #sparseMatrix(i=i,j=j, x=as.vector(newdata), dims=dim(x))
    }
  }
}

#' @export
projection_fun.Matrix <- function(x) {
  function(newdata, subind=NULL) {
    if (is.null(subind)) {
      assert_that(ncol(newdata) == ncol(x))
      newdata
    } else {
      assert_that(max(subind) <= ncol(x))
      assert_that(ncol(newdata) == length(subind))
      #i <- rep(1:nrow(x), length(subind))
      #j <- rep(subind, each=nrow(x))
      #sparseMatrix(i=i,j=j, x=as.vector(newdata), dims=dim(x))
      newdata
    }
  }
}



#' @export
#' 
ncomp.matrix <- function(x) ncol(x)

#' @export
#' 
ncomp.Matrix <- function(x) ncol(x)


#' @export
project.matrix <- function(x, newdata=NULL, subind=NULL) {
  if (is.null(newdata)) {
    if (is.null(subind)) {
      x
    } else {
      assert_that(max(subind) <= ncol(x))
      assert_that(all(subind > 0))
      #i <- rep(1:nrow(x), length(subind))
      #j <- rep(subind, each=nrow(x))
      #sparseMatrix(i=i,j=j, x=as.vector(x[,subind]), dims=dim(x))
      x[,subind]
    }
  } else {
    if (is.null(subind)) {
      assert_that(dim(newdata) == dim(x))
      as.matrix(newdata)
    } else {
      assert_that(max(subind) <= ncol(x))
      assert_that(ncol(newdata) == length(subind))
      newdata
    }
  }
}

#' @export
project.block_matrix <- function(x, newdata=NULL) {
  stop()
  #if (is.null(newdata)) {
  #  as.matrix(x)
  #} else {
  #  as.matrix(newdata)
  #}
}


