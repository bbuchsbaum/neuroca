#' @export
ncomp.projector <- function(x) {
  x$ncomp
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
  function(newdata, colind=NULL) {
    newdata <- as.matrix(newdata)
    if (is.null(colind)) {
      assert_that(ncol(newdata) == ncol(x))
      newdata
    } else {
      assert_that(max(colind) <= ncol(x))
      assert_that(ncol(newdata) == length(colind))
      newdata
      #i <- rep(1:nrow(x), length(colind))
      #j <- rep(colind, each=nrow(x))
      #sparseMatrix(i=i,j=j, x=as.vector(newdata), dims=dim(x))
    }
  }
}

#' @export
projection_fun.Matrix <- function(x) {
  function(newdata, colind=NULL) {
    if (is.null(colind)) {
      assert_that(ncol(newdata) == ncol(x))
      newdata
    } else {
      assert_that(max(colind) <= ncol(x))
      assert_that(ncol(newdata) == length(colind))
      #i <- rep(1:nrow(x), length(colind))
      #j <- rep(colind, each=nrow(x))
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
project.matrix <- function(x, newdata=NULL, colind=NULL) {
  if (is.null(newdata)) {
    if (is.null(colind)) {
      x
    } else {
      assert_that(max(colind) <= ncol(x))
      assert_that(all(colind > 0))
      #i <- rep(1:nrow(x), length(colind))
      #j <- rep(colind, each=nrow(x))
      #sparseMatrix(i=i,j=j, x=as.vector(x[,colind]), dims=dim(x))
      x[,colind]
    }
  } else {
    if (is.null(colind)) {
      assert_that(dim(newdata) == dim(x))
      as.matrix(newdata)
    } else {
      assert_that(max(colind) <= ncol(x))
      assert_that(ncol(newdata) == length(colind))
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

compose_all <- function(...) {
  args <- list(...)
  assert_that(length(args) > 1)
  
  out <- compose.projector(args[[1]], args[[2]])
  if (length(args) > 2) {
    for (i in 3:length(args)) {
      out <- compose.projector(out, args[[i]])
    }
  }
  
  out
}




#' @export 
compose.projector <- function(x,y) {
  assert_that(inherits(y, "projector"), msg=paste("y does not inherit from class 'projector': ", class(y)))
  
  proj <- function(newdata) {
    project(y, project(x, newdata))
  }
  recon <- function(newdata=NULL) {
    reconstruct(x, reconstruct(y, newdata))
  }
  
  structure(
    list(
      a=x,
      b=y,
      d=y$d,
      u=y$u,
      xdim=dim(x),
      ncomp=ncomp(y),   
      proj=proj,
      scores=scores(y),
      recon=recon),
    class=c("composed_projector", "projector")
  )
}

#' @export
loadings.composed_projector <- function(x) {
  t(x$u %*% diag(1/x$d)) %*% reconstruct(x)
}

#' @export
scores.composed_projector <- function(x) { x$scores }

#' @export
dim.composed_projector <- function(x) c(nrow(x), ncol(x))

#' @export
ncol.composed_projector <- function(x) x$xdim[2]

#' @export
nrow.composed_projector <- function(x) x$xdim[1]

#' @export
ncomp.composed_projector <- function(x) ncol(scores(x))

#' @export
project.composed_projector <- function(x, newdata=NULL) {
  if (is.null(newdata)) {
    scores(x)
  } else {
    x$proj(newdata)
  }
}

#' @export
reconstruct.composed_projector <- function(x, newdata=NULL) {
  x$recon(newdata)
}

#' @export
combine.projector <- function(x,..., orthogonalize=FALSE, nretain=NULL) {
  xl <- c(list(x), list(...))
  
  scores <- do.call(cbind, map(xl, ~ scores(.)))
  
  if (orthogonalize || !is.null(nretain)) {
    if (is.null(nretain)) {
      nretain <- ncol(scores)
    }
    
    res <- RSpectra::svds(scores, k=nretain)
    pseudo_svd(res$u,res$v,res$d)
  } else {
    v <- do.call(cbind, map(xl, ~ loadings(.)))
    d <- apply(scores, 2, function(y) sqrt(sum(y^2)))

    u <- scores %*% Matrix::Diagonal(x=1/d)
    pseudo_svd(u,v,d)
  }
}



