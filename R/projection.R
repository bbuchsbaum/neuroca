


#' @export
scores.projector <- function(x) {
  if (is.null(x$scores)) {
    project(x, x$procres$Xp)
  } else {
    x$scores
  }
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
  assert_that(inherits(y, "projector"), 
              msg=paste("y does not inherit from class 'projector': ", class(y)))
  
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



