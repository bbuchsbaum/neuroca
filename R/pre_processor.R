
.uncenter_scale <- function(newdata, center_vec, scale_vec, center, scale, subind=NULL) {
  if (is.null(subind)) {
    subind <- 1:ncol(newdata)
  }
  assert_that(length(subind) == ncol(newdata))
  
  if (!center && !scale) {
    newdata
  } else if (center && scale) {
    m1 <- sweep(newdata, 2, scale_vec[subind], "*")
    sweep(m1, 2, center_vec[subind], "+")
  } else if (!center && scale) {
    sweep(newdata, 2, scale_vec[subind], "*")
  } else if (center && !scale) {
    sweep(newdata, 2, center_vec[subind], "+")
  } else {
    stop()
  }
}

.center_scale <- function(newdata, center_vec, scale_vec, center, scale, subind=NULL) {
    if (is.null(subind)) {
      subind <- 1:ncol(newdata)
    }
    assert_that(length(subind) == ncol(newdata))
    
    if (!center && !scale) {
      newdata
    } else if (center && scale) {
      m1 <- sweep(newdata, 2, center_vec[subind], "-")
      sweep(m1, 2, scale_vec[subind], "/")
    } else if (!center && scale) {
      sweep(newdata, 2, scale_vec[subind], "/")
    } else if (center && !scale) {
      sweep(newdata, 2, center_vec[subind], "-")
    } else {
      stop()
    }
}


pre_process.matrix_pre_processor <- function(x, newdata, subind=1:length(x$center_vec)) {
  .center_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, subind)  
}

reverse_pre_process.matrix_pre_processor <- function(x, newdata, subind=1:length(x$center_vec)) {
  .uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, subind)  
}

pre_process.projector_pre_processor <- function(x, newdata, subind=NULL) {
  .center_scale(x$projfun(newdata, subind=subind), x$center_vec, x$scale_vec, x$center, x$scale, subind)  
}

reverse_pre_process.projector_pre_processor <- function(x, newdata, subind=NULL) {
  .uncenter_scale(x$projfun(newdata, subind=subind), x$center_vec, x$scale_vec, x$center, x$scale, subind)  
}


#' @importFrom matrixStats colSds
#' @export
pre_processor.matrix <- function(X, center=TRUE, scale=FALSE) {
  center_vec <- if (center) colMeans(X) else rep(0, ncol(X))
  scale_vec <- if (scale) matrixStats::colSds(X) else rep(1, ncol(X))

  structure(
    list(center=center,
    scale=scale,
    center_vec=center_vec,
    scale_vec=scale_vec),
    class="matrix_pre_processor")
  
}
  

#' @export
pre_processor.projector <- function(X, center=TRUE, scale=FALSE) {
  projfun <- projection_fun(X)
  Xp <- project(X)
  center_vec <- if (center) colMeans(Xp) else rep(0, ncol(Xp))
  scale_vec <- if (scale) matrixStats::colSds(Xp) else rep(1, ncol(Xp))
  
  structure(
    list(center=center,
         scale=scale,
         center_vec=center_vec,
         scale_vec=scale_vec,
         projfun=projfun),
    class="projector_pre_processor")
  
}


