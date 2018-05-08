
.uncenter_scale <- function(newdata, center_vec, scale_vec, center, scale) {
  if (!center && !scale) {
    newdata
  } else if (center && scale) {
    m1 <- sweep(newdata, 2, scale_vec, "*")
    sweep(m1, 2, center_vec, "+")
  } else if (!center && scale) {
    sweep(newdata, 2, scale_vec, "*")
  } else if (center && !scale) {
    sweep(newdata, 2, center_vec, "+")
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

#' @export
pre_process.matrix_pre_processor <- function(x, newdata=NULL, subind=1:length(x$center_vec)) {
  if (is.null(newdata)) {
    x$Xp[,subind]
  } else {
    .center_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, subind)  
  }
}

#' @export
reverse_pre_process.matrix_pre_processor <- function(x, newdata) {
  .uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale)  
}

#' @export
pre_process.block_matrix_pre_processor <- function(x, newdata=NULL, block_index=NULL) {
  assert_that(is.null(block_index) || length(block_index) == 1)

  if (is.null(newdata)) {
    if (!is.null(block_index)) {
      subind <- block_index_list(x)[[block_index]]
      .center_scale(x$Xp[,subind], x$center_vec, x$scale_vec, x$center, x$scale, subind)
    } else {
      .center_scale(x$Xp, x$center_vec, x$scale_vec, x$center, x$scale)
    }
  } else {
    if (!is.null(block_index)) {
      subind <- block_index_list(x$Xp)[[block_index]]
      .center_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, subind)  
    } else {
      .center_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale)  
    }
  }
}

#'@export
reverse_pre_process.block_matrix_pre_processor <- function(x, newdata) {
  .uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale)  
}


#' @export
pre_process.projector_pre_processor <- function(x, newdata=NULL, subind=NULL) {
  if (is.null(newdata)) {
    .center_scale(x$Xp, x$center_vec, x$scale_vec, x$center, x$scale, subind)
  } else {
    #.center_scale(x$projfun(newdata, subind=subind), x$center_vec, x$scale_vec, x$center, x$scale, subind) 
    .center_scale(x$projfun(newdata, subind=subind), x$center_vec, x$scale_vec, x$center, x$scale) 
  }
}

#' @export
reverse_pre_process.projector_pre_processor <- function(x, newdata, subind=NULL) {
  x$reconfun(.uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale), subind=subind) 
}

#' @export
pre_process.block_projector_pre_processor <- function(x, newdata=NULL, block_index=NULL) {
  assert_that(is.null(block_index) || length(block_index) == 1)
  if (is.null(newdata)) {
    if (!is.null(block_index)) {
      subind <- block_index_list(x)[[block_index]]
      .center_scale(x$Xp[,subind], x$center_vec, x$scale_vec, x$center, x$scale, subind)
    } else {
      .center_scale(x$Xp, x$center_vec, x$scale_vec, x$center, x$scale)
    }
  } else {
    #.center_scale(x$projfun(newdata, subind=subind), x$center_vec, x$scale_vec, x$center, x$scale, subind) 
    .center_scale(x$projfun(newdata, block_index=block_index), x$center_vec, x$scale_vec, x$center, x$scale) 
  }
}

#' @export
reverse_pre_process.block_projector_pre_processor <- function(x, newdata, block_index=NULL) {
  if (!is.null(block_index)) {
    subind <- block_index_list(x)[[block_index]]
    x$reconfun(.uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, subind))
  } else {
    .uncenter_scale(x$reconfun(newdata), x$center_vec, x$scale_vec, x$center, x$scale) 
  }
}

#' @importFrom matrixStats colSds
#' @export
pre_processor.matrix <- function(X, center=TRUE, scale=FALSE) {
  center_vec <- if (center) colMeans(X) else rep(0, ncol(X))
  scale_vec <- if (scale) matrixStats::colSds(X) else rep(1, ncol(X))
  
  structure(
    list(
      Xp=.center_scale(X, center_vec, scale_vec, center,scale),
      center=center,
      scale=scale,
      center_vec=center_vec,
      scale_vec=scale_vec),
      class="matrix_pre_processor")
  
}

#' @importFrom matrixStats colSds
#' @export
pre_processor.block_matrix <- function(X, center=TRUE, scale=FALSE) {
  center_vec <- if (center) colMeans(X) else rep(0, ncol(X))
  scale_vec <- if (scale) matrixStats::colSds(X) else rep(1, ncol(X))
  
  structure(
    list(
      Xp=.center_scale(X, center_vec, scale_vec, center,scale),
      center=center,
      scale=scale,
      center_vec=center_vec,
      scale_vec=scale_vec),
    class="block_matrix_pre_processor")
  
}
  

#' @export
pre_processor.projector <- function(X, center=TRUE, scale=FALSE) {
  projfun <- projection_fun(X)
  reconfun <- function(newdata) inverse_project(X, newdata)
  
  Xp <- project(X)
  
  center_vec <- if (center) colMeans(Xp) else rep(0, ncol(Xp))
  scale_vec <- if (scale) matrixStats::colSds(Xp) else rep(1, ncol(Xp))
  
  structure(
    list(Xp=.center_scale(Xp,center_vec, scale_vec, center,scale),
         center=center,
         scale=scale,
         center_vec=center_vec,
         scale_vec=scale_vec,
         reconfun=reconfun,
         projfun=projfun),
    class="projector_pre_processor")
  
}

#' @export
pre_processor.block_projection_matrix <- function(X, center=TRUE, scale=FALSE) {
  projfun <- projection_fun(X)
  reconfun <- function(newdata) inverse_project(X, newdata)
  
  Xp <- project(X)
  
  center_vec <- if (center) colMeans(Xp) else rep(0, ncol(Xp))
  scale_vec <- if (scale) matrixStats::colSds(Xp) else rep(1, ncol(Xp))
  
  structure(
    list(Xp=.center_scale(Xp,center_vec, scale_vec, center,scale),
         center=center,
         scale=scale,
         center_vec=center_vec,
         scale_vec=scale_vec,
         reconfun=reconfun,
         projfun=projfun),
    class="block_projector_pre_processor")
  
}



