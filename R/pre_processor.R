
.uncenter_scale <- function(newdata, center_vec, scale_vec, center, scale, colind=NULL) {
  if (is.null(colind)) {
    colind <- 1:ncol(newdata)
  }
  
  if (!center && !scale) {
    return(newdata[,colind,drop=FALSE])
  }
  
  if (is.null(colind)) {
    colind <- 1:length(center_vec)
  }
  
  if (center && scale) {
    m1 <- sweep(newdata, 2, scale_vec[colind], "*")
    sweep(m1, 2, center_vec[colind], "+")
  } else if (!center && scale) {
    sweep(newdata, 2, scale_vec[colind], "*")
  } else if (center && !scale) {
    sweep(newdata, 2, center_vec[colind], "+")
  } else {
    stop()
  }
}

.center_scale <- function(newdata, center_vec, scale_vec, center, scale, colind=NULL) {
  if (is.null(colind)) {
    colind <- 1:ncol(newdata)
  }
  
  if (!center && !scale) {
    return(newdata[,colind,drop=FALSE])
  }
  
  assert_that(length(colind) == ncol(newdata), 
              msg=paste("center_scale: length(colind) = ", length(colind), 
                        "and ncol(newdata) = ", ncol(newdata)))
  
  if (center && scale) {
    m1 <- sweep(newdata, 2, center_vec[colind], "-")
    sweep(m1, 2, scale_vec[colind], "/")
  } else if (!center && scale) {
    sweep(newdata, 2, scale_vec[colind], "/")
  } else if (center && !scale) {
    sweep(newdata, 2, center_vec[colind], "-")
  } else {
    stop()
  }
}

#' @export
pre_process.no_op_matrix_pre_processor <- function(x, newdata, colind=NULL) {
  if (is.null(colind)) {
    newdata
  } else {
    newdata[,colind]
  }
}

#' @export
reverse_pre_process.no_op_matrix_pre_processor <- function(x, newdata, colind=NULL) {
  if (is.null(colind)) {
    newdata
  } else {
    newdata[,colind,drop=FALSE]
  }
}

#' @export
pre_process.matrix_pre_processor <- function(x, newdata, colind=1:length(x$center_vec)) {
  .center_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, colind)  
}

#' @export
reverse_pre_process.matrix_pre_processor <- function(x, newdata, colind=NULL) {
  .uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, colind)  
}

#' @export
pre_process.block_matrix_pre_processor <- function(x, newdata=NULL,colind=NULL) {
  if (is.null(newdata) && is.null(colind)) {
    x$Xp
  } else {
    ## return a block_matrix
    .center_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, colind)  
  }
}

#'@export
reverse_pre_process.block_matrix_pre_processor <- function(x, newdata,colind=NULL) {
  .uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale)  
}


#' @export
pre_process.projector_pre_processor <- function(x, newdata=NULL, colind=NULL) {
  if (is.null(newdata)) {
    .center_scale(x$Xp, x$center_vec, x$scale_vec, x$center, x$scale, colind)
  } else {
    #.center_scale(x$projfun(newdata, colind=colind), x$center_vec, x$scale_vec, x$center, x$scale, colind) 
    .center_scale(x$projfun(newdata, colind=colind), x$center_vec, x$scale_vec, x$center, x$scale) 
  }
}

#' @export
reverse_pre_process.projector_pre_processor <- function(x, newdata, colind=NULL) {
  xrecon <- x$reconfun(newdata)
  .uncenter_scale(xrecon, x$center_vec, x$scale_vec, x$center, x$scale) 
}

#' @export
pre_process.block_projector_pre_processor <- function(x, newdata=NULL, block_index=NULL) {
  assert_that(is.null(block_index) || length(block_index) == 1)
  if (is.null(newdata)) {
    if (!is.null(block_index)) {
      colind <- block_index_list(x)[[block_index]]
      .center_scale(x$Xp[,colind], x$center_vec, x$scale_vec, x$center, x$scale, colind)
    } else {
      .center_scale(x$Xp, x$center_vec, x$scale_vec, x$center, x$scale)
    }
  } else {
    #.center_scale(x$projfun(newdata, colind=colind), x$center_vec, x$scale_vec, x$center, x$scale, colind) 
    .center_scale(x$projfun(newdata, block_index=block_index), x$center_vec, x$scale_vec, x$center, x$scale) 
  }
}

#' @export
reverse_pre_process.block_projector_pre_processor <- function(x, newdata, block_index=NULL) {
  if (!is.null(block_index)) {
    colind <- block_index_list(x)[[block_index]]
    x$reconfun(.uncenter_scale(newdata, x$center_vec, x$scale_vec, x$center, x$scale, colind))
  } else {
    .uncenter_scale(x$reconfun(newdata), x$center_vec, x$scale_vec, x$center, x$scale) 
  }
}


#' @importFrom matrixStats colSds
#' @export
pre_processor.matrix <- function(X, center=TRUE, scale=FALSE) {
  
  if (center == FALSE && scale == FALSE) {
    structure(
      list(
        center=center,
        scale=scale,
        center_vec=0,
        scale_vec=1),
      class="no_op_matrix_pre_processor")
  } else {
    structure(
      list(
        center=center,
        scale=scale,
        center_vec=colMeans(X),
        scale_vec=matrixStats::colSds(X)),
        class="matrix_pre_processor")
  }
  
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
  reconfun <- function(newdata, .colind=NULL) reconstruct(X, newdata, colind=.colind)
  
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
pre_processor.block_projector <- function(X, center=TRUE, scale=FALSE) {
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



