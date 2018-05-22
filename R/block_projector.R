
#' #' projector_list
#' projector_list <- function(Xs) {
#'   assert_that(all(sapply(Xs, function(x) inherits(x, "projector"))))
#'   structure(Xs,
#'             block_indices=block_indices(Xs),
#'             class=c("projector_list", "list"))
#' }
#' 
#' #' @export
#' project.projector_list(x, newdata, block_index=NULL) {
#'   if (!is.null(block_index)) {
#'     assert_that(length(block_index) == 1 || length(block_index) == length(x))
#'   }
#'   if (is.null(block_index)) {
#'     
#'   }
#' }


#' block_projector 
#' 
#' a blockwise concatenation of a list of \code{projectors} each with the same row-dimension
#' 
#' @param Xs a list of \code{projector} objects.
#' 
block_projector <- function(Xs) {
  assertthat::assert_that(all(sapply(Xs, function(x) inherits(x, "projector"))))
  assertthat::assert_that(all(sapply(Xs, nrow) == nrow(Xs[[1]])))
  
  ## get the fitted projections
  Xproj <- lapply(Xs, project)
  
  ## the projected blocks
  proj_ind <- block_indices(Xproj)
  
  ## the original blocks
  block_ind <- block_indices(Xs)
  
  P <- sum(sapply(Xs, ncol))
  projP <- sum(sapply(Xproj, ncol))
  
  attr(Xs, "block_data") <- block_matrix(Xproj)
  attr(Xs, "block_indices") <- block_ind
  attr(Xs, "proj_indices") <- proj_ind
  attr(Xs, "nblock") <- length(Xs)
  attr(Xs, "nrow") <- nrow(Xs[[1]])
  attr(Xs, "ncol") <- P
  attr(Xs, "proj_ncol") <- projP
  attr(Xs, "block_names") <- names(Xs)
  class(Xs) <- c("block_projector", "block_matrix", "projector", "list") 
  Xs
}

#' @export
scores.block_projector <- function(x) {
  att(x, "block_data")
}


#' @export
get_block.block_projector <- function(x, i) {
  get_block(attr(x, "block_data"), i)
}


#' @export
dim.block_projector <- function(x) {
  c(attr(x, "nrow"), attr(x, "ncol"))
}

#' @export
nrow.block_projector <- function(x) {
  attr(x, "nrow")
}

#' @export
ncol.block_projector <- function(x) {
  attr(x, "ncol")
}

#' @export
ncomp.block_projector <- function(x) {
  attr(x, "proj_ncol") 
}

#' @export
as.matrix.block_projector <- function(x) {
  attr(x, "block_data")
}


#' @export
projection_fun.block_projector <- function(x, pre_process=TRUE) {
  function(newdata, block_index=NULL) {
    project(x, newdata, block_index=block_index)
  }
}

#' @export
project.block_projector <- function(x, newdata=NULL, block_index=NULL, ...) {
  if (!is.null(block_index)) {
    assert_that(length(block_index) == 1, msg="block_index must be of length 1")
  }
  
  if (is.null(newdata)) {
    if (is.null(block_index)) {
      return(attr(x, "block_data"))
    } else {
      get_block(x,block_index)
    }
  } else {
    assert_that(is.matrix(newdata) || inherits(newdata, "Matrix"))
    oind <- attr(x, "block_indices")
    if (is.null(block_index)) {
      block_index <- 1:nblocks(x)
      len <- sum(sapply(block_index, function(i) diff(oind[i,])+1))
      assert_that(ncol(newdata) == len)
    
      oind <- attr(x, "block_indices")
      block_index <- 1:nblocks(x)
      block_matrix(lapply(block_index, function(i) {
        ind <- seq(oind[i,1], oind[i,2])
        project(x[[i]], newdata=newdata[,ind],...)
      }))
    } else {
      assert_that(block_lengths(x)[block_index] == ncol(newdata))
      block_matrix(list(project(x[[block_index]], newdata=newdata,...)))
    }
  }
  
}


#' @export
print.block_projector <- function(object) {
  bind <- attr(object, "block_indices")
  
  cat("block_projector", "\n")
  cat("  nblocks: ", attr(object, "nblock"), "\n")
  cat("  nrows: ", nrow(object), "\n")
  cat("  ncols (input): ", ncol(object), "\n")
  cat("  ncols (output): ", attr(object, "proj_ncol"), "\n")
  cat("  block cols (input): ", apply(bind, 1, diff)+1, "\n")
  cat("  block cols (output): ", sapply(attr(object, "block_data"), ncol), "\n")
  cat("  block names: ", attr(object, "block_names"))
}

#' @export
compose.block_projector <- function(x, y) {
  projector_list <- lapply(1:nblocks(x), function(i) {
    function(newdata=NULL) {
      oind <- attr(x, "proj_indices")[i,]
      oind <- seq(oind[1], oind[2])
      project(y, project(x, newdata, block_index=i), subind=oind)
    }
  })

  
  global_projector <- function(newdata) {
    project(y, project(x, newdata))
  }
  
  recon_list <- lapply(1:nblocks(x), function(i) {
    function(newdata=NULL) {
      reconstruct(x[[i]], reconstruct(y, newdata))
    }
  })
  
  global_recon <- function(newdata=NULL) {
    reconstruct(x, reconstruct(y, newdata))
  }
  
  ret <- structure(
    list(
      a=x,
      b=y,
      d=y$d,
      u=y$u,
      ncomp=ncomp(y),   
      scores=scores(y),
      global_projector=global_projector,
      projector_list=projector_list),
    class=c("composed_block_projector", "block_projector", "projector"))
    
    attr(ret, "block_indices") <- block_indices(x)
    ret
}
  

