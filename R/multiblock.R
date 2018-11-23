
#' @export
nblocks.multiblock <- function(x) length(block_index_list(x))


#' @export
project.multiblock <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, block_index=NULL) {
  if (missing(newdata)) {
    if (is.null(block_index)) {
      return(scores(x)[,comp])
    } else {
      return(partial_scores(x,block_index=block_index)[,comp])
    }
  }
  
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  
  assert_that(is.matrix(newdata))
  
  if (is.null(block_index)) {
    # new data must have same number of columns as original data
    assert_that(ncol(newdata) == x$nvars)
    xnewdat <- if (pre_process) reprocess(x, newdata) else newdata
    project(x$fit, xnewdat, comp=comp)
  } else if (length(block_index) == 1) {
    ind <- x$block_indices[[block_index]]
    assert_that(length(ind) == ncol(newdata))
    xnewdat <- if (pre_process) reprocess(x, newdata, block_index) else newdata
    x$ntables * project(x$fit, xnewdat, comp=comp, colind=ind)
    
  } else {
    stop("block_index must have length of 1 or length equal to ntables")
  }
  
}


#' @export
predict.multiblock <- function(x, newdata, ncomp=x$ncomp, block_index=1:x$ntables, pre_process=TRUE) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  
  assert_that(is.matrix(newdata))
  assert_that(length(block_index) == 1 || length(block_index) == x$ntables)
  
  fscores <- if (length(block_index) == x$ntables) {
    assert_that(ncol(newdata) == ncol(x$X))
    Reduce("+", lapply(block_index, function(i) {
      project(x, newdata[,ind,drop=FALSE], comp=1:ncomp, pre_process=pre_process, block_index=i) 
    }))
  } else if (length(block_index) == 1) {
    fscores <- project(x, newdata, comp=1:ncomp, block_index=block_index)
    
  }
  
}


#' @export
contributions.multiblock <- function(x, type=c("column", "row", "table")) {
  contr <- contributions(x$fit)
  type <- match.arg(type)
  
  if (type == "table") {
    out <- do.call(cbind, lapply(1:x$ncomp, function(i) {
      sapply(1:x$ntables, function(j) {
        ind <- x$block_indices[[j]]
        sum(contr[i,ind])
      })
    }))
    
  } else if (type == "row") {
    contributions(x$fit, type="row")
  } else {
    contributions(x$fit, type="column")
  }
}

#' @export
ncol.multiblock <- function(x) x$xdim[2]

#' @export
nrow.multiblock <- function(x) x$xdim[1]

#' @export
ncomp.multiblock <- function(x) ncomp(x$fit)

#' @export
scores.multiblock <- function(x) scores(x$fit)

#' @export
loadings.multiblock <- function(x) loadings(x$fit)

#' @export
block_indices <- function(x,i) {
  x$block_indices[[i]]
}

#' @export 
block_index_list.multiblock <- function(x) x$block_indices

#' @export
partial_scores.multiblock <- function(x, block_index=1:x$ntables) {
  bind <- block_index_list(x)
  res <- lapply(block_index, function(i) {
    x$ntables * project(x$fit, get_block(Xr, i), colind=bind[[i]])
  })
  
  names(res) <- paste0("Block_", block_index)
  res
}


#' @export
reprocess.multiblock <- function(x, newdata, colind=NULL) {
  if (is.null(colind)) {
    assert_that(ncol(newdata) == nrow(loadings(x)))
    x$preproc$transform(newdata)
  } else {
    assert_that(length(colind) == ncol(newdata), 
                msg=paste("length of colind not equal to number of columns of newdata", 
                          length(colind), "!=", ncol(newdata)))
    x$preproc$transform(newdata, colind)
  }
  
}



