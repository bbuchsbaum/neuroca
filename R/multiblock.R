
#' @export
nblocks.multiblock <- function(x) length(block_index_list(x))


#' @export
project.multiblock <- function(x, newdata, comp=1:ncomp(x), pre_process=TRUE, table_index=NULL) {
  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol=length(newdata))
  }
  
  assert_that(is.matrix(newdata))
  
  if (is.null(table_index)) {
    # new data must have same number of columns as original data
    assert_that(ncol(newdata) == x$nvars)
    xnewdat <- x$reprocess(newdata)
    project(x$fit, xnewdat, comp=comp)
  } else if (length(table_index) == 1) {
    ind <- x$block_indices[[table_index]]
    assert_that(length(ind) == ncol(newdata))
    xnewdat <- x$reprocess(newdata, table_index)
    x$ntables * project(x$fit, xnewdat, comp=comp, subind=x$block_indices[[table_index]])
    
  } else {
    stop("table_index must have length of 1")
  }
  
}


## project from existing table
#' @export
predict.multiblock <- function(x, newdata, ncomp=x$ncomp, table_index=1:x$ntables, pre_process=TRUE) {
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  fscores <- if (length(table_index) == x$ntables) {
    assert_that(ncol(newdata) == ncol(x$X))
    Reduce("+", lapply(table_index, function(i) {
      ind <- x$block_indices[[i]]
      
      Xp <- if (pre_process) {
        x$reprocess(newdata[, ind], i)
      } else {
        newdata[, ind]
      }
      
      project(x, Xp, comp=1:ncomp) * x$ntables
    }))
  } else if (length(table_index) == 1) {
    ind <- x$block_indices[[table_index]]
    
    Xp <- if (pre_process) {
      x$reprocess(newdata, table_index)
    } else {
      newdata
    }
    
    fscores <- project(x, Xp, comp=1:ncomp, subind=ind) * x$ntables
    
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
ncomp.multiblock <- function(x) ncomp(x$fit)

#' @export
scores.multiblock <- function(x) scores(x$fit)

#' @export
loadings.multiblock <- function(x) loadings(x$fit)

#' @export 
block_index_list.multiblock <- function(x) x$block_indices

#' @export
partial_scores.multiblock <- function(x, table_index=1:x$ntables) {
  bind <- block_index_list(x)
  res <- lapply(table_index, function(i) {
    x$ntables * project(x$fit, get_block(Xr, i), subind=bind[[i]])
  })
  
  names(res) <- paste0("Block_", table_index)
  res
}

reprocess.multiblock <- function(x, newdat, table_index=NULL) {
  ## given a new observation(s), pre-process it in the same way the original observations were processed
  
  prep <- attr(x$X, "pre_process")
  
  newdat <- if (is.null(table_index)) {
    prep(newdat)
  } else {
    prep(newdat, bind[[table_index]])
  }
  
  if (x$is_reduced && !is.null(table_index)) {
    project(x$reducer, newdat, table_index)
  } else if (is_reduced & is.null(table_index)) {
    project(x$reducer, newdat)
  } else {
    newdat
  }
  
}



