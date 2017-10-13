
#' @export
super_block_matrix <- function(bxlist) {
  assertthat::assert_that(all(sapply(bxlist, is.block_matrix)))
  assertthat::assert_that(all(sapply(bxlist, nrow) == nrow(bxlist[[1]])))
  
  blockInd <- blockIndices(bxlist)
  P <- sum(sapply(bxlist, ncol))
  
  if (is.null(bxlist)) {
    names(bxlist) <- paste0("SB_", 1:length(bxlist))
  }
  
  attr(bxlist, "block_indices") <- blockInd
  attr(bxlist, "nblock") <- length(bxlist)
  attr(bxlist, "block_names") <- names(bxlist)
  attr(bxlist, "nsub_block") <- sapply(bxlist, nblocks)
  attr(bxlist, "nrow") <- nrow(bxlist[[1]])
  attr(bxlist, "ncol") <- P
  class(bxlist) <- c("super_block_matrix", "list") 
  
  bxlist
}

nrow.super_block_matrix <- function(x) {
  attr(x, "nrow")
  
}

ncol.super_block_matrix <- function(x) {
  attr(x, "ncol")
}

dim.super_block_matrix <- function(x) {
  c(attr(x, "nrow"), attr(x, "ncol"))
}

#' @export
nblocks.super_block_matrix <- function(x) {
  attr(x, "nblock")
}

block_lengths.super_block_matrix <- function(object) {
  bind <- attr(object, "block_indices")
  apply(bind, 1, diff)+1
}

block_index_list.super_block_matrix <- function(object) {
  bind <- attr(object, "block_indices")
  lapply(1:nrow(bind), function(i) seq(bind[i,1], bind[i,2]))
}

#' @export
get_block.super_block_matrix <- function(x, i, j=NULL) {
  ind <- attr(x, "block_indices")
  if (!is.null(j)) {
    get_block(x[[i]], j)
  } else {
    x[[i]]
  }
}

#' @export
as.list.super_block_matrix <- function(x, recursive=FALSE) {
  if (recursive) {
    lapply(1:attr(x, "nblock"), function(i) {
      as.list(x[[i]])
    })
  } else {
    lapply(1:attr(x, "nblock"), function(i) get_block(x, i))
  }
}

#' @export
as.matrix.super_block_matrix <- function(x) {
  do.call(cbind, lapply(x, as.matrix))
}

print.super_block_matrix <- function(object) {
  bind <- attr(object, "block_indices")
  
  cat("super_block_matrix", "\n")
  cat("  n super blocks: ", attr(object, "nblock"), "\n")
  cat("  n sub blocks: ", attr(object, "nsub_block"), "\n")
  cat("  nrows: ", dim(object)[1], "\n")
  cat("  ncols: ", dim(object)[2], "\n")
  cat("  super block cols: ", apply(bind, 1, diff)+1, "\n")
  cat("  super block names: ", attr(object, "block_names"))
}


#' @export
block_apply.super_block_matrix <- function(x, f, descend=FALSE) {
  ret <- if (descend) {
    lapply(1:nblocks(x), function(i) {
      b <- get_block(x,i)
      block_apply(b, f)
    })
  } else {
    lapply(1:nblocks(x), function(i) {
      f(get_block(x,i), i)
    })
  }
  
  super_block_matrix(ret)
  
}

#' @export
names.block_matrix <- function(x) attr(x, "block_names")



