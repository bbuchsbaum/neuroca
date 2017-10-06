
blockIndices <- function(Xlist) {
  ncols <- sapply(Xlist, ncol)
  csum <- cumsum(ncols)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  m
}



#' block_matrix_list
#' 
#' @param Xs a list of k matrices each with with N_k rows and M_k columns 
#' @export
#' @importFrom assertthat assert_that
block_matrix_list <- function(Xs) {
  assertthat::assert_that(all(sapply(Xs, is.matrix)))
  assertthat::assert_that(all(sapply(Xs, nrow) == nrow(Xs[[1]])))
  
  blockInd <- blockIndices(Xs)
  P <- sum(sapply(Xs, ncol))
  
  attr(Xs, "block_indices") <- blockInd
  attr(Xs, "nblock") <- length(Xs)
  attr(Xs, "nrow") <- nrow(Xs[[1]])
  attr(Xs, "ncol") <- P
  attr(Xs, "block_names") <- names(Xs)
  class(Xs) <- c("block_matrix_list", "block_matrix", "list") 
  
  Xs
}

#' @export
to_block_matrix <- function(X, block_lengths) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
  assertthat::assert_that(sum(block_lengths) == ncol(X))
  
  csum <- cumsum(block_lengths)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  blockInd <- m
  
  attr(X, "block_indices") <- blockInd
  attr(X, "nblock") <- length(block_lengths)
  attr(X, "block_names") <- paste0("B", 1:length(block_lengths))
  class(X) <- c("block_matrix", "matrix") 
  
  X
  
}


#' block_matrix
#' 
#' @param Xs a list of k matrices with N rows and M_k columns 
#' @export
#' @importFrom assertthat assert_that
block_matrix <- function(Xs) {
  assertthat::assert_that(all(sapply(Xs, is.matrix)))
  assertthat::assert_that(all(sapply(Xs, nrow) == nrow(Xs[[1]])))
  
  blockInd <- blockIndices(Xs)
  P <- sum(sapply(Xs, ncol))
  
  X <- do.call(cbind, Xs)
  attr(X, "block_indices") <- blockInd
  attr(X, "nblock") <- length(Xs)
  attr(X, "block_names") <- names(Xs)
  class(X) <- c("block_matrix", "matrix") 
  
  X
}


#' @export
matrix_to_block_matrix <- function(X, groups) {
  assert_that(length(groups) == ncol(X))
  glevs <- sort(unique(groups))
  Xlist <- lapply(glevs, function(i) {
    idx <- which(groups==i)
    x <- X[, idx]
  })
  
  block_matrix(Xlist)
}

nrow.block_matrix_list <- function(x) {
  attr(x, "nrow")
  
}

ncol.block_matrix_list <- function(x) {
  attr(x, "ncol")
}

dim.block_matrix_list <- function(x) {
  c(attr(x, "nrow"), attr(x, "ncol"))
}



block_lengths.block_matrix <- function(object) {
  bind <- attr(object, "block_indices")
  apply(bind, 1, diff)+1
}

block_index_list.block_matrix <- function(object) {
  bind <- attr(object, "block_indices")
  lapply(1:nrow(bind), function(i) seq(bind[i,1], bind[i,2]))
}


print.block_matrix <- function(object) {
  bind <- attr(object, "block_indices")
  
  cat("block_matrix", "\n")
  cat("  nblocks: ", attr(object, "nblock"), "\n")
  cat("  nrows: ", dim(object)[1], "\n")
  cat("  ncols: ", dim(object)[2], "\n")
  cat("  block cols: ", apply(bind, 1, diff)+1, "\n")
  cat("  block names: ", attr(object, "block_names"))
}

#' get_block
#' 
#' @export
#' @param x the object to retrive the block of variables from.
#' @param i the index of the block.
get_block <- function (x, i,...) { UseMethod("get_block") }


#' @export
get_block.block_matrix <- function(x, i) {
  ind <- attr(x, "block_indices")
  x[, seq(ind[i,1], ind[i,2]) ]
}

#' @export
get_block.block_matrix_list <- function(x, i) {
  x[[i]]
}

#' @export
as.list.block_matrix <- function(x) {
  lapply(1:attr(x, "nblock"), function(i) get_block(x, i))
}

#' @export
as.list.block_matrix_list <- function(x) {
  x
}

#' @export
as.matrix.block_matrix_list <- function(x) {
  block_matrix(x)
}


#' @export
nblocks.block_matrix <- function(x) {
  attr(x, "nblock")
}


#' @export
rbind.block_matrix <- function(...) {
  mlist <- list(...)
  nb <- unlist(lapply(mlist, nblocks))
  
  assert_that(all(nb[1] == nb))
  res <- lapply(1:nb[1], function(bnum) {
    do.call(rbind, lapply(mlist, function(x) get_block(x, bnum)))
  })

  block_matrix(res)
}

#' @export
block_apply.block_matrix <- function(x, f) {
  ret <- lapply(1:nblocks(x), function(i) {
    f(get_block(x,i), i)
  })
  
  block_matrix_list(ret)
}

#' @export
names.block_matrix <- function(x) attr(x, "block_names")

#' @export
is.block_matrix <- function(x) { inherits(x, "block_matrix") }


#' @export
reduce_rank.block_matrix <- function(x, k, center=TRUE, scale=FALSE) {
  
  nb <- nblocks(x)
  bind <- attr(x, "block_indices")
  
  if (length(k) == 1) {
    k <- rep(k, nblocks(x))
  }
  
  assert_that(length(k) == nblocks(x))
  
  pcres <- lapply(1:nblocks(x), function(i) {
    print(i)
    xb <- get_block(x, i)
    k <- min(min(dim(xb)), k[i])
    pca(get_block(x, i), k, center=center, scale=scale, svd.method="propack")
  })
  
  
  components <- lapply(pcres, function(x) scores(x))
  bm <- block_matrix_list(components)
  
  projector <- function(x,i) {
    project(pcres[[i]], x)
  }
  
  global_projector <- function(x) {
    out <- lapply(1:nb, function(i) {
      xi <- x[, bind[i,1]:bind[i,2]]
      projector(xi,i)
    })
    
    block_matrix_list(out)
  }
  
  ret <- list(x=bm, block_projector=projector, global_projector=global_projector)
  class(ret) <- c("reduced_rank_block_matrix", "list")
  ret
}

nblocks.reduced_rank_block_matrix <- function(x) {
  nblocks(x$x)
}


#' @export
project.reduced_rank_block_matrix <- function(x, newX, i) {
  if (missing(i)) {
    x$global_projector(newX)
  } else {
    assert_that(length(i) == 1 && i >0 && i <= nblocks(x))
    x$block_projector(newX,i)
  }
}





