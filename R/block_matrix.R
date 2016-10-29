
blockIndices <- function(Xlist) {
  ncols <- sapply(Xlist, ncol)
  csum <- cumsum(ncols)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  m
}



#' block_matrix_list
#' @param Xs a list of k matrices with N_k rows and M_k columns 
#' @export
#' @importFrom assert assertthat
block_matrix_list <- function(Xs) {
  assertthat::assert_that(all(sapply(Xs, is.matrix)))
  assertthat::assert_that(all(sapply(Xs, nrow) == nrow(Xs[[1]])))
  
  blockInd <- blockIndices(Xs)
  P <- sum(sapply(Xs, ncol))
  
  attr(Xs, "block_indices") <- blockInd
  attr(Xs, "nblock") <- length(Xs)
  class(Xs) <- c("block_matrix_list", "list") 
  
  Xs
}

#' block_matrix
#' @param Xs a list of k matrices with N rows and M_k columns 
#' @export
#' @importFrom assert assertthat
block_matrix <- function(Xs) {
  assertthat::assert_that(all(sapply(Xs, is.matrix)))
  assertthat::assert_that(all(sapply(Xs, nrow) == nrow(Xs[[1]])))
  
  blockInd <- blockIndices(Xs)
  P <- sum(sapply(Xs, ncol))
  
  X <- do.call(cbind, Xs)
  attr(X, "block_indices") <- blockInd
  attr(X, "nblock") <- length(Xs)
  class(X) <- c("block_matrix", "matrix") 
  
  X
}

#' @export
#' @param x the object to retrive the block of variables from.
#' @param i the index of the block.
get_block <- function (x, i) { UseMethod("get_block") }


get_block.block_matrix <- function(x, i) {
  ind <- attr(x, "block_indices")
  X[, seq(ind[i,1], ind[i,2]) ]
}

get_block.block_matrix_list <- function(x, i) {
  x[[i]]
}

as.list.block_matrix <- function(x) {
  lapply(1:attr(x, "nblock"), function(i) get_block(x, i))
}

as.list.block_matrix_list <- function(x) {
  x
}

split_matrix <- function(X, fac) {
  idx <- split(1:nrow(X), fac)
  lapply(idx, function(i) X[i,])
}

