
#' @export
nblocks.multiblock <- function(x) length(block_index_list(x))


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

# predict.multiblock <- function(x, newdata, ncomp=ncomp(x), table_index=1:nblocks(x), pre_process=TRUE) {
#   assert_that(is.matrix(newdata))
#   assert_that(length(table_index) == 1 || length(table_index) == nblocks(x))
#   
#   indlist <- block_index_list(x)
#   
#   fscores <- if (length(table_index) == nblocks(x)) {
#     ncols <- sum(sapply(indlist, length))
#     
#     assert_that(ncol(newdata) == ncols)
#     Reduce("+", lapply(table_index, function(i) {
#       ind <- x$block_indices[[i]]
#       
#       Xp <- if (pre_process) {
#         x$reprocess(newdata[, ind], i)
#       } else {
#         newdata[, ind]
#       }
#       
#       project(x, Xp, ncomp=ncomp, table_index=i) 
#     }))
#   } else if (length(table_index) == 1) {
#     ind <- x$block_indices[[table_index]]
#     
#     Xp <- if (pre_process) {
#       x$reprocess(newdata, table_index)
#     } else {
#       newdata[, ind]
#     }
#     
#     fscores <- project(x, Xp, ncomp=ncomp, table_index=table_index) 
#     
#   }
#   
# }