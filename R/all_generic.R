


#' @export
scores <- function(x,...) UseMethod("scores")

#' @export
reduce <- function(x, Y, ...) UseMethod("reduce")

#' @export
loadings <- function(x,...) UseMethod("loadings")


#' @export
correlations <- function(x,...) UseMethod("correlations")

#' @export
cross_validate <- function(x, ...) UseMethod("cross_validate")

#' @export
nested_cv <- function(x, ...) UseMethod("nested_cv")

#' @export
bootstrap <- function(x, niter, ...) UseMethod("bootstrap")

#' @export
resample <- function(x, ...) UseMethod("resample")

#' @export
block_apply <- function(x, f, ...) UseMethod("block_apply")

#' @export
nblocks <- function(x) UseMethod("nblocks")



#' @export
jackstraw <- function(x, nsynth, niter, ...) UseMethod("jackstraw")

#' @export
permutation <- function(x, ...) UseMethod("permutation")

#' @export
permute_refit <- function(x, ...) UseMethod("permute_refit")

#' @export
project_cols <- function(x, ncomp,...) UseMethod("project_cols")

#' @export
split_half_reliability <- function(x, ...) UseMethod("split_half_reliability")

#' @export
supplementary_predictor <- function(x, ...) UseMethod("supplementary_predictor")

#' @export
optimal_components <- function(x, ...) UseMethod("optimal_components")

#' @export
singular_values <- function(x) UseMethod("singular_values")

#' @export
eigen_values <- function(x) UseMethod("eigen_values")

#' @export
partial_scores <- function(x, ...) UseMethod("partial_scores")

contributions <- function(x, ...) UseMethod("contributions")

#' @export
reproducibility <- function(x, folds, metric, ...) UseMethod("reproducibility")

#' @export
reconstruct <- function(x, ncomp) UseMethod("reconstruct")

#' @export
project <- function(x, newX, ...) UseMethod("project")

#' @export
reduce_rank <- function(x, k, ...) UseMethod("reduce_rank")
# pre_process <- function(obj, X, ...) UseMethod("pre_process")



