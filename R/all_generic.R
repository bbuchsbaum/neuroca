
#' scores
#' 
#' extract the factor score matrix for a multivariate model
#' @param ... extra args
#' @param x the model fit
#' @export
scores <- function(x,...) UseMethod("scores")


#reduce <- function(x, Y, ...) UseMethod("reduce")

#' loadings
#' 
#' extract the loadings matrix (the variable coefficients) for a multivariate model.
#' @param x the model fit
#' @param ... extra args
#' @export
loadings <- function(x,...) UseMethod("loadings")

#' ncomp get the number of components in the estimated model
#' 
#' @param x the model fit
#' @export
ncomp <- function(x) UseMethod("ncomp")


#ncol          <- function(x) UseMethod("ncol") 
#ncol.default  <- base::ncol


#nrow          <- function(x) UseMethod("nrow") 
#nrow.default  <- base::nrow

#dim          <- function(x) UseMethod("dim") 
#dim.default  <- base::dim



#' @export
correlations <- function(x,...) UseMethod("correlations")


#' cross_validate a model
#' 
#' @param the model fit
#' @param ... extra args
#' @export
cross_validate <- function(x, ...) UseMethod("cross_validate")


#' @export
nested_cv <- function(x, ...) UseMethod("nested_cv")

#' bootstrap a model
#' @param x the model fit
#' @param nboot the number of bootstrap resamples
#' @param ... extra args
#' @export
bootstrap <- function(x, nboot, ...) UseMethod("bootstrap")

#' resample data from a model fit
#' @param x the model fit
#' @param ... extra args
#' @export
resample <- function(x, ...) UseMethod("resample")

#' block_apply
#' apply a function to each block of a multi-block data structure.
#' @param x the multi-block data
#' @param f the function to apply
#' @param ... extra args
#' @export
block_apply <- function(x, f, ...) UseMethod("block_apply")

#' nblocks
#' extract the number of blocks in a mutli-block data structure or model
#' 
#' @param x the object
#' @export
nblocks <- function(x) UseMethod("nblocks")


#' @export
#' @param x the model fit object
#' @param ... extra args
permutation <- function(x, ...) UseMethod("permutation")

#' @export
#' @param x the model fit object
#' @param type the tpye of variance to account for.
#' @param ncomp the number of components to include.
vaf <- function(x, type, ncomp) UseMethod("vaf")


#' @param x the model fit object
#' @export 
performance <- function(x, yobs, ncomp, folds, metric, ...) UseMethod("performance")

#' @param x the model fit object
#' @export
permute_refit <- function(x, ...) UseMethod("permute_refit")

#' @param x the model fit object
#' @export 
permute <- function(x, ...) UseMethod("permute")

#' @param x the model fit object
#' @export
refit <- function(x, ...) UseMethod("refit")

#' @param x the model fit object
#' @export
reprocess <- function(x, ...) UseMethod("reprocess")

#' @param x the model fit object
#' @export
supplementary_predictor <- function(x, ...) UseMethod("supplementary_predictor")

#' @param x the model fit object
#' @export
supplementary_loadings <- function(x,...) UseMethod("supplementary_loadings")


#' @param x the model fit object
#' @export
summarize_loadings <- function(x, stat, ncomp, ...) UseMethod("summarize_loadings")

#' @param x the model fit object
#' @export
compose <- function(x,y) UseMethod("compose")

#' @param x the model fit object
#' @export
concat <- function(x,...) UseMethod("concat")

#' @param x the model fit object
#' @export
combine <- function(x,...) UseMethod("combine")

#' @export
singular_values <- function(x) UseMethod("singular_values")

#' @param x the model fit object
#' @export
truncate <- function(x, ncomp) UseMethod("truncate")

#' @param x the model fit object
#' @export 
subset_rows <- function(x, idx) UseMethod("subset_rows")


#' partial_scores
#' 
#' compute the partial scores from a multivariate model using a subset of the input
#' 
#' @param x the model fit
#' @export
partial_scores <- function(x, ...) UseMethod("partial_scores")

#' contributions
#' 
#' compute the contributions (of observations, variables, tables) to a model. 
#' @param x the model fit
#' @param ... extra args
#' @export
contributions <- function(x, ...) UseMethod("contributions")


#' reproducibility
#' 
#' compute a measure of the reproducibility of a model under replication.
#' @param x the model fit
#' @param folds
#' @param metric
#' @param ...
#' @export
reproducibility <- function(x, folds, metric, ...) UseMethod("reproducibility")

#' reconstruct the data with some number of components
#' 
#' @param x the model fit
#' @param newdata newdata to be inverse projected (optional)
#' @param comp the components to use
#' @param ... extra args
#' @export
reconstruct <- function(x, newdata, comp,...) UseMethod("reconstruct")

#' get the residuals of a model, after removing the first \code{ncomp} components
#' 
#' @param x the model fit
#' @param ncomp the number of components
#' @param ... extra arguments
residuals <- function(x, ncomp, ...) UseMethod("residuals")


#' project_cols
#' 
#' project supplementary variables in to the subspace defined by the model
#' 
#' @param x the model fit
#' @param newdata a matrix or vector of new variables(s)
#' @param ... extra args
#' @export
project_cols <- function(x, newdata, ...) UseMethod("project_cols")


#' project
#' 
#' project supplementary observations in to the subspace defined by the model
#' 
#' @param x the model fit
#' @param newdata a matrix or vector of new obervations(s)
#' @param ... extra args
#' @export
project <- function(x, newdata, ...) UseMethod("project")


#' project a selected subset of indices onto the subspace
#' 
#' @inheritParams project
#' @param the column indices to select in the projection matrix
#' @export
partial_project <- function(x, newdata, colind) UseMetod("partial_project")


#' project a single 'block' of data onto the subspace
#' 
#' @inheritParams project
#' @param block the block id to select in the block projection matrix
#' @export
block_project <- function(x, newdata, block,...) UseMethod("block_project")


#' partial_projector
#' 
#' extract a partial projector from a fitted model
#' 
#' @param x the model fit
#' @param colind a vector of unique subindices
#' @param ... extra args
#' @export
partial_projector <- function(x, colind, ...) UseMethod("partial_projector")

#' projection_fun
#' 
#' return a function that projects data to lower-dimensional space
#' 
#' @export
projection_fun <- function(x, colind, ...) UseMethod("projection_fun")

#' project_table
#' 
#' project a block of data into the subspace defined by the model.
#' 
#' @export
project_table <- function(x, supY, supX, ncomp, ...) UseMethod("project_table")

#' @export
procrusteanize <- function(x,...) UseMethod("procrusteanize")


#' @export
pre_process <- function(obj, X, ...) UseMethod("pre_process")

#' @export
reverse_pre_process <- function(obj, X, ...) UseMethod("reverse_pre_process")

#' block_lengths
#' 
#' extract the vector of lengths of each block in a multi-block object
#' @param object the object
#' @export
block_lengths <- function(object) UseMethod("block_lengths")

#' block_index_list
#' 
#' extract the list of indices associated with each block in a multi-block object
#' @param object the object
#' @export
block_index_list <- function(object) UseMethod("block_index_list")


#' @export
project_copy <- function(x, ...) UseMethod("project_copy")


#' @param center whether to mean center the columns
#' @param scale whether to scale the columns
#' @export
pre_processor <- function(x, center, scale) UseMethod("pre_processor")

#' turn a model object into a classifier
#' 
#' @param x the fitted object
#' @export
classifier <- function(x, ...) UseMethod("classifier")


#' apply a rotation matrix to a solution
#' 
#' @param x the object to rotate
#' @param rot the rotation matrix to apply
#' @param ... extra arguments
#' @export
rotate <- function(x, rot,...) UseMethod("rotate")


#' get_block
#' 
#' @export
#' @param x the object to retrive the block of variables from.
#' @param i the index of the block.
get_block <- function (x, i,...) { UseMethod("get_block") }

