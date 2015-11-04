#' PLS_Config
#' Configuration information for PLS and related mutivariate models
#' 
#' @param data_files the names of the image data files 
#' @param design the file or list of files containing the design matrix information
#' @param folding_dimension name of dimension will be folded into the columns of the X matrix.
#' @param block_ids 
#' @param blocking_dimension this dimension will be treated as experimental units or blocks (e.g. subject, session)
#' @param model_name the name of the statistical model to be applied to dataset.
#' @param model_control a list of control parameters for the model algorithm.
#' @export
PLS_Config <- function(
  base_path="",
  design_list,
  data_format="image",
  folding_dimension=NULL,
  blocking_dimension=NULL,
  block_ids=NULL,
  #stacking_dimension,
  model_name="Mean_Centered_PLS",
  model_control=list()) {
  obj <- list(
    data_list=design_list,
    design=design,
    data_format=data_format,
    #replicate_design=replicate_design,
    folding_dimension=folding_dimension,
    blocking_dimension=blocking_dimension,
    base_path=base_path,
    #block_ids=block_ids,
    #stacking_dimension=stacking_dimension,
    model_control=model_control)
  
  class(obj) <- c("list", "PLS_Config")
  obj
}

#' @export
configure_model <- function(conf) UseMethod("configure_model")

configure_model.PLSConfig <- function(conf) {
  subject_ids <- lapply(conf$design_list, "[[", "id")
  mask_names <- lapply(conf$design_list, "[[", "mask")
  
  design_list <- lapply(conf$design_list, function(dname) {
    read.table(dname$design, header=TRUE)
  })
    
  data_list <- lapply(seq_along(conf$design_list), function(i) {
    mask <- mask_names[[i]]
    maskvol <- neuroim::loadVolume(mask)
    lapply(data_files, function(fname) {
      neuroim::loadVector(fname, mask=mask)
    })
  })
  
  
      
}




