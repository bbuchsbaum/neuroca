



#' scada
#' 
#' simulataneous component discriminant analysis
#' 
#' @importFrom assertthat assert_that 
#' 
#' @param Y dependent \code{factor} variable. If All X matrices have same number of rows, Y can be a single factor.
#'        If there are a different number of rows (e.g. different numbers of replications per subject), Y can be a list of factors.
#' @param Xlist a list of X matrices, one per subject. 
#' @param ncomp number of common components to estimate
#' @param center whether to center the variables
#' @param scale whether to scale the variables by 1/sd
#' @param normalization the type of normalization
#' @param rank_k reduce data to k components per block via pca
#' @export
scada <- function(Y, Xlist, ncomp=2, center=TRUE, scale=FALSE,rank_k=NULL,
                  type=c("sca-p","sca-pf2","sca-ind","sca-ecp")) {


  type <- match.arg(type)

  assertthat::assert_that(all(sapply(Xlist, is.matrix)))
  
  cstruc <- prep_multiblock_da(Y, Xlist)
  Yl <- cstruc$Yl
  Xlist <- cstruc$Xlist
  conditions <- cstruc$conditions
  Xr <- cstruc$Xr

  fit <- sca(Xr, ncomp=ncomp, center=center, scale=scale, type=type, rank_k=rank_k)
 
  result <- list(
    Xlist=Xlist,
    Y=Yl,
    Y_reps=cstruc$Y_reps,
    conditions=conditions,
    Xr=Xr,
    scores=scores(fit),
    ntables=length(Xlist),
    ncond=nrow(Xr),
    fit=fit,
    center=center,
    scale=scale,
    ncomp=fit$ncomp,
    block_indices=fit$block_indices,
    refit=refit,
    table_names=cstruc$table_names,
    reprocess=fit$reprocess,
    rank_k=rank_k,
    type=type,
    permute_refit=permute_refit
  )
  
  class(result) <- c("scada", "multiblock_da", "list")
  result
}


#' @export
refit.scada <- function(x, Y, Xlist, ncomp=x$ncomp) { 
  scada(Y, Xlist, ncomp=ncomp, center=x$center, scale=x$scale, type=x$type, rank_k=x$rank_k) 
}



