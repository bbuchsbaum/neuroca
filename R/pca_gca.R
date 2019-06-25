pca_gca <- function(X, ncomp=rep(2, length(X)), preproc=center(), 
                normalization=c("MFA", "RV", "None", "RV-MFA", "custom"), M=NULL, A=NULL, ...) {
  
  
  assertthat::assert_that(inherits(X, "block_matrix"), msg="X must be a 'block_matrix'")
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom") {
    assert_that(!is.null(A))
  }
  
  ## pre-process the variables.
  procres <- prep(preproc, X)
  Xp <- procres$Xp
  
  ## normalize the matrices 
  
  if (normalization != "custom") {
    alpha <- normalization_factors(Xp, type=normalization)
    A <- rep(alpha, block_lengths(X))
  } else {
    alpha <- rep(1, nblocks(X))
  }
  
  bind <- block_index_list(X)
  
  fit <- genpca(unclass(Xp), 
                preproc=pass(),
                A=A, 
                M=M,
                ncomp=ncomp,
                ...)
  
  
  result <- list(
    X=Xp,
    preproc=procres,
    ntables=nblocks(X),
    fit=fit,
    ncomp=fit$ncomp,
    block_indices=bind,
    alpha=alpha,
    normalization=normalization,
    table_names=names(X),
    nvars=ncol(X),
    ntables=length(block_lengths(X)),
    A=A,
    M=M
  )
  
  class(result) <- c("mfa", "multiblock", "bi_projector", "list")
  result
}
