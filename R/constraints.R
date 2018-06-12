#' compute a spatial contraint matrix that can span across blocks
#' 
#' @param a matrix of coordinates for one block that must replicate exactly across \code{nblocks} blocks
#' @param feature_scores a vector of feature weights used to weight each individual feature
#' @param sigma the within block smoothing radius
#' @param sigma_between the across block smoothing radius
#' @param shrinkage_factor the amount to weight the across block smoothing matrix (heiger leads to more shrinkae)
#' @export
spatial_constraints <- function(coords, nblocks=1, feature_scores=rep(1, nrow(coords)*nblocks), 
                                                                  sigma=5, sigma_between=5, 
                                                                  shrinkage_factor=.1) {
  
  Wg <- Diagonal(x=sqrt(feature_scores))   
  
  if (nblocks > 1 && sigma_between > 0) {
    print(paste("between smoothing, sigma = ", sigma_between))
    cds2 <- do.call(rbind, lapply(1:nblocks, function(i) cbind(coords, i*50)))
    ## within subject spatial smoother
    Swithin <- neighborweights::spatial_smoother(cds2,sigma=sigma,nnk=27,stochastic = TRUE)
    
    cds3 <- do.call(rbind, lapply(1:nblocks, function(i) coords))
    ## between subject spatial smoother
    Sbetween <- neighborweights::spatial_adjacency(cds3,weight_mode="binary", 
                                                 normalize=FALSE,
                                                 sigma=sigma_between,
                                                 dthresh=sigma_between*2,
                                                 nnk=27*nblocks)
  
    Sbetween <- Sbetween/(RSpectra::eigs_sym(Sbetween, k=1, which="LA")$values[1])
  
    
    S <- Wg %*% Swithin %*% Wg
    S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
    S <- S + Sbetween*shrinkage_factor
    S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
  
  } else {
    print("within smoothing only")
    Swithin <- Matrix::bdiag(replicate(nblocks, neighborweights::spatial_smoother(coords,sigma=sigma,nnk=27,stochastic = TRUE), simplify=FALSE))
    S <- Wg %*% Swithin %*% Wg
    S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
    S
  }
  
}