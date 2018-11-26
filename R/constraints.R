
#' @param coords the spatial coordinates as a matrix with rows as objects and columns as dimensions
#' @param nblocks the number of repeated coordinate blocks
#' @param sigma the standard deviation of the within-block smoother
#' @param shrinkage_factor 
spatial_constraints <- function(coords, nblocks=1, sigma=5, shrinkage_factor=.1, shrinkage_radius=1, 
                                global_weights=rep(1, ncol(coords)*nblocks)) {
  cds2 <- do.call(rbind, lapply(nblocks, function(i) cbind(coords, i*50)))
  
  ##Sw <- neighborweights::spatial_smoother(coords, sigma=sigma, nnk=27, stochastic = TRUE)
  Swithin <- neighborweights::spatial_smoother(cds2, sigma=sigma, nnk=27, stochastic = TRUE)
  #Matrix::bdiag(as.list(rep(Sw, nblocks))
  
  indices <- rep(1:nrow(coords), nblocks)
  
  Sbetween <- neighborweights::spatial_adjacency(as.matrix(indices),weight_mode="binary", 
                                                 normalize=FALSE,dthresh=1, include_diagonal=FALSE)/nblocks
  diag(Sbetween) <- 0
  
  Wg <- Diagonal(x=sqrt(global_weights))
  
  S <- Wg %*% Swithin %*% Wg
  S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
  S <- S + Sbetween*bw
  S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
  
 
}

# Swithin <- neighborweights::spatial_smoother(cds2,sigma=sigma,nnk=27,stochastic = TRUE)
# Sbetween <- neighborweights::spatial_adjacency(as.matrix(indices),weight_mode="binary", 
#                                                normalize=FALSE,dthresh=bdthresh, include_diagonal=FALSE)/length(perc_dat)
# diag(Sbetween) <- 0
# 
# Wg <- Diagonal(x=sqrt(wg))
# 
# S <- Wg %*% Swithin %*% Wg
# S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
# S <- S + Sbetween*bw
# S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
