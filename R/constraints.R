
#' @param coords the spatial coordinates as a matrix with rows as objects and columns as dimensions
#' @param nblocks the number of repeated coordinate blocks
#' @param sigma the standard deviation of the within-block smoother
#' @param shrinkage_factor 
spatial_constraints <- function(coords, nblocks=1, 
                                sigma_within=5, 
                                sigma_between=1,
                                shrinkage_factor=.1, 
                                nnk_within=27,
                                nnk_between=1,
                                weight_mode_between="binary",
                                weight_mode_within="heat",
                                global_weights=rep(1, ncol(coords)*nblocks)) {
  
  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)
  
  #cds2 <- do.call(rbind, lapply(1:nblocks, function(i) cbind(coords)))
  coords <- as.matrix(coords)
  
  Sw <- neighborweights::spatial_adjacency(coords, sigma=sigma_within, 
                                           weight_mode=weight_mode_within,
                                           nnk=nnk_within, stochastic = TRUE)
  #Swithin <- neighborweights::spatial_smoother(cds2, sigma=sigma, nnk=50, stochastic = TRUE)
  Swithin <- Matrix::bdiag(replicate(nblocks, Sw, simplify=FALSE))
  #indices <- rep(1:nrow(coords), nblocks)
  
  Sb <- neighborweights::spatial_adjacency(coords,
                                                 weight_mode="binary", 
                                                 sigma=sigma_between,
                                                 normalize=FALSE,
                                                 nnk=nblocks*nnk_between,
                                                 dthresh=2.5,
                                                 include_diagonal=FALSE)
  Sbt <- as(Sb, "dgTMatrix")
  nvox <- nrow(coords)
  offsets <- cumsum(c(0, rep(nvox, nblocks-1)))
  
  out <- do.call(rbind, lapply(1:nblocks, function(i) {
    print(i)
    do.call(rbind, lapply(1:nblocks, function(j) {
      if (i == j) {
        NULL
      } else {
        cbind(Sbt@i + offsets[i] + 1, Sbt@j + offsets[j] +1, Sbt@x)
      }
    }))
  }))
  
  Sbfin <- sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(nvox*nblocks, nvox*nblocks))
  Sbfin <- neighborweights::make_doubly_stochastic(Sbfin)
  
 
  ## scale within matrix by variable weights
  if (any(global_weights[1] != global_weights)) {
    Wg <- Diagonal(x=sqrt(global_weights))
    Swithin <- Wg %*% Swithin %*% Wg
  }
  
  ## compute ration of witin to between weights
  rat <- sum(Swithin)/sum(Sbfin)
  sfac <- 1/shrinkage_factor * rat
  
  ## balance within and between weights
  Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin
  
  ## ratio of 
  
  S <- Stot/(RSpectra::eigs_sym(Stot, k=1, which="LA")$values[1])
  #S <- S + Sbetween*shrinkage_factor
  #S <- S/(RSpectra::eigs_sym(S, k=1, which="LA")$values[1])
  
 
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
