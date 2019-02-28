

#' construct a sparse matrix of spatial constraints for a set of blocks
#' 
#' @param coords the spatial coordinates as a matrix with rows as objects and columns as dimensions
#' @param nblocks the number of repeated coordinate blocks
#' @param sigma_within the bandwidth of the within-block smoother
#' @param sigma_between the bandwidth of the between-block smoother
#' @param shrinkage_factor the amount of shrinkage towards the spatially localized average
#' @param nnk_within the maximum number of nearest neighbors for within block smoother 
#' @param nnk_between the maximum number of nearest neighbors for between block smoother 
#' @param weight_mode_within
#' @param weight_mode_between
#' @param variable_weights
spatial_constraints <- function(coords, nblocks=1, 
                                sigma_within=5, 
                                sigma_between=1,
                                shrinkage_factor=.1, 
                                nnk_within=27,
                                nnk_between=1,
                                weight_mode_within="heat",
                                weight_mode_between="binary",
                                variable_weights=rep(1, ncol(coords)*nblocks)) {
  
  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)
  
 
  coords <- as.matrix(coords)
  
  Sw <- neighborweights::spatial_adjacency(coords, sigma=sigma_within, 
                                           weight_mode=weight_mode_within,
                                           nnk=nnk_within, stochastic = TRUE)
 
  Swithin <- Matrix::bdiag(replicate(nblocks, Sw, simplify=FALSE))
  #indices <- rep(1:nrow(coords), nblocks)
  
  Sb <- neighborweights::spatial_adjacency(coords,sigma=sigma_between,
                                                 weight_mode=weight_mode_between, 
                                                 nnk=nnk_between, stochastic=TRUE)
  
  browser()
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
  if (any(variable_weights[1] != variable_weights)) {
    Wg <- Diagonal(x=sqrt(variable_weights))
    Swithin <- Wg %*% Swithin %*% Wg
  }
  
  ## compute ratio of within to between weights
  rat <- sum(Swithin)/sum(Sbfin)
  sfac <- 1/shrinkage_factor * rat
  
  ## balance within and between weights
  Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin
  S <- Stot/(RSpectra::eigs_sym(Stot, k=1, which="LA")$values[1])
  S
 
}

