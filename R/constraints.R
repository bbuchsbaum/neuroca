
#' spatial constraints
#' 
#' 
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
                                variable_weights=rep(1, ncol(coords)*nblocks), verbose=FALSE) {
  
  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)
  
 
  coords <- as.matrix(coords)
  
  if (verbose) {
    message("spatial_contraints: computing spatial adjacency (within)")
  }
  
  Sw <- neighborweights::spatial_adjacency(coords, sigma=sigma_within, 
                                           weight_mode=weight_mode_within,
                                           nnk=nnk_within, stochastic = TRUE)
 
  if (verbose) {
    message("spatial_contraints: replicating within blocks")
  }
  Swithin <- Matrix::bdiag(replicate(nblocks, Sw, simplify=FALSE))
  #indices <- rep(1:nrow(coords), nblocks)
  
  if (verbose) {
    message("spatial_contraints: computing spatial adjacency (between)")
  }
  
  Sb <- neighborweights::spatial_adjacency(coords,sigma=sigma_between,
                                                 weight_mode=weight_mode_between, 
                                                 nnk=nnk_between, stochastic=TRUE)
  
 
  Sbt <- as(Sb, "dgTMatrix")
  nvox <- nrow(coords)
  offsets <- cumsum(c(0, rep(nvox, nblocks-1)))
  
  out <- unlist(lapply(1:nblocks, function(i) {
    print(i)
    lapply(1:nblocks, function(j) {
      if (i == j) {
        NULL
      } else {
        cbind(Sbt@i + offsets[i] + 1, Sbt@j + offsets[j] +1, Sbt@x)
      }
    })
  }), recursive=FALSE)
  
  out <- out[!sapply(out, is.null)]
  out <- do.call(rbind, out)
  
  
  if (verbose) {
    message("spatial_contraints: building between matrix")
  }
  
  
  Sbfin <- sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(nvox*nblocks, nvox*nblocks))
  
  if (verbose) {
    message("spatial_contraints: making doubly stochastic")
  }
  
  Sbfin <- neighborweights::make_doubly_stochastic(Sbfin)
  
 
  ## scale within matrix by variable weights
  if (any(variable_weights[1] != variable_weights)) {
    Wg <- Diagonal(x=sqrt(variable_weights))
    Swithin <- Wg %*% Swithin %*% Wg
  }
  
  if (verbose) {
    message("spatial_contraints: weighting by sfac")
  }
  
  
  ## compute ratio of within to between weights
  rat <- sum(Swithin)/sum(Sbfin)
  sfac <- 1/shrinkage_factor * rat
  
  if (verbose) {
    message("spatial_contraints: balancing within and between by shrinkage_factor")
  }
  
  ## balance within and between weights
  Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin
  
  if (verbose) {
    message("spatial_contraints: normalizing by first eigenvalue")
  }
  
  S <- Stot/(RSpectra::eigs_sym(Stot, k=1, which="LA")$values[1])
  S
 
}

feature_weighted_spatial_constraints <- function(coords, 
                                                 feature_mats,
                                                 nblocks=length(feature_mats), 
                                                 sigma_within=5, 
                                                 alpha=.5,
                                                 sigma_between=3,
                                                 wsigma_within=.73,
                                                 wsigma_between=.73,
                                                 shrinkage_factor=.1, 
                                                 nnk_within=27,
                                                 nnk_between=5,
                                                 weight_mode_within="heat",
                                                 weight_mode_between="binary",
                                                 variable_weights=rep(1, ncol(coords)*nblocks), verbose=FALSE) {
  
  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)
  
  
  coords <- as.matrix(coords)
  
  if (verbose) {
    message("spatial_contraints: computing spatial adjacency (within)")
  }
  
  Sw <- neighborweights::weighted_spatial_adjacency(coords, feature_mats[[1]], 
                                           wsigma=wsigma_within,
                                           sigma=sigma_within, 
                                           weight_mode=weight_mode_within,
                                           nnk=nnk_within, stochastic = TRUE)
  
  if (verbose) {
    message("spatial_contraints: replicating within blocks")
  }
  Swithin <- Matrix::bdiag(replicate(nblocks, Sw, simplify=FALSE))
  #indices <- rep(1:nrow(coords), nblocks)
  
  if (verbose) {
    message("spatial_contraints: computing spatial adjacency (between)")
  }
  
  Sb <- neighborweights::spatial_adjacency(coords,sigma=sigma_between,
                                           weight_mode=weight_mode_between, 
                                           nnk=nnk_between, stochastic=TRUE)
  
  
  Sbt <- as(Sb, "dgTMatrix")
  nvox <- nrow(coords)
  offsets <- cumsum(c(0, rep(nvox, nblocks-1)))
  
  
  out <- unlist(furrr::future_map(1:nblocks, function(i) {
    print(i)
    lapply(1:nblocks, function(j) {
      if (i == j) {
        NULL
      } else {
        cbind(Sbt@i + offsets[i] + 1, Sbt@j + offsets[j] +1, Sbt@x)
      }
    })
  }), recursive=FALSE)
  
  out <- out[!sapply(out, is.null)]
  out <- do.call(rbind, out)
  
  
  if (verbose) {
    message("spatial_contraints: building between matrix")
  }
  
  
  Sbfin <- sparseMatrix(i=out[,1], j=out[,2], x=out[,3], dims=c(nvox*nblocks, nvox*nblocks))
  
  if (verbose) {
    message("spatial_contraints: making doubly stochastic")
  }
  
  Sbfin <- neighborweights::make_doubly_stochastic(Sbfin)
  
  
  ## scale within matrix by variable weights
  if (any(variable_weights[1] != variable_weights)) {
    Wg <- Diagonal(x=sqrt(variable_weights))
    Swithin <- Wg %*% Swithin %*% Wg
  }
  
  if (verbose) {
    message("spatial_contraints: weighting by sfac")
  }
  
  
  ## compute ratio of within to between weights
  rat <- sum(Swithin)/sum(Sbfin)
  sfac <- 1/shrinkage_factor * rat
  
  if (verbose) {
    message("spatial_contraints: balancing within and between by shrinkage_factor")
  }
  
  ## balance within and between weights
  Stot <- (1-shrinkage_factor)*(1/rat)*Swithin + shrinkage_factor*Sbfin
  
  if (verbose) {
    message("spatial_contraints: normalizing by first eigenvalue")
  }
  
  S <- Stot/(RSpectra::eigs_sym(Stot, k=1, which="LA")$values[1])
  S
  
}

