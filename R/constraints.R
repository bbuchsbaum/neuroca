
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
                                           nnk=nnk_within,  normalized=FALSE,
                                           stochastic = TRUE)
 
  if (nblocks == 1) {
    return(Sw)
  }
  
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
                                                 normalized=FALSE,
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


#' @export
#' @examples 
#' 
#' coords <- as.matrix(expand.grid(1:10, 1:10))
#' fmats <- replicate(20, matrix(rnorm(100*10), 10, 100), simplify=FALSE)
#' conmat <- feature_weighted_spatial_constraints(coords, fmats)
#' 
#' conmat <- feature_weighted_spatial_constraints(coords, fmats, maxk_between=4, maxk_within=2,sigma_between=5, nnk_between=60)
#' 
#' #future::plan("multiprocess")
#' #conmat <- feature_weighted_spatial_constraints(coords, fmats)
#' @importFrom furrr future_map
feature_weighted_spatial_constraints <- function(coords, 
                                                 feature_mats,
                                                 sigma_within=5, 
                                                 sigma_between=3,
                                                 wsigma_within=.73,
                                                 wsigma_between=.73,
                                                 alpha_within=.5,
                                                 alpha_between=.5,
                                                 shrinkage_factor=.1, 
                                                 nnk_within=27,
                                                 nnk_between=27,
                                                 maxk_within=nnk_within,
                                                 maxk_between=nnk_between,
                                                 weight_mode_within="heat",
                                                 weight_mode_between="binary",
                                                 variable_weights=rep(1, ncol(coords)*length(feature_mats)), verbose=FALSE) {
  
  assert_that(shrinkage_factor > 0 & shrinkage_factor <= 1)

  coords <- as.matrix(coords)
  nvox <- nrow(coords)
  nblocks <- length(feature_mats)

  Swl <- furrr::future_map(seq_along(feature_mats), function(i) {
    sw <- neighborweights::weighted_spatial_adjacency(coords, t(feature_mats[[i]]), 
                                           alpha=alpha_within,
                                           wsigma=wsigma_within,
                                           sigma=sigma_within, 
                                           weight_mode=weight_mode_within,
                                           nnk=nnk_within, normalized=TRUE, stochastic=FALSE)
    neighborweights::make_doubly_stochastic(sw)
  })
  
 
  Swithin <- Matrix::bdiag(Swl)
  
  cmb <- t(combn(1:length(feature_mats), 2))

  offsets <- cumsum(c(0, rep(nvox, nblocks-1)))
  
  
  bet <- do.call(rbind, furrr::future_map(1:nrow(cmb), function(i) {
    a <- cmb[i,1]
    b <- cmb[i,2]
    sm <- neighborweights::cross_weighted_spatial_adjacency(
      coords, coords, t(feature_mats[[a]]), t(feature_mats[[b]]),
      wsigma=wsigma_between, weight_mode=weight_mode_between,
      alpha=alpha_between,
      nnk=nnk_between, maxk=maxk_between,
      sigma=sigma_between, normalized=FALSE)
    sm <- neighborweights::make_doubly_stochastic(sm)
    
    sm_nc <- as (sm, "dgTMatrix")
    r1 <- cbind (i = sm_nc@i + 1 + offsets[a], j = sm_nc@j + 1 + offsets[b], x = sm_nc@x)
    r2 <- cbind (i = sm_nc@j + 1 + offsets[b], j = sm_nc@i + 1 + offsets[a], x = sm_nc@x)
    rbind(r1,r2)
  }))
  

  Sbfin <- sparseMatrix(i=bet[,1], j=bet[,2], x=bet[,3], dims=c(nvox*nblocks, nvox*nblocks))

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

