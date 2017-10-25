
get_perm <- function(G, strata) {
  if (!is.null(x$strata)) {
    Gperm <- do.call(rbind, lapply(levels(x$strata), function(lev) {
      Gs <- G[strata==lev,]
      Gs <- Gs[sample(1:nrow(Gs)),]
    }))
  } else {
    Gperm <- x$G[sample(1:nrow(x$G)),]
  }
}



permutation_ <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL, deflate=FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  d <- singular_values(x)[1:ncomp]^2
  
  if (ncomp > x$ncomp) {
    ncomp <- x$ncomp
  }
  
  dstat0 <- matrix(0, nperms, ncomp)
  
  for (i in 1:nperms) {
    print(i)
    message("permutation ", i)
    fit <- x$permute_refit()
    dp <- singular_values(fit)[1:ncomp]^2
    dstat0[i, ] <- dp
  }
  
  p <- rep(0, ncomp)
  for (i in 1:ncomp) {
    p[i] <- mean(dstat0[, i] >= d[i])
  }
  
  for (i in 2:ncomp) {
    p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  return(list(perm_eig=dstat0, eig=d, pval=p))
  
}


permutation.musu_asca <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL,term="*") {
  if (term == "*") {
    term <- names(x$results)
  }
  
  res <- lapply(term, function(tname) {
    message("permutation for term: ", tname)
    permutation(x$results[[tname]]$bada_result, nperms=nperms, threshold=threshold,ncomp=ncomp, verbose=verbose, seed=seed)
  })
  
  names(res) <- term
  res
}


permutation.sca <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  permutation_(x,nperms,threshold, ncomp, verbose, seed)
}

permutation.mfa <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  permutation_(x,nperms,threshold, ncomp, verbose, seed)
}

permutation.mubada <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  permutation_(x,nperms,threshold, ncomp, verbose, seed)
}

permutation.bada <- function(x, nperms=100, threshold=.05, ncomp=1, verbose=TRUE, seed=NULL) {
  permutation_(x,nperms,threshold, ncomp, verbose, seed)
}

#' #' @export
#' permutation.pls_result_contrast <- function(x, nperms=100, threshold=.05, ncomp=2, verbose=TRUE, seed=NULL) {
#'   if (!is.null(seed)) set.seed(seed)
#'   
#'   dstat <- x$d[1:x$ncomp]^2/sum(x$d[1:x$ncomp]^2)
#'   
#'   dummy_matrix <- turner::factor_to_dummy(strata)
#'   
#'   subtract_recon <- function(recon) {
#'     ret <- lapply(1:length(Xblocks), function(i) {
#'       Xblocks[[i]] - recon
#'     })
#'   }
#'   
#'   permute_blocks <- function(xb, permset) {
#'     lapply(1:length(xb), function(i) {
#'       xb[[i]][permset[[i]],]
#'     })
#'   }
#'   
#'   Xblocks <- lapply(1:ncol(dummy_matrix), function(i) {
#'     ind <- dummy_matrix[,i] == 1
#'     scale(t(t(x$X[ind,]) %*% x$G[ind,]), center=x$center, scale=x$scale) 
#'   })
#'   
#'   if (ncomp > x$ncomp) {
#'     ncomp <- x$ncomp
#'   }
#'   
#'   dstat0 <- matrix(0, nperms, ncomp)
#'   
#'   for (i in 1:nperms) {
#'     
#'     message("permutation", i)
#'     for (j in 1:ncomp) {
#'       permset <- lapply(Xblocks, function(x) sample(1:nrow(x)))
#'       message("factor: ", j)
#'       if (j > 1) {
#'         message("recon")
#'         recon <- t(x$u[,1:(j-1),drop=FALSE] %*% t(x$scores[,1:(j-1),drop=FALSE]))
#'         message("subtract recon")
#'         Xresid <- subtract_recon(recon)
#'         message("permute blocks")
#'         Xperm <- permute_blocks(Xresid, permset)
#'         message("compute avg")
#'         Xpermavg <- Reduce("+", Xperm)/length(Xperm)
#'         #Xresidavg <- Reduce("+", Xresid)/length(Xresid)
#'         message("svd")
#'         dstat0[i,j] <- svd(Xpermavg)$d[1]
#'       } else {
#'         
#'         Xperm <- permute_blocks(Xblocks, permset)
#'         Xpermavg <- Reduce("+", Xperm)/length(Xperm)
#'         dstat0[i,j] <- svd(Xpermavg)$d[1]
#'       }
#'     }
#'   }
#'   
#'   p <- rep(0, x$ncomp)
#'   for (i in 1:x$ncomp) {
#'     p[i] <- mean(dstat0[, i] >= dstat[i])
#'   }
#'   for (i in 2:x$ncomp) {
#'     p[i] <- max(p[(i - 1)], p[i])
#'   }
#'   r <- sum(p <= threshold)
#'   return(list(r = r, p = p))
#'   
#' }

