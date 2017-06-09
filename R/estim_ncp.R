
#' fast_estim_ncomp
#' 
#' @param X
#' @param ncp.min
#' @param ncp.max
#' @param scale
#' @export
fast_estim_ncomp <- function(X, ncp.min=0,ncp.max=NULL, scale=FALSE) {
  
  p=ncol(X)
  n=nrow(X)
  
  if (is.null(ncp.max)) 
    ncp.max <- ncol(X)
  
  ncp.max <- min(nrow(X),ncol(X),ncp.max)
  
  crit <- NULL
  
  X <- scale(X,scale)
  
  if (ncp.min==0)  crit0 = mean(X^2, na.rm = TRUE)*(n*p)/(p*(n-1))
  
  rr <- svd(X)
  
  chash <- list()
  chash[["0"]] <- crit0
  
  compute_gcv <- function(q) {
    rec <- tcrossprod(sweep(rr$u[, 1:q, drop = F], 2, rr$d[1:q], FUN = "*"), rr$v[, 1:q, drop = F])
    mean(( n*p*(X-rec)/ ((n-1)*(p)- q*(n+p-q-1)))^2, na.rm=T)
  }
  
  run_gcv <- function(cand) {
    crit <- NULL
    for (i in seq_along(cand)) {
      q <- cand[i]
      message("q: ", q)
      
      if (!is.null(chash[[as.character(q)]])) {
        cval <- chash[[as.character(q)]]
      } else  {
        cval <- compute_gcv(q)
        chash[[as.character(q)]] <<- cval
      } 
      
      crit <- c(crit, cval)
      if (i > 1 && crit[i] > crit[i-1]) {
        break
      }
    }
    
    cbind(cand[1:length(crit)], crit)
  }
  
  if (ncp.max < 10) {
    cand <- seq(ncp.min, ncp.max)
  } else {
    cand <- sort(unique(round(c(seq(ncp.min, ncp.max, by=sqrt(ncp.max)), ncp.max))))
  }
  
  if (length(cand) < 10) {
    for (i in cand) {
      if (i == 0) {
        next
      }
      chash[[as.character(i)]] <- compute_gcv(i)
    }
  } else {
    
    while (TRUE) {
      crit <- run_gcv(cand)
      wmin <- which.min(crit[,2])
      opt <-crit[wmin,1]
      
      low <- if (wmin == 1) {
        opt
      } else {
        max(ncp.min, crit[wmin-1,1] + 1)
      }
      
      
      high <- if (wmin == nrow(crit)) {
        opt
      } else {
        min(ncp.max, crit[wmin+1,1] - 1)
      }
      
      
      if ((high - low) < 5) { 
        cand <- seq(low, high, by=1)
        crit <- run_gcv(cand)
        wmin <- which.min(crit[,2])
        opt <- crit[wmin,1]
        break
      } else if (low == high) {
        break
      } else {
        cand <- sort(unique(c(round(seq(low, high, by=sqrt(high-low))), high)))
      }
    }
  }
  
  ord <- order(as.numeric(names(chash)))
  list(ncomp=as.integer(names(chash)[ord]), criterion=as.numeric(chash[ord]), bestcomp=opt)
  
}
