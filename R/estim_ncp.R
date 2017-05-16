
#' fast_estim_ncomp
#' 
#' @param X
#' @param ncp.min
#' @param ncp.max
#' @param scale
#' @export
fast_estim_ncomp <- function(X, ncp.min=0,ncp.max=NULL, scale=TRUE){
  
  p=ncol(X)
  n=nrow(X)
  
  if (is.null(ncp.max)) 
    ncp.max <- ncol(X)
  
  ncp.max <- min(nrow(X)-2,ncol(X)-1,ncp.max)
  
  crit <- NULL
  
  X <- scale(X,scale)
  
  if (ncp.min==0)  crit = mean(X^2, na.rm = TRUE)*(n*p)/(p*(n-1))
  
  rr <- RSpectra::svds(X,k=ncp.max)
  
  chash <- list()
  
  run_gcv <- function(cand) {
    crit <- NULL
    for (i in seq_along(cand)) {
      q <- cand[i]
      message("q: ", q)
      
      if (!is.null(chash[[as.character(q)]])) {
        cval <- chash[[as.character(q)]]
      } else if (q>1)  {
        rec <- tcrossprod(sweep(rr$u[, 1:q, drop = F], 2, rr$d[1:q], FUN = "*"), rr$v[, 1:q, drop = F])
        cval <- mean(( n*p*(X-rec)/ ((n-1)*(p)- q*(n+p-q-1)))^2, na.rm=T)
        chash[[as.character(q)]] <<- cval
        print(chash)
      } else {
        rec <- tcrossprod(rr$u[,1,drop=FALSE]*rr$d[1],rr$v[,1,drop=FALSE])
        cval <- mean(( n*p*(X-rec)/ ((n-1)*(p)- q*(n+p-q-1)))^2, na.rm=T)
        chash[[as.character(q)]] <<- cval
      }
      
      crit <- c(crit, cval)
      if (i > 1 && crit[i] > crit[i-1]) {
          break
      }
    }
    
    cbind(cand[1:length(crit)], crit)
  }
  
  cand <- sort(unique(round(c(seq(ncp.min, ncp.max, by=sqrt(ncp.max)), ncp.max))))
  
  while (TRUE) {
    crit <- run_gcv(cand)
    wmin <- which.min(crit[,2])
    low <- max(ncp.min, crit[wmin-1,1] + 1)
    high <- min(ncp.max, crit[wmin+1,1] - 1)
    opt <-crit[wmin,1]
    
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
  
  ord <- order(as.numeric(names(chash)))
  list(ncomp=as.integer(names(chash)[ord]), criterion=as.numeric(chash[ord]), bestcomp=opt)
                    
}
