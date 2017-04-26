

#' pre_processor
#' 
#' @param X
#' @param scale
#' @param center
#' @export
pre_processor <- function(X, scale=FALSE, center=TRUE) {
  Xs <- scale(X, scale=scale, center=center)
  
  center_vec <- attr(Xs, "scaled:center")
  scale_vec <- attr(Xs, "scaled:scale")
  
  f <- function(M) {
    if (!center && !scale) {
      M
    } else if (center && scale) {
      sweep(M, 2, cen, "-")
    } else if (!center && scale) {
      sweep(M, 2, stdev, "/")
    } else if (center && scale) {
      M <- sweep(M, 2, cen)
      sweep(M, 2, stdev, "/")
    } else {
      stop()
    }
  }
  
  attr(Xs, "pre_process") <- f
  Xs
}


