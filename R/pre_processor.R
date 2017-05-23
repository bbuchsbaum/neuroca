

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
      sweep(M, 2, center_vec, "-")
    } else if (!center && scale) {
      sweep(M, 2, scale_vec, "/")
    } else if (center && scale) {
      M <- sweep(M, 2, center_vec)
      sweep(M, 2, scale_vec, "/")
    } else {
      stop()
    }
  }
  
  attr(Xs, "pre_process") <- f
  Xs
}


