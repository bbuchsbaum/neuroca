

#' pre_processor
#' 
#' @param X
#' @param scale
#' @param center
#' @export
pre_processor <- function(X, center=TRUE, scale=FALSE) {
  Xs <- scale(X, scale=scale, center=center)
  
  center_vec <- attr(Xs, "scaled:center")
  scale_vec <- attr(Xs, "scaled:scale")
  
  if (is.null(scale_vec)) {
    scale_vec <- rep(1, ncol(X))
  }
  
  if (is.null(center_vec)) {
    center_vec <- rep(0, ncol(X))
  }
  
  
  rf <- function(M) {
    if (!center && !scale) {
      M
    } else if (center && scale) {
      m1 <- sweep(M, 2, center_vec, "+")
      sweep(m1, 2, scale_vec, "*")
    } else if (!center && scale) {
      sweep(M, 2, scale_vec, "*")
    } else if (center && !scale) {
      sweep(M, 2, center_vec, "+")
    } else {
      stop()
    }
  }
  
  f <- function(M) {

    if (!center && !scale) {
      M
    } else if (center && scale) {
      m1 <- sweep(M, 2, center_vec, "-")
      sweep(m1, 2, scale_vec, "/")
    } else if (!center && scale) {
      sweep(M, 2, scale_vec, "/")
    } else if (center && !scale) {
      sweep(M, 2, center_vec, "-")
    } else {
      stop()
    }
  }
  
  attr(Xs, "pre_process") <- f
  attr(Xs, "reverse") <- rf
  attr(Xs, "center_vec") <- center_vec
  attr(Xs, "scale_vec") <- scale_vec
  Xs
}


