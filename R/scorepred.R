
#' Given a set of factor scores and a set of reference scores, compute one of several performance metrics.
#' 
#' @param fscores the projected scores to be compared to the reference scores
#' @param scores the reference scores
#' @param type the type of metric to compute
#' @param ncomp the number of dimensions to use
#' @export
scorepred <- function(fscores, scores, type=c("class", "prob", "scores", "cosine", "distance", "r2"), ncomp=ncol(fscores), class_name=TRUE) {
  if (type == "scores") {
    fscores
  } else if (type == "cosine") {
    proxy::simil(as.matrix(fscores[,1:ncomp,drop=FALSE]), as.matrix(scores[,1:ncomp,drop=FALSE]), method="cosine")
  } else if (type == "distance") {
    D <- fields::rdist(as.matrix(fscores[,1:ncomp,drop=FALSE]), as.matrix(scores[,1:ncomp,drop=FALSE]))
  } else if (type =="class") {
    D <- fields::rdist(as.matrix(fscores[,1:ncomp,drop=FALSE]), as.matrix(scores[,1:ncomp,drop=FALSE]))
    D2 <- D^2
    min.d <- apply(D2, 1, which.min)
    if (class_name) {
      row.names(scores)[min.d]
    } else {
      min.d
    }
  
  #} #else if (type == "r2") {
    ## distance of prediction to the category barycenters
    #DBetween <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    #DTotal <- (scores[, 1:ncomp] - rowMeans(scores[, 1:ncomp]))^2
    #DBetween/DTotal
  } else if (type == "prob"){
    D <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    Dmin <- D - min.col.value(D)
    dsd <- apply(Dmin, 1, sd)
    Dmin <- sweep(Dmin, 1, dsd, "/")
    probs <- exp(-Dmin)
    probs <- zapsmall(probs/rowSums(probs))
  } else {
    stop(paste("illegal 'type' argument: ", type))
  }
}

#' @keywords internal
max.col.value <- function (x) {
  return(x[cbind(1:nrow(x), max.col(x, ties.method = "first"))])
}

#' @keywords internal
min.col.value <- function (x) {
  return(x[cbind(1:nrow(x), max.col(-x, ties.method = "first"))])
}