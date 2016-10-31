#' @param fscores
#' @param scores
#' @param type
#' @param ncomp
#' @export
scorepred <- function(fscores, scores, type=c("class", "prob", "scores", "cosine", "distance", "r2"), ncomp=1) {
  if (type == "scores") {
    fscores
  } else if (type == "cosine") {
    proxy::simil(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE], method="cosine")
  } else if (type == "distance") {
    D <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
  } else if (type =="class") {
    D <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    D2 <- D^2
    min.d <- apply(D2, 1, which.min)
    row.names(scores)[min.d]
  #} #else if (type == "r2") {
    ## distance of prediction to the category barycenters
    #DBetween <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    #DTotal <- (scores[, 1:ncomp] - rowMeans(scores[, 1:ncomp]))^2
    #DBetween/DTotal
  } else if (type == "prob"){
    D <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])^2
    t(apply(D, 1, function(x) 1/x / sum(1/x)))
  } else {
    stop(paste("illegal 'type' argument: ", type))
  }
}
