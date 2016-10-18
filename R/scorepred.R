#' @param fscores
#' @param scores
#' @param type
#' @param ncomp
#' @export
scorepred <- function(fscores, scores, type=c("class", "prob", "scores", "cosine", "distance"), ncomp=1) {
  if (type == "scores") {
    fscores
  } else if (type == "cosine") {
    proxy::simil(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE], method="cosine")
  } else if (type == "distance") {
    D <- fields::rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
  } else if (type =="class") {
    D <- rdist(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    D2 <- D^2
    min.d <- apply(D2, 1, which.min)
    row.names(scores)[min.d]
  } else if (type == "prob"){
    ## type is 'prob'
    probs <- tcrossprod(fscores[,1:ncomp,drop=FALSE], scores[,1:ncomp,drop=FALSE])
    maxid <- apply(probs,1,which.max)
    maxp <- probs[cbind(1:nrow(probs), maxid)]
    probs <- exp(probs - maxp)
    probs <- zapsmall(probs/rowSums(probs))
  } else {
    stop(paste("illegal 'type' argument: ", type))
  }
}
