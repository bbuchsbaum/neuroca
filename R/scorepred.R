




#' @param x the model fit
#' @param labels the labels associated with the rows of the fitted or projected data (see newdata)
#' @param newdata an optional supplementary training set that is projected in to the fitted space.
#' @param colind the subset of column indices in the fitted model to use.
#' @param knn the number of nearest neighbors when classifying a new point. 
#' @export
classifier.projector <- function(x, labels, newdata=NULL, colind=NULL, knn=1) {

  
  if (is.null(newdata)) {
    newdata <- scores(x)
  } else {
    newdata <- project(x, newdata, colind=colind, ...)
  }
  
  assert_that(length(labels) == nrow(newdata))
  
  structure(
    list(
      fit=x,
      labels=labels,
      scores=newdata,
      colind=colind,
      knn=knn),
    class="classifier"
  )

}

rank_score <- function(prob, observed) {
  pnames <- colnames(prob)
  assert_that(all(observed %in% pnames))
  prank <- apply(prob, 1, function(p) {
    rp <- rank(p, ties.method="random")
    rp/length(rp)
  })
  
  mids <- match(observed, pnames)
  pp <- prank[cbind(mids, 1:length(observed))]
  
  data.frame(prank=pp, observed=observed)
}


#' @keywords internal
normalize_probs <- function(p) {
  apply(p, 2, function(v) {
    v2 <- v - min(v)
    v2/sum(v2)
  })
}

#' @keywords internal
avg_probs <- function(prob, labels) {
  pmeans <- t(group_means(labels, prob))
  t(apply(pmeans, 1, function(v) v/sum(v)))
}


#' @keywords internal
nearest_class <- function(prob, labels,knn=1) {
  
  apply(prob, 2, function(v) {
    ord <- order(v, decreasing=TRUE)[1:knn]
    l <- labels[ord]
    table(l)
    names(which.max(table(l)))
  })
  
}


#' @export
predict.classifier <- function(object, newdata, ncomp=object$fit$ncomp, metric=c("cosine", "euclidean"),
                               colind=NULL) {
  metric <- match.arg(metric)
  
  proj <- project(object$fit, newdata, comp=1:ncomp, colind=object$colind)
  
  doit <- function(p) {
    prob <- normalize_probs(p)
    pmeans <- avg_probs(prob, object$labels)
    cls <- nearest_class(prob, object$labels, object$knn)
    
    list(class=cls, prob=pmeans)
  }
    

  if (metric == "cosine") {
    p <- proxy::simil(as.matrix(object$scores)[,1:ncomp,drop=FALSE], as.matrix(proj), method="cosine")
    doit(p)

  } else if (metric == "euclidean") {
    D <- proxy::dist(as.matrix(object$scores)[,1:ncomp,drop=FALSE], as.matrix(proj), method="euclidean")
    doit(exp(-D))
  }
  
}


  
#' Given a set of projected scores and a set of reference scores, compute one of several performance metrics.
#' 
#' @param fscores the projected scores to be compared to the reference scores
#' @param scores the original scores from the fitted model
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