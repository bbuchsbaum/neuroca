#' @export
#' @param Y \code{factor} variable defining the groups
#' @param X \code{matrix} defining the matrix data to be group-wise averaged
group_means <- function(Y, X) {
  if (all(table(Y) == 1)) {
    row.names(X) <- names(table(Y))
    X
  } else {
    Rs <- rowsum(X,Y)
    yt <- table(Y)
    ret <- sweep(Rs, 1, yt, "/")
    row.names(ret) <- names(yt)
    ret
  }
}

#group_means1 <- function(Y, X) {
# G <- model.matrix(~ Y - 1)
# colnames(G) <- levels(Y)
# GW <- G/colSums(G)
# R <- t(crossprod(X, GW))
#centroid <-  rowMeans(R)
#Rcent <- sweep(R, 1, centroid)
#list(Rcent=Rcent, Ymat=G, centroid=centroid)
#}\\



