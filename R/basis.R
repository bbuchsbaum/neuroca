

# gen_basis <- function(X, reducer=pca, ...) {
#   res <- reducer(as.matrix(X), ...)
#   ret <- list(projector=res, X=scores(res), Xorig=X, ncomp=ncomp(res))
#   class(ret) <- c("basis", "list")
#   ret
# }

