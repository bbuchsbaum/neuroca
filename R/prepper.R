
#' add a pre-processing stage
#' @param x the processing pipeline
#' 
#' @export
add_node <- function(x,...) UseMethod("add_node")

#' prepare a dataset by applying a pre-processing pipeline
#' 
#' @param the pipeline
prep <- function(x, ...) UseMethod("prep")

#' contains a series of pre-processing steps
#' 
#' @keywords internal
prepper <- function() {
  steps <- list()
  ret <- list(steps=steps)
  class(ret) <- c("prepper", "list")
  ret
}

#' @export
add_node.prepper <- function(preproc, step) {
  preproc$steps[[length(preproc$steps)+1]] <- step
  preproc
}

#' @export
prep.prepper <- function(x, X) {
  ff=do.call(purrr::compose, lapply(x$steps, "[[", "forward"))
  Xp <- ff(X)
  
  
  ret <- list(
       preproc=x,
       Xp=Xp,
       transform=do.call(compose, lapply(x$steps, "[[", "apply")),
       reverse_transform=do.call(compose, rev(lapply(x$steps, "[[", "reverse"))))
  
  class(ret) <- "pre_processor"
  
}

#' @export
center <- function(preproc=prepper()) {
  env=new.env()
  
  forward <- function(X) {
    cmeans <- colMeans(X)
    env[["cmeans"]] <- cmeans
    sweep(X, 2, cmeans, "-")
  }
  
  apply <- function(X, colind=NULL) {
    cmeans <- env[["cmeans"]] 
    if (is.null(colind)) {
      sweep(X, 2, cmeans, "-")
    } else {
      assert_that(ncol(X) == length(colind))
      sweep(X, 2, cmeans[colind], "-")
    }
  }
    
  reverse=function(X, colind=NULL) {
    if (is.null(colind)) {
      sweep(X, 2, env[["cmeans"]], "+")
    } else {
      assert_that(ncol(X) == length(colind))
      sweep(X, 2, env[["cmeans"]][colind], "+")
    }
  }
  
  ret <- list(name="center",
              forward=forward,
              reverse=reverse,
              apply=apply)
  class(ret) <- c("center", "pre_processor")
  add_node(preproc, ret)
}

#' @export
colscale <- function(preproc=prepper(), type=c("unit", "z")) {
  type <- match.arg(type)
  env=new.env()
  
  forward <- function(X) {
    sds <- matrixStats::colSds(X)
    
    if (type == "unit") {
      sds <- sds * sqrt(nrow(X)-1)
    }
    
    env[["sds"]] <- sds
    sweep(X, 2, sds, "/")
  }
  
  apply <- function(X, colind=NULL) {
    if (is.null(colind)) {
      sweep(X, 2, env[["sids"]], "/")
    } else {
      assert_that(ncol(X) == length(colind))
      sweep(X, 2, env[["sids"]][colind], "/")
    }
  }
  
  reverse=function(X, colind=NULL) {
    if (is.null(colind)) {
      sweep(X, 2, env[["sds"]], "*")
    } else {
      assert_that(ncol(X) == length(colind))
      sweep(X, 2, env[["sds"]][colind], "*")
    }
  }
  
  ret <- list(forward=forward,
              reverse=reverse,
              apply=apply)
  
  class(ret) <- c("colscale", "pre_processor")
  add_node(preproc, ret)
}

#' dimension reduction as a pre-processing stage
#' @param preproc the pipeline
#' @param method the dimension reduction method (e.g. pca)
dim_reduce <- function(preproc=prepper(), method=pca, ...) {
  env=new.env()
  args <- list(...)
  
  forward <- function(X) {
    fit <- do.call(method, append(list(X), args))
    env[["fit"]] <- fit
    env[["X"]] <- X
    scores(fit)
  }
  
  apply <- function(X, colind=NULL) {
    if (is.null(colind)) {
      project(env[["fit"]], X)
    } else {
      assert_that(ncol(X) == length(colind))
      project(env[["fit"]], X, colind=colind)
    }
  }
  
  reverse=function(X, colind=NULL) {
    if (is.null(colind)) {
      env[["X"]]
    } else {
      assert_that(ncol(X) == length(colind))
      env[["X"]][,colind,drop=FALSE]
    }
  }
  
  ret <- list(name="dim_reduce",
              forward=forward,
              reverse=reverse,
              apply=apply)
  class(ret) <- c("dim_reduce", "pre_processor")
  add_node(preproc, ret)
}



# standardize <- function(f=identity) {
#   print(f)
#   force(f)
#   function(X) {
#     sweep(f(X), 2, apply(X,2,sd), "/")
#   }
# }


