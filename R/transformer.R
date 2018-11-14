
add <- function(x,...) UseMethod("add")

prep <- function(x, ...) UseMethod("prep")

pipeline <- function() {
  steps <- list()
  ret <- list(steps=steps)
  class(ret) <- c("pipeline", "list")
  ret
}

add.pipeline <- function(preproc, step) {
  preproc$steps[[length(preproc$steps)+1]] <- step
  preproc
}

prep.pipeline <- function(x, X) {
  ff=do.call(purrr::compose, lapply(x$steps, "[[", "forward"))
  Xp <- ff(X)
  
  list(Xp=Xp,
       transform=do.call(compose, lapply(x$steps, "[[", "apply")),
       reverse_transform=do.call(compose, rev(lapply(x$steps, "[[", "reverse"))))
  
}

transform.pipeline <- function() {
  lapply()
}

center <- function(preproc=pipeline()) {
  env=new.env()
  
  forward <- function(X) {
    cmeans <- colMeans(X)
    env[["cmeans"]] <- cmeans
    sweep(X, 2, cmeans, "-")
  }
  
  apply <- function(X, colind=NULL) {
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
  add(preproc, ret)
}


colscale <- function(preproc=pipeline()) {
  env=new.env()
  
  forward <- function(X) {
    sds <- matrixStats::colSds(X)
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
  add(preproc, ret)
}

dim_reduce <- function(preproc=pipeline(), method=pca, ...) {
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
  add(preproc, ret)
}



# standardize <- function(f=identity) {
#   print(f)
#   force(f)
#   function(X) {
#     sweep(f(X), 2, apply(X,2,sd), "/")
#   }
# }


