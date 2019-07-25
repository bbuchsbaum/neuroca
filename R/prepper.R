#' get a fresh pre-processing node cleared of any cached data
#' @param x the processing pipeline
#' 
#' 
fresh <- function(x,...) UseMethod("fresh")



#' add a pre-processing stage
#' @param x the processing pipeline
#' 
#' @export
add_node <- function(x,...) UseMethod("add_node")

#' prepare a dataset by applying a pre-processing pipeline
#' 
#' @param the pipeline
#' @param ... extra args
#' @export
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


#' @importFrom purrr compose
#' @export
prep.prepper <- function(x, X) {
  ff=do.call(purrr::compose, lapply(x$steps, "[[", "forward"))
  Xp <- ff(X)
  
  tform <- function(X, colind=NULL) {
    xin <- X
    for (i in 1:length(x$steps)) {
      xin <- x$steps[[i]]$apply(xin, colind)
    }
    
    xin
  }
  
  rtform <- function(X, colind=NULL) {
    xin <- X
    for (i in length(x$steps):1) {
      xin <- x$steps[[i]]$reverse(xin, colind)
    }
    
    xin
  }
  ret <- list(
       preproc=x,
       Xp=Xp,
       transform=tform,
       reverse_transform=rtform)
       #transform=do.call(purrr::compose, lapply(x$steps, "[[", "apply"), .dir="forward"),
       #reverse_transform=do.call(purrr::compose, lapply(x$steps, "[[", "reverse"), .dir="backward"))
  
  class(ret) <- "pre_processor"
  ret
  
}

fresh.prepper <- function(x) {
  p <- prepper()
  for (step in x$steps) {
    p <- prep_node(p, step$name, step$create)
  }
  p
}

fresh.pre_processor <- function(x, preproc=prepper()) {
  p <- x$create()
}

prep_node <- function(pipeline, name, create,  ...) {
  node <- create()
  ret <- list(name=name,
              create=create,
              forward=node$forward,
              reverse=node$reverse,
              apply=node$apply,
              ...)
  class(ret) <- c(name, "pre_processor")
  add_node(pipeline, ret)
}

#' @export
pass <- function(preproc=prepper()) {
  
  create <- function() {
    list(
      forward = function(X, colind=NULL) {
        X
      },
  
      reverse = function(X, colind=NULL) {
        X
      },
  
      apply = function(X, colind=NULL) {
        X
      }
    )
  }
  
  prep_node(preproc, "pass", create)
  
}

#' @export
center <- function(preproc = prepper()) {
  create <- function() {
    env = new.env()
    
    list(
      forward = function(X) {
        cmeans <- colMeans(X)
        env[["cmeans"]] <- cmeans
        sweep(X, 2, cmeans, "-")
      },
      
      apply = function(X, colind = NULL) {
        cmeans <- env[["cmeans"]]
        if (is.null(colind)) {
          sweep(X, 2, cmeans, "-")
        } else {
          assert_that(ncol(X) == length(colind))
          sweep(X, 2, cmeans[colind], "-")
        }
      },
      
      reverse = function(X, colind = NULL) {
        if (is.null(colind)) {
          sweep(X, 2, env[["cmeans"]], "+")
        } else {
          assert_that(ncol(X) == length(colind))
          sweep(X, 2, env[["cmeans"]][colind], "+")
        }
      }
    )
  }
  
  prep_node(preproc, "center", create)
}

#' @export
colscale <- function(preproc = prepper(),
           type = c("unit", "z", "weights"),
           weights = NULL) {
    type <- match.arg(type)
    
    if (type != "weights" && !is.null(weights)) {
      warning("colscale: weight ignored because type != 'weights'")
    }
    if (type == "weights") {
      assert_that(!is.null(weights))
    }
    
    create <- function() {
      env = new.env()
      list(
        forward = function(X) {
          wts <- if (type == "weights") {
            assert_that(length(weights) == ncol(X))
            weights
          } else {
            sds <- matrixStats::colSds(X)
            
            if (type == "unit") {
              sds <- sds * sqrt(nrow(X) - 1)
            }
            
            sds[sds == 0] <- median(sds)
            1 / sds
          }
          env[["weights"]] <- wts
          sweep(X, 2, wts, "*")
          
        },
        
        apply = function(X, colind = NULL) {
          if (is.null(colind)) {
            sweep(X, 2, env[["weights"]], "*")
          } else {
            assert_that(ncol(X) == length(colind))
            sweep(X, 2, env[["weights"]][colind], "*")
          }
        },
        
        reverse = function(X, colind = NULL) {
          if (is.null(colind)) {
            sweep(X, 2, env[["weights"]], "/")
          } else {
            assert_that(ncol(X) == length(colind))
            sweep(X, 2, env[["weights"]][colind], "/")
          }
        }
      )
    }
    
    prep_node(preproc, "colscale", create)
  }

#' dimension reduction as a pre-processing stage
#' @param preproc the pipeline
#' @param method the dimension reduction method (e.g. pca)
#' @param ... args to be passed to the dimension reduction method
#' @export
dim_reduce <- function(preproc = prepper(),
                       method = pca,
                       ...) {
  args <- list(...)
  
  create = function() {
    env = new.env()
    
    list(
      forward = function(X) {
        fit <- do.call(method, append(list(X), args))
        env[["fit"]] <- fit
        env[["X"]] <- X
        scores(fit)
      },
      
      apply = function(X, colind = NULL) {
        if (is.null(colind)) {
          project(env[["fit"]], X)
        } else {
          assert_that(ncol(X) == length(colind))
          project(env[["fit"]], X, colind = colind)
        }
      },
      
      reverse = function(X, colind = NULL) {
        if (is.null(colind)) {
          env[["X"]]
        } else {
          assert_that(ncol(X) == length(colind))
          env[["X"]][, colind, drop = FALSE]
        }
      }
    )
  }
  
  prep_node(preproc, "dim_reduce", create)
}




#' @export
standardize <- function(preproc = prepper()) {
  create <- function() {
    env = new.env()
    
    list(
      forward = function(X) {
        sds <- matrixStats::colSds(X)
        sds[sds == 0] <- median(sds)
        cmeans <- colMeans(X)
        env[["sds"]] <- sds
        env[["cmeans"]] <- cmeans
        x1 <- sweep(X, 2, cmeans, "-")
        sweep(x1, 2, sds, "/")
      },
      
      apply = function(X, colind = NULL) {
        if (is.null(colind)) {
          x1 <- sweep(X, 2, env[["cmeans"]], "-")
          sweep(x1, 2, env[["sds"]], "/")
        } else {
          assert_that(ncol(X) == length(colind))
          x1 <- sweep(X, 2, env[["cmeans"]][colind], "-")
          sweep(x1, 2, env[["sds"]][colind], "/")
        }
      },
      
      reverse = function(X, colind = NULL) {
        if (is.null(colind)) {
          x0 <- sweep(X, 2, env[["sds"]], "*")
          sweep(x0, 2, env[["cmeans"]], "+")
        } else {
          assert_that(ncol(X) == length(colind))
          x0 <- sweep(X, 2, env[["sds"]][colind], "*")
          sweep(x0, 2, env[["cmeans"]][colind], "+")
        }
      }
    )
  }
  prep_node(preproc, "standardize", create)
}

