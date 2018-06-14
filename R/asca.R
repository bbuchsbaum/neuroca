
gen_threeway_mat <- function(n1,n2,n3, nvox, index) {
  nk <- n1*n2*n3

  grid <- expand.grid(A=paste0("A", 1:n1), B=paste0("B", 1:n2), C=paste0("C", 1:n3))
  grid$sid <- index
  mat <- matrix(rnorm(nk*nvox), nk, nvox)
  list(mat=mat, design=grid)

}

#' @keywords internal
construct_design_table <- function(termorder, facs, Ymaximal, all_terms, designL) {
  ## dgrid contains the factors for each term in the ANOVA
  dgrid <- expand.grid(lapply(all_terms, function(tname) levels(designL[[1]][[tname]])))
  names(dgrid) <- all_terms
  
  ## add higher-order terms to dgrid
  if (any(termorder > 1)) {
    high_facs <- which(termorder > 1)
    
    for (i in high_facs) {
      nam <- colnames(facs)[i]
      idx <- which(facs[,i] > 0)
      cols <- dgrid[,idx]
      dgrid[[nam]] <- do.call(paste, c(cols, list(sep=":")))
    }
  }
  
  
  dgrid[[paste0(all_terms, collapse=":")]] <- levels(Ymaximal[[1]])
  dgrid
}


#' Column-wise average a matrix of variables, collapsing over some set of factors
#' 
#' @param form the formula defining the model to collapse over
#' @param X the matrix to collpased over rows
#' @param design the \code{data.frame} containing the variables reference in \code{form}
#' @export
#' @examples 
#' X <- matrix(rnorm(20*10), 20, 10)
#' des <- data.frame(a=rep(letters[1:4], 5), b=factor(rep(1:5, each=4)))
#' xcoll <- collapse(~ a+b, X, design=des)
collapse <- function(form, X, design) {
  tt <- terms(form)
  facnames <- colnames(attr(tt, "factors"))
  facs <- lapply(facnames, function(fname) as.factor(design[[fname]]))
  names(facs) <- facnames
  crossed <- do.call(interaction, facs)
  Xbar <- group_means(crossed, X)
  
  levlist <- lapply(facs, levels)
  rdes <- do.call(expand.grid, levlist)
  list(X=Xbar, design=rdes)
}


#' Compute a regression model for each column in a matrix and return residual matrix
#' 
#' @param form the formula defining the model to fit for residuals
#' @param X the response matrix
#' @param design the \code{data.frame} containing the design variables specified in \code{form} argument.
#' @examples 
#' 
#' X <- matrix(rnorm(20*10), 20, 10)
#' des <- data.frame(a=rep(letters[1:4], 5), b=factor(rep(1:5, each=4)))
#' xresid <- residualize(~ a+b, X, design=des)
#' 
#' ## design is saturated, residuals should be zero
#' xresid2 <- residualize(~ a*b, X, design=des)
#' sum(xresid2) == 0
#' @export
residualize <- function(form, X, design) {
  options(contrasts = c("contr.sum", "contr.poly"))
  modmat <- model.matrix(form, data=design)
  resid(lsfit(modmat, X, intercept=FALSE))
}
  

#' muasca
#' 
#' @param formula a formula specifying the ANOVA design
#' @param Xlist a \code{list} of X matrices, one per subject.
#' @param design a \code{data.frame} providing the variables provided in \code{formula} argument.
#' @param scale whether to scale the variables
#' @export
#' 
#' @examples 
#' 
#' f1 <- factor(rep(letters[1:4], 25))
#' f2 <- factor(rep(1:4, each=25))
#' 
#' des <- data.frame(f1=f1,f2=f2)
#' 
#' Xlist <- replicate(5, matrix(rnorm(100*100), 100, 100), simplify=FALSE)
#' 
#' res <- muasca(~ f1*f2, Xlist, design=des, ncomp=3)
muasca <- function(formula, Xlist, ncomp=2, design, scale=FALSE, A=NULL) {
  assert_that(inherits(Xlist, "list"))
  assertthat::assert_that(all(sapply(Xlist, function(x) is.matrix(x))))
  
  block_indices <- block_indices(Xlist)
  
  tform <- terms(formula)
  facs <- attr(tform, "factors")
  
  ## get the order number of each term (1 = main effect, 2 = two-way interaction, etc.)
  termorder <- apply(facs,2, function(x) sum(x > 0))
  terms <- names(termorder)
  
  ## the set of orders
  orders <- sort(unique(termorder))
  
  if (length(ncomp) == 1) {
    ncomp <- rep(ncomp, length(terms))
  }
  
  assert_that(length(ncomp) == length(terms))
  
  if (is.data.frame(design) && nrow(design) == nrow(Xlist[[1]])) {
    designL <- replicate(length(Xlist), design, simplify=FALSE)
  } else if (is.list(design) && length(design) == length(Xlist)) {
    designL <- design
  } else {
    stop("musu_asca: design must be a data.frame or a list of data.frames whose rows match the corresponding entries of Xlist")
  }
  
  assert_that(length(designL) == length(Xlist))
  assert_that(all(sapply(1:length(Xlist), function(i) nrow(Xlist[[i]]) == nrow(designL[[i]]))))
  

  ## construct a list of conditions for each factor combination
  Yfacl <- lapply(1:ncol(facs), function(i) {
    find <- which(facs[,i] > 0)
    facnames <- row.names(facs)[find]
    lapply(designL, function(d) {
      if (length(facnames) > 1) {
        do.call(function(...) interaction(..., sep=":"), d[, facnames])
      } else {
        d[, facnames]
      }
    })
  })
  
  names(Yfacl) <- terms
  
  ## the highrest order factor combination
  Ymaximal <- lapply(designL, function(d) {
    if (length(row.names(facs)) > 1) {
      do.call(function(...) interaction(..., sep=":"), d[, row.names(facs)])
    } else {
      d[, row.names(facs)]
    }
  })
  
  
  lens <- sapply(Ymaximal, function(x) length(levels(x)))
  
  assert_that(all(lens[1] == lens))
  
  main_terms <- terms[termorder == 1]
  all_terms <- row.names(facs)
  
  dgrid <- construct_design_table(termorder, facs, Ymaximal, all_terms, designL)

  Yfacl[[names(dgrid)[ncol(dgrid)]]] <- Ymaximal
  
  get_lower_form <- function(ord) {
    tnames <- colnames(facs)[termorder < ord]
    meat <-paste0("(", paste(tnames, collapse= " + "),")")
    as.formula(paste(" ~ ", meat))
  }
  
  expand_scores <- function(scores, y) {
    mind <- match(as.character(y), row.names(scores))
    scores[mind,]
  }
  
  res <- lapply(1:length(termorder), function(i) {
    ord <- termorder[i]
    message("pca decomposition for factor ", colnames(facs)[i])
    
    if (ord == 1) {
      form <- as.formula(paste("~ ", colnames(facs)[i], "-1"))
      
      Gl <- lapply(designL, function(des) {
        model.matrix(form, data=des)
      })
      
      Yl <- Yfacl[[i]]
      
      bres <- mubada(Yl, Xlist, ncomp=min(ncomp[i], length(levels(Yl[[1]]))), 
                     center=TRUE, scale=scale, normalization="None", A=A)
      

      list(G=Gl, Y=Yl, form=form, lower_form = ~ 1, term=colnames(facs)[i], 
           ex_scores=expand_scores(scores(bres), dgrid[,terms[i]]), 
           scores=scores(bres), bada_result=bres)
    
    } else {
      lower_form <- get_lower_form(ord)
      
      XresidL <- lapply(1:length(Xlist), function(i) {
        residualize(lower_form, Xlist[[i]], designL[[i]])
      })
      
      form <- as.formula(paste("~ ", colnames(facs)[i], "-1"))
      
      Gl <- lapply(designL, function(des) {
        G <- model.matrix(form, data=des)
      })
      
      Yl <- Yfacl[[i]]
      
      bres <- mubada(Yl, XresidL, ncomp=min(ncomp[i], length(levels(Yl[[1]]))), center=TRUE, scale=scale, 
                     normalization="None")
      
      list(G=Gl, Y=Yl, form=form, lower_form = lower_form, term=colnames(facs)[i], 
           ex_scores=expand_scores(bres$scores, dgrid[,terms[i]]), scores=scores(bres), bada_result=bres)
    }
  })
  
  #permute <- function(obj) {
  #  idx <- sample(1:nrow(obj$X))
  #  list(X=obj$X, Y=obj$Y[idx,], idx=idx, group=group[idx])
  #}
  
  expanded_scores <- lapply(res, "[[", "ex_scores")
  scores <- lapply(res, "[[", "scores")
  
  names(scores) <- colnames(facs)
  names(expanded_scores) <- colnames(facs)
  names(res) <- colnames(facs)

  sc <- do.call(cbind, expanded_scores)
  row.names(sc) <- levels(Ymaximal[[1]]) 
    
  ret <- list(
    ntables=length(Xlist),
    ncomp=ncomp,
    Xlist=Xlist,
    Yl=Yfacl,
    Ymaximal=Ymaximal,
    fac_design=dgrid,
    results=res,
    reduced_scores=scores,
    scores=sc,
    loadings=do.call(cbind, lapply(res, function(x) loadings(x$bada_result))),
    formula=formula,
    design=designL,
    refit=refit,
    terms=colnames(facs)
  )
  
  class(ret) <- c("muasca", "projector")
  ret
  
}

#' @export
loadings.musu_asca <- function(x) {
  do.call(cbind, lapply(res, function(x) loadings(x$bada_result)))
}

#' @export
scores.musu_asca <- function(x) {
  x$scores
}

#' @export
bootstrap.musu_asca <- function(x, niter=100, term=res$terms[1], ncomp=x$ncomp[1]) {
  bootstrap(x$results[[term]]$bada_result, niter=niter, ncomp=ncomp)
}


#' @export
predict.musu_asca <- function(x, newdata, type=c("class", "prob", "scores", "crossprod", "distance", "cosine"), 
                              block_index=1:x$ntables) {
  type <- match.arg(type)
  assert_that(is.matrix(newdata))
  assert_that(length(block_index) == 1 || length(block_index) == x$ntables)
  
  if (length(block_index) == x$ntables) {
    assert_that(ncol(newdata) == sum(sapply(x$Xlist, ncol)))
  }
  
  terms <- x$terms
  
  fscores <- do.call(cbind, lapply(terms, function(tname) {
      xres <- x$results[[tname]]
      predict(xres$bada_result, newdata, type="scores",block_index=block_index)
  }))
  
  preds <- scorepred(fscores, x$scores, type=type, class_name=FALSE)
  
  if (type == "class") {
    ret <- x$fac_design[preds,]
    row.names(ret) <- NULL
    ret
  } else {
    preds
  }
}


asca_subset <- function(x, fidx) {
  .Xlist <- lapply(1:length(x$Xlist), function(i) {
    xi <- x$Xlist[[i]]
    xi[fidx[[i]],,drop=FALSE]
  })
  
  .d <- lapply(1:length(x$Xlist), function(i) {
    des <- x$design[[i]]
    des[fidx[[i]],]
  })
  
  list(x=.Xlist, design=.d)
  
}


#' @export
performance.musu_asca <- function(x, ncomp=x$ncomp, blocks, term=names(x$fac_design)[ncol(x$fac_design)], metric=c("ACC", "AUC")) {
  assertthat::assert_that(term %in% names(x$fac_design))
  
  metric <- match.arg(metric)
  
  folds <- lapply(blocks, function(bind) split(1:length(bind), bind))
  
  yobs <- x$Yl[[term]]

  res <- lapply(seq_along(folds[[1]]), function(fnum) {
    message("musu_asca: performance fold: ", fnum)
    
    fidx <- lapply(folds, "[[", fnum)
    xsub <- asca_subset(x, lapply(fidx, "*", -1))
    
    rfit <- refit(x, xsub$x, .design=xsub$design,.ncomp=ncomp) 
    
    xsubout <- asca_subset(x, fidx)
    
    preds <- lapply(1:rfit$ntables, function(i) {
      if (metric == "ACC") {
        predict.musu_asca(rfit, xsubout$x[[i]], type="class", block_index=i)
      } else {
        predict.musu_asca(rfit, xsubout$x[[i]], type="prob", block_index=i)
      }
    })
    
  })
  
  if (metric == "ACC") {
    perf <- lapply(1:x$ntables, function(tind) {
      ptabs <- lapply(res, "[[", tind)
      p <- unlist(lapply(ptabs, function(x) x[, term]))
      
      ord <- unlist(folds[[tind]])
      data.frame(block_index=tind, pred=p, observed=yobs[[tind]][ord])
    })
    
    acc <- sapply(perf, function(x) sum(as.character(x$pred)==as.character(x$observed))/nrow(x))
    total_acc <- unlist(lapply(perf, "[[","pred")) == unlist(lapply(perf, "[[","observed"))
    total_acc <- sum(total_acc)/length(total_acc)
    list(acc=acc, total_acc=total_acc, pred_table=do.call(rbind, perf))
  } else {

    perf <- lapply(1:x$ntables, function(tind) {
      ptabs <- lapply(res, "[[", tind)
      ptab <- do.call(rbind, ptabs)
      rowlabs <- as.character(x$fac_design[,term])
      colnames(ptab) <- rowlabs
      
      ord <- unlist(folds[[tind]])
      auc <- .combinedAUC(ptab, yobs[[tind]][ord], rowlabs)
      list(block_index=tind, pred=ptab, observed=yobs[[tind]][ord])
    })
   
    auc <- sapply(perf, function(x) .combinedAUC(x$pred, x$observed, colnames(x$pred)))
    ptab_all <- do.call(rbind, lapply(perf, function(x) x$pred))
    yobs_all <- unlist(lapply(perf, "[[", "observed"))
    auc_all <- .combinedAUC(ptab_all, yobs_all, colnames(ptab_all))
    
    list(auc=auc, total_auc=auc_all, pred_list=perf)
  }
  
}



.combinedAUC <- function(Pred, Obs, predlabels) {

  mean(sapply(levels(Obs), function(lev) {
    pos <- Obs == lev
    pclass <- rowMeans(Pred[,predlabels==lev, drop=FALSE])
    pother <- rowMeans(Pred[,predlabels!=lev, drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)
  }))
}

