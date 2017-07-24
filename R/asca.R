
gen_threeway_mat <- function(n1,n2,n3, nvox, index) {
  nk <- n1*n2*n3

  grid <- expand.grid(A=paste0("A", 1:n1), B=paste0("B", 1:n2), C=paste0("C", 1:n3))
  grid$sid <- index
  mat <- matrix(rnorm(nk*nvox), nk, nvox)
  list(mat=mat, design=grid)

}

#Xs <- lapply(1:10, function(i) gen_threeway_mat(3,3,5, 10, i))
#Xlist <- lapply(Xs, "[[", "mat")
#design <- Xs[[1]]$design
#formula = ~ A + B + C + A:B + A:C + B:C + A:B:C
#strata <- "sid"


#' collapse
#' @param form the formula defining the model to collapse over
#' @param X the matrix to collpased ober rows
#' @param design the \code{data.frame} containing the variables reference in \code{form}
#' 
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


#' residualize
#' @param form the formula defining the model to fit for residuals
#' @param X the response matrix
#' @param design the \code{data.frame} containing the design variables specified in \code{form} argument.
#' @export
residualize <- function(form, X, design) {
  options(contrasts = c("contr.sum", "contr.poly"))
  modmat <- model.matrix(form, data=design)
  resid(lsfit(modmat, X, intercept=FALSE))
}
  
  
#' musu_asca
#' 
#' 
#' @param Xlist the list of data matrices
#' @param formula a formula specifying the ANOVA design
#' @param design a \code{data.frame} providing the variables provided in \code{formula} argument.
#' @param center whether to center the variables
#' @param scale whether to scale the variables
#' @param svd.method the svd method to use.
#' @export
musu_asca <- function(Xlist, formula, ncomp=2, design, center=TRUE, scale=FALSE, svd.method="fast") {
  Xlist <- lapply(Xlist, as.matrix)
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
    designL <- replicate(design, length(Xlist))
  } else if (is.list(design) && length(design) == length(Xlist)) {
    designL <- design
  } else {
    stop("musu_asca: design must be a data.frame or a list of data.frames whose rows match the correspinding entries of Xlist")
  }
  
  assert_that(length(designL) == length(Xlist))
  assert_that(all(sapply(1:length(Xlist), function(i) nrow(Xlist[[i]]) == nrow(designL[[i]]))))
  
  blockInd <- blockIndices(Xlist)
  
  
  
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
  
  dgrid <- expand.grid(lapply(all_terms, function(tname) levels(designL[[1]][[tname]])))
  
  names(dgrid) <- all_terms
  

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
  
  Yfacl[[names(dgrid)[ncol(dgrid)]]] <- Ymaximal
  
  refit <- function(.Xlist, .design, .ncomp=ncomp) { 
    musu_asca(.Xlist, formula, ncomp=.ncomp, design=.design, center=center, scale=scale, svd.method=svd.method) 
  }
  

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
      
      bres <- musu_bada(Yl, Xlist, ncomp=min(ncomp[i], length(levels(Yl[[1]]))), center=TRUE, svd.method=svd.method, normalization="None")
      
      list(G=Gl, Y=Yl, form=form, lower_form = ~ 1, term=colnames(facs)[i], ex_scores=expand_scores(bres$scores, dgrid[,terms[i]]), scores=bres$scores, bada_result=bres)
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
      
      bres <- musu_bada(Yl, XresidL, ncomp=min(ncomp[i], length(levels(Yl[[1]]))), center=TRUE, svd.method=svd.method, normalization="None")
      list(G=Gl, Y=Yl, form=form, lower_form = lower_form, term=colnames(facs)[i], ex_scores=expand_scores(bres$scores, dgrid[,terms[i]]), scores=bres$scores, bada_result=bres)
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
    loadings=do.call(cbind, lapply(res, function(x) loadings(x$bada_result$pca_fit))),
    formula=formula,
    design=designL,
    refit=refit,
    terms=colnames(facs)
  )
  
  
  class(ret) <- c("projector", "musu_asca")
  ret
  
}

#' @export
loadings.musu_asca <- function(x) {
  do.call(cbind, lapply(res, function(x) loadings(x$bada_result$pca_fit)))
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
                              table_index=1:x$ntables) {
  type <- match.arg(type)
  assert_that(is.matrix(newdata))
  assert_that(length(table_index) == 1 || length(table_index) == x$ntables)
  
  if (length(table_index) == x$ntables) {
    assert_that(ncol(newdata) == sum(sapply(x$Xlist, ncol)))
  }
  
  terms <- x$terms
  
  fscores <- do.call(cbind, lapply(terms, function(tname) {
      xres <- x$results[[tname]]
      predict(xres$bada_result, newdata, type="scores",table_index=table_index)
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
    
    rfit <- x$refit(xsub$x, .design=xsub$design,.ncomp=ncomp) 
    
    xsubout <- asca_subset(x, fidx)
    
    preds <- lapply(1:rfit$ntables, function(i) {
      if (metric == "ACC") {
        predict.musu_asca(rfit, xsubout$x[[i]], type="class", table_index=i)
      } else {
        predict.musu_asca(rfit, xsubout$x[[i]], type="prob", table_index=i)
      }
    })
    
  })
  
  if (metric == "ACC") {
    perf <- lapply(1:x$ntables, function(tind) {
      ptabs <- lapply(res, "[[", tind)
      p <- unlist(lapply(ptabs, function(x) x[, term]))
      
      ord <- unlist(folds[[tind]])
      data.frame(table_index=tind, pred=p, observed=yobs[[tind]][ord])
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
      list(table_index=tind, pred=ptab, observed=yobs[[tind]][ord])
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

