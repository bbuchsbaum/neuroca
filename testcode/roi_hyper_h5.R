
path <- "/Users/bbuchsbaum/analysis/hyper/all_surface_h5/"

use_s5 <- TRUE
impute_nback <- TRUE
use_video <- TRUE
roi_thresh <- .1

library(hdf5r)
library(softImpute)
library(lme4)
library(softImpute)
library(matrixStats)
library(missMDA)

sids <- scan(paste0(path, "sids"), "")[-1]

h5files <- paste0(path, "/", sids, "_surface_data.h5")
hfset <- lapply(h5files, function(fname) H5File$new(fname, mode="r+"))
names(hfset) <- sids

zerostr <- function(vals) {
  ndigits <- log(max(vals), 10) + 1
  nzeros <- ndigits - as.integer(log(vals,10)) -1
  prefix <- sapply(nzeros, function(nz) paste(rep("0", times=nz), collapse=""))
  paste(prefix, vals, sep="")  
}


load_vdat <- function(roinum) {
  lapply(sids, function(sid) {
    dat <- hfset[[paste0(as.character(sid))]]
    
    dat <- lapply(roinum, function(rnum) { dat[[as.character(rnum)]] })
    
    prefix <- if (use_s5) "smoothed_" else ""
    
    v1 <- lapply(1:length(roinum), function(i) dat[[i]][[paste0(prefix, "video_1")]])
    v2 <- lapply(1:length(roinum), function(i) dat[[i]][[paste0(prefix, "video_2")]])
    v3 <- lapply(1:length(roinum), function(i) dat[[i]][[paste0(prefix, "video_3")]])
   
    des <- data.frame(labels=c(v1[[1]][["labels"]][], v2[[1]][["labels"]][], v3[[1]][["labels"]][]),
                      run=factor(rep(1:3, each=364)))
    
    matlist <- lapply(1:length(roinum), function(i) {
      rbind(t(v1[[i]][["data"]][,]), t(v2[[i]][["data"]][,]), t(v3[[i]][["data"]][,]))
    })
    
    mat <- do.call(cbind, matlist)
    
  
    keep <- apply(mat,2, function(x) sum(x==0)) == 0
    mat <- mat[,keep]
  
    nodes <- unlist(lapply(v1, function(v) v[["nodes"]][]))
    nodes <- nodes[keep]
    
    rmat <- residualize(~ run, mat, des)
    rmat <- t(scale(t(rmat)))
  
    list(sid=sid, mat=rmat, design=des, keep=keep, roi=roinum, nodes=nodes)
  })
}

load_ndat <- function(roinum) {
  lapply(1:length(sids), function(i) {
    sid <- sids[i]
    dat <- hfset[[paste0(as.character(sid))]]
    dat <- lapply(roinum, function(rnum) { dat[[as.character(rnum)]] })

    prefix <- if (use_s5) "smoothed_" else ""
    v1 <- lapply(1:length(roinum), function(i) dat[[i]][[paste0(prefix, "nback")]])
    
    des <- read.table(paste0(path, sid, "_nback_design.txt"))
    des$label_version <- factor(paste0(as.character(des$Label), "_", des$Version))
  
    matlist <- lapply(1:length(roinum), function(i) {
      t(v1[[i]][["data"]][,])
    })
    
    mat <- do.call(cbind, matlist)
  
    keep <- apply(mat,2, function(x) sum(x==0)) == 0
    mat <- mat[,keep]
    
    nodes <- unlist(lapply(v1, function(v) v[["nodes"]][]))
    nodes <- nodes[keep]
    
 
    des$run <- factor(des$run)
    rmat <- residualize(~ run, mat, des)
    rmat <- t(scale(t(rmat)))
    
    mat <- group_means(des$label_version, rmat)
    mat <- t(scale(t(mat), center=FALSE))
    
    mat_miss <- matrix(NA, nrow(mat)*2, ncol(mat))
  
    levs <- levels(des$label_version)
    labs <- levels(des$Label)
    
    ver <- sapply(strsplit(levs, "_"), function(x) {
      x[length(x)]
    })
  
    for (i in 1:length(ver)) {
      if (ver[i] == "1") {
        mat_miss[i*2-1,] <- mat[i,]
      } else {
        mat_miss[i*2,] <- mat[i,]
      }
    }
    
    des2 <- data.frame(label=paste0(rep(labs, each=2), "_", rep(1:2, 90)), version=rep(1:2, 90))
    list(sid=sid, mat=mat_miss, design=des2, keep=keep, len=ncol(mat), roi=roinum, nodes=nodes)
  })
}

load_idat <- function(roinum) {
  lapply(1:length(sids), function(i) {
    sid <- sids[i]
    print(sid)
    dat <- hfset[[paste0(as.character(sid))]]
    
    dat <- lapply(roinum, function(rnum) { dat[[as.character(rnum)]] })

    
    prefix <- if (use_s5) "smoothed_" else ""
    v1 <- lapply(1:length(roinum), function(i) dat[[i]][[paste0(prefix, "imagery")]])
    
    des <- read.table(paste0(path, sid, "_imagery_design.txt"))
     
    label_version <- factor(paste0(as.character(des$label), "_", des$version))
    
    matlist <- lapply(1:length(roinum), function(i) {
      t(v1[[i]][["data"]][,])
    })
    
    mat <- do.call(cbind, matlist)
    keep <-  keep <- apply(mat,2, function(x) sum(x==0)) == 0
    mat <- mat[, keep]
    
    nodes <- unlist(lapply(v1, function(v) v[["nodes"]][]))
    nodes <- nodes[keep]
    
    des$run <- factor(des$run)
    rmat <- residualize(~ run, mat, des)
    rmat <- t(scale(t(rmat)))
    
    spmat <- split_matrix(rmat, ordered(des$time))
    
    spmat <- lapply(spmat, function(m) {
      sdavg <- mean(apply(m,1,sd))
      m <- m/sdavg
      m <- scale(m, scale=FALSE)
      
      
    })

    list(sid=sid, spmat=spmat, design=des, keep=keep, len=ncol(mat), roi=roinum, nodes=nodes)
    
  })
}

roinum <- c(11102, 11119,11120, 11121,11122,11123,11159, 11160, 11161,11162,11145)
roinum <- c(roinum, roinum+1000)

vdat <- load_vdat(roinum)
ndat <- load_ndat(roinum)
idat <- load_idat(roinum)

bothmat <- do.call(cbind, lapply(1:length(vdat), function(i) rbind(vdat[[i]]$mat, ndat[[i]]$mat)))

group <- sapply(vdat, function(x) length(x$nodes))

#both_imp <- imputeMFA(as.data.frame(bothmat), group, ncp=12, method="Regularized", maxiter=100, threshold = 1e-05)
#both_mat <- to_block_matrix(as.matrix(both_imp$completeObs), group)
bothmat_blocked <- to_block_matrix(bothmat, group)
both_imp <- impute_mfa(bothmat_blocked, ncomp=25, normalization="MFA")
both_mat <- to_block_matrix(both_imp, group)


both_xlist <- as.list(both_mat)
Y <- factor(c(zerostr(1:length(vdat[[1]]$design$labels)), as.character(ndat[[1]]$design$label)))

mres1 <- mubada(Y, both_xlist, ncomp=12, scale=TRUE, normalization="MFA")
#mres2 <-  mubada(Y, both_xlist, ncomp=12, scale=TRUE, normalization="MFA")
#mres3 <-  mubada(Y, both_xlist, ncomp=12, scale=TRUE, normalization="RV")

get_predictions <- function(fit, Xim, nback_ind, nc,i) {
  
  if (class(fit)[1] == "sca") {
    pred <- predict(fit, Xim, ncomp=nc, table_index=i)
    cpred <- scorepred(pred, scores(fit), type="cosine")
    cpred[,nback_ind]
  } else {
    
    cpred <- predict(fit, Xim, type="cosine", ncomp=nc, table_index=i)
    cpred <- cpred[,nback_ind]
  }
}


do_projection <- function(mfit, idatlist, nc=6) {
  #nback_ind <- 274:453
  nback_ind <- 1093:(1093+179)
  
  tvals <- sort(unique(idatlist[[1]]$design$time))
  
  do.call(rbind, lapply(1:length(tvals), function(tind) {
    Xim <- lapply(idatlist, function(x) x$spmat[[tind]])
    
    pres_cos <- lapply(1:length(Xim), function(i) {
      print(i)
      cpred <- get_predictions(mfit, Xim[[i]], nback_ind, nc, i)
      
      des <- subset(idatlist[[i]]$design, time==tvals[tind])
      lver <- paste0(des$label, "_", des$version)
      ind <- match(as.character(lver), as.character(Y[nback_ind]))
      
      react <- cpred[cbind(1:length(lver), ind)]
      
      other_version <- ifelse(des$version == 1, 2, 1)
      other_label <- paste0(des$label, "_", other_version)
      ind_other <- match(as.character(other_label), Y[nback_ind])
      
      react_other <- cpred[cbind(1:length(lver), ind_other)]
      
      
      dfx1 <- do.call(rbind, replicate(3, data.frame(time = tvals[tind],
                                                     cond = ifelse(des$cresp == 1, "same", "different"),
                                                     acc = des$acc,
                                                     zvivid = scale(des$vivid),
                                                     zconf = scale(des$confidence),
                                                     sid=sids[i],
                                                     label=des$label,
                                                     label_version=lver), simplify=FALSE))
      
      dfx1$type <- rep(c("react", "react_other", "react_diff"), each=nrow(dfx1)/3)
      dfx1$react <- c(react, react_other, react-react_other)
      dfx1
    })
    
    do.call(rbind, pres_cos)
  }))
}

#allfits <- list(mres1, mres3, mres4, sres2, sres3, sres4)
#names(allfits) <- c("mfa", "mfa-scaled", "pmfa-scaled", "sca-ind", "sca-pf2", "sca-ecp")

allfits <- list(mres1, mres2, mres3)
names(allfits) <- c("mfa", "mfa_scaled", "mfa_rv")

pres_cos <- lapply(1:length(allfits), function(i) {
  pres_cos = do_projection(allfits[[i]], idat, nc=6)
  pres_cos$model <- names(allfits)[i]
  pres_cos
})

pres_cos <- do.call(rbind, pres_cos)

pres_cos_comp <- do.call(rbind, lapply(3:12, function(i) {
  pres_cos <- do_projection(allfits[[3]], idat, nc=i)
  pres_cos$model <- names(allfits)[1]
  pres_cos$ncomp <- i
  pres_cos
}))

library(dplyr)
library(ggplot2)

reac_by_time_model <- pres_cos %>% group_by(time, type, model) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_acc_model <- pres_cos %>% group_by(time, acc, type,model) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_cond_model <- pres_cos %>% group_by(time, cond, type,model) %>% summarize(react=mean(react, na.rm=TRUE))

qplot(time, react, colour=model, data=subset(reac_by_time_model), colour=model, geom=c("point", "line")) + facet_wrap( ~ type)
qplot(time, react, colour=factor(acc), data=subset(reac_by_acc_model,type=="react"), geom=c("point", "line")) + facet_wrap( ~ model)
qplot(time, react, colour=factor(cond), data=subset(reac_by_cond_model,type=="react"), geom=c("point", "line")) + facet_wrap( ~ model)


reac_by_time_comp <- pres_cos_comp %>% group_by(time, type, ncomp) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_acc_comp <- pres_cos_comp %>% group_by(time, acc, type,ncomp) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_viv_comp <- pres_cos_comp %>% filter(acc==1) %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 4)) %>% ungroup() %>% group_by(time, qvivid, type,ncomp) %>% 
  summarize(react=mean(react, na.rm=TRUE))

reac_by_cond_comp <- pres_cos_comp %>%  group_by(time, cond, type,ncomp) %>% summarize(react=mean(react, na.rm=TRUE)) 
reac_by_cond_acc_comp <- pres_cos_comp %>%  group_by(time, cond, acc, type,ncomp) %>% summarize(react=mean(react, na.rm=TRUE)) 


reac_by_conf_comp <- pres_cos_comp %>%  group_by(sid) %>% mutate(qconf=ntile(zconf, 4)) %>% ungroup() %>% group_by(time, qconf, type,ncomp) %>% 
  summarize(react=mean(react, na.rm=TRUE))



qplot(time, react, colour=factor(acc), data=subset(reac_by_acc_comp,type=="react"), geom=c("point", "line")) + facet_wrap( ~ ncomp)
qplot(time, react, colour=factor(ncomp), data=subset(reac_by_time_comp), colour=model, geom=c("point", "line")) + facet_wrap( ~ type)
qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv_comp,type=="react" & !is.na(qvivid)), geom=c("point", "line")) + facet_wrap( ~ ncomp)
qplot(time, react, colour=factor(qconf), data=subset(reac_by_conf_comp,type=="react" & !is.na(qconf)), geom=c("point", "line")) + facet_wrap( ~ ncomp)


qplot(time, react, colour=factor(cond), data=subset(reac_by_cond_comp, type=="react"), geom=c("point", "line")) + facet_wrap(~ ncomp)
qplot(time, react, colour=factor(acc), data=subset(reac_by_cond_acc_comp, type=="react_diff" & ncomp==6), geom=c("point", "line")) + 
  facet_wrap(~ cond)



