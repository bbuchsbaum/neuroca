library(missMDA)


zerostr <- function(vals) {
  ndigits <- log(max(vals), 10) + 1
  nzeros <- ndigits - as.integer(log(vals,10)) -1
  prefix <- sapply(nzeros, function(nz) paste(rep("0", times=nz), collapse=""))
  paste(prefix, vals, sep="")  
}


alldat <- readRDS("/Users/bbuchsbaum/analysis/hyper/ventral_surface/ROI_11121_alldat.RDS")
sids <- sapply(alldat$vdat, function(x) x$sid)
vmatlist <- lapply(alldat$vdat, function(x) x$mat)
ndatlist <- lapply(alldat$ndat, function(x) x$mat)
idatlist <- lapply(alldat$idat, function(x) x)

bothmat <- do.call(cbind, lapply(1:length(vmatlist), function(i) rbind(vmatlist[[i]], ndatlist[[i]])))
group <- sapply(vmatlist, ncol)

#bothmat_blocked <- to_block_matrix(bothmat, group)
#both_imp <- impute_mfa(bothmat_blocked, ncomp=25, normalization="MFA")
both_imp <- imputeMFA(as.data.frame(bothmat), group, ncp=12, method="Regularized", maxiter=100, threshold = 1e-05)
both_mat <- to_block_matrix(as.matrix(both_imp$completeObs), group)
#both_mat <- to_block_matrix(both_imp, group)


both_xlist <- as.list(both_mat)
both_xlist <- lapply(both_xlist, function(x) t(scale(t(x))))

Y <- factor(c(zerostr(1:nrow(vmatlist[[1]])), as.character(alldat$ndat[[1]]$design$label)))

mres1 <- mubada(Y, both_xlist, ncomp=12, normalization="MFA")
mres2 <- procrusteanize.mubada(mres1, ncomp=12)
mres3 <-  mubada(Y, both_xlist, ncomp=12, scale=TRUE, normalization="MFA")
mres4 <- procrusteanize.mubada(mres3, ncomp=12)
mres5 <- mubada(Y, both_xlist, ncomp=12, scale=TRUE, normalization="None")
mres6 <- mubada(Y, both_xlist, ncomp=12, scale=TRUE, normalization="RV")

sres1 <- sca(block_matrix(both_xlist), scale=TRUE, ncomp=12, type="sca-p")
sres2 <- sca(block_matrix(both_xlist), scale=TRUE, ncomp=12, type="sca-ind")
sres3 <- sca(block_matrix(both_xlist), scale=TRUE, ncomp=12, type="sca-pf2")
sres4 <- sca(block_matrix(both_xlist), scale=TRUE, ncomp=12, type="sca-ecp")



get_predictions <- function(fit, Xim, nback_ind, nc,i) {

  if (class(fit)[1] == "sca") {
    pred <- predict(fit, Xim, ncomp=nc, table_index=i)[[1]]
    cpred <- scorepred(pred, scores(fit), type="cosine")
    cpred[,nback_ind]
  } else {
    
    cpred <- predict(fit, Xim, type="cosine", ncomp=nc, table_index=i)
    cpred <- cpred[,nback_ind]
  }
}


do_projection <- function(mfit, idatlist, nc=6) {
  nback_ind <- 274:453
  
  tvals <- sort(unique(idatlist[[1]]$design$time))

  do.call(rbind, lapply(1:length(tvals), function(tind) {
    Xim <- lapply(idatlist, function(x) x$spmat[[tind]])
  
    pres_cos <- lapply(1:length(Xim), function(i) {
      print(i)
      cpred <- get_predictions(mfit, Xim[[i]], nback_ind, nc, i)
      
      des <- subset(idatlist[[i]]$design, time==tvals[tind])
      lver <- des$LabelVersion
      ind <- match(as.character(lver), Y[nback_ind])
  
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
                       label_version=des$LabelVersion), simplify=FALSE))
    
      dfx1$type <- rep(c("react", "react_other", "react_diff"), each=nrow(dfx1)/3)
      dfx1$react <- c(react, react_other, react-react_other)
      dfx1
    })
  
    do.call(rbind, pres_cos)
  }))
}

allfits <- list(mres1, mres3, mres4, sres2, sres3, sres4)
names(allfits) <- c("mfa", "mfa-scaled", "pmfa-scaled", "sca-ind", "sca-pf2", "sca-ecp")

pres_cos <- lapply(1:length(allfits), function(i) {
  pres_cos = do_projection(allfits[[i]], idatlist, nc=8)
  pres_cos$model <- names(allfits)[i]
  pres_cos
})

pres_cos <- do.call(rbind, pres_cos)

pres_cos_comp <- do.call(rbind, lapply(2:12, function(i) {
  pres_cos = do_projection(allfits[[4]], idatlist, nc=i)
  pres_cos$model <- names(allfits)[4]
  pres_cos$ncomp <- i
  pres_cos
}))

library(dplyr)
library(ggplot2)

reac_by_time_model <- pres_cos %>% group_by(time, type, model) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_acc_model <- pres_cos %>% group_by(time, acc, type,model) %>% summarize(react=mean(react, na.rm=TRUE))

qplot(time, react, colour=model, data=subset(reac_by_time_model), colour=model, geom=c("point", "line")) + facet_wrap( ~ type)
qplot(time, react, colour=factor(acc), data=subset(reac_by_acc_model,type=="react"), geom=c("point", "line")) + facet_wrap( ~ model)

reac_by_time_comp <- pres_cos_comp %>% group_by(time, type, ncomp) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_acc_comp <- pres_cos_comp %>% group_by(time, acc, type,ncomp) %>% summarize(react=mean(react, na.rm=TRUE))
reac_by_viv_comp <- pres_cos_comp %>% filter(acc==1) %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 4)) %>% ungroup() %>% group_by(time, qvivid, type,ncomp) %>% 
  summarize(react=mean(react, na.rm=TRUE))

reac_by_cond_comp <- pres_cos_comp %>%  group_by(time, cond, type,ncomp) %>% summarize(react=mean(react, na.rm=TRUE)) 
reac_by_conf_comp <- pres_cos_comp %>%  group_by(sid) %>% mutate(qconf=ntile(zconf, 4)) %>% ungroup() %>% group_by(time, qconf, type,ncomp) %>% 
  summarize(react=mean(react, na.rm=TRUE))



qplot(time, react, colour=factor(ncomp), data=subset(reac_by_time_comp), colour=model, geom=c("point", "line")) + facet_wrap( ~ type)
qplot(time, react, colour=factor(acc), data=subset(reac_by_acc_comp,type=="react"), geom=c("point", "line")) + facet_wrap( ~ ncomp)
qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv_comp,type=="react" & !is.na(qvivid)), geom=c("point", "line")) + facet_wrap( ~ ncomp)
qplot(time, react, colour=factor(qconf), data=subset(reac_by_conf_comp,type=="react" & !is.na(qconf)), geom=c("point", "line")) + facet_wrap( ~ ncomp)


qplot(time, react, colour=factor(cond), data=subset(reac_by_cond_comp, type=="react_diff"), geom=c("point", "line")) + facet_wrap(~ ncomp)





# 
reac_by_time <- pres_cos %>% group_by(time, type) %>% summarize(react=mean(react, na.rm=TRUE))
                                                       
qplot(time, react, data=reac_by_time, geom=c("point", "line")) + facet_wrap( ~ type)
# 
reac_by_acc <- pres_cos %>% group_by(time, acc, type) %>% summarize(react=mean(react, na.rm=TRUE))
qplot(time, react, colour=factor(acc), data=reac_by_acc, geom=c("point", "line")) + facet_wrap( ~ type)
# 
# 
reac_by_viv <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3)) %>% ungroup() %>% group_by(time, qvivid, type) %>% 
   summarize(react=mean(react, na.rm=TRUE))
#   
qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv, !is.na(qvivid)), geom=c("point", "line")) + facet_wrap( ~ type)

reac_by_viv_cond <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3)) %>% ungroup() %>% group_by(time, qvivid, cond, type) %>% 
summarize(react=mean(react, na.rm=TRUE))
 
qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv_cond, !is.na(qvivid)), geom=c("point", "line")) + facet_wrap(cond ~ type)
# 
# 
reac_by_viv_acc <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3)) %>% ungroup() %>% group_by(time, qvivid, acc, type) %>% 
   summarize(react=mean(react, na.rm=TRUE))
 
qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv_acc, !is.na(qvivid) & acc==1), geom=c("point", "line")) + 
  facet_wrap( ~ type)

qplot(time, react, colour=factor(acc), data=subset(reac_by_viv_acc, !is.na(qvivid) & qvivid > 2), geom=c("point", "line")) + 
  facet_wrap( ~ type)

reac_by_cond <- pres_cos %>%  group_by(time, cond, type) %>% summarize(react=mean(react, na.rm=TRUE)) 
qplot(time, react, colour=factor(cond), data=reac_by_cond, geom=c("point", "line")) + facet_wrap(~ type)


reac_by_cond_acc <- pres_cos %>%  group_by(time, cond,acc,type) %>% summarize(react=mean(react, na.rm=TRUE))
qplot(time, react, colour=factor(acc), facets = . ~ cond, data=reac_by_cond_acc, geom=c("point", "line")) + facet_wrap(cond ~ type)

reac_by_conf <- pres_cos %>%  group_by(sid) %>% mutate(qconf=ntile(zconf, 3)) %>% ungroup() %>% group_by(time, qconf, type) %>% 
summarize(react=mean(react, na.rm=TRUE))

qplot(time, react, colour=factor(qconf), data=subset(reac_by_conf, !is.na(qconf)), geom=c("point", "line")) + facet_wrap(~ type)

reac_by_conf_cond <- pres_cos %>% group_by(sid) %>% mutate(qconf=ntile(zconf, 3)) %>% ungroup()  %>% group_by(time, qconf, cond, type) %>% 
   summarize(react=mean(react, na.rm=TRUE))
 
qplot(time, react, colour=factor(qconf), facets = . ~ cond, data=subset(reac_by_conf_cond, !is.na(qconf)), geom=c("point", "line")) +
   facet_wrap(cond ~ type)
 
reac_by_conf_acc <- pres_cos %>%  group_by(sid) %>% mutate(qconf=ntile(zconf, 3)) %>% ungroup() %>% group_by(time, qconf, acc, type) %>% 
   summarize(react=mean(react, na.rm=TRUE))
 
qplot(time, react, colour=factor(qconf), data=subset(reac_by_conf_acc, !is.na(qconf)), geom=c("point", "line")) +
facet_wrap(acc ~ type)
 
reac_by_viv_conf <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3), qconf=ntile(zconf, 3)) %>% ungroup() %>% 
   group_by(time, qvivid, qconf, type) %>% 
   summarize(react=mean(react, na.rm=TRUE))

pres_cos$zmult <- (pres_cos$zvivid - min(pres_cos$zvivid,na.rm=TRUE)) * (pres_cos$zconf - min(pres_cos$zconf, na.rm=TRUE))
reac_by_qmult <- pres_cos %>% group_by(sid) %>% mutate(qmult=ntile(zmult, 3)) %>% ungroup() %>% 
  group_by(time, qmult, type,acc) %>% 
  summarize(react=mean(react, na.rm=TRUE))


 qplot(time, react, colour=factor(qmult), data=subset(reac_by_qmult, !is.na(qmult)), geom=c("point", "line")) + 
   facet_wrap(acc ~ type)

