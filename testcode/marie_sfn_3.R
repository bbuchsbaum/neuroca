library(tibble)
library(dplyr)
library(neuroim)
library(vegan)
library(ggrepel)
library(assertthat)
library(devtools)
library(energy)

devtools::load_all()



sids <- scan(paste0(PATH, "/sids.txt"), "")
PATH <- "~/Dropbox/Brad/neuropls/video_marie/SFN/"


mask <- loadVolume(paste0(PATH, "/clus_mask.nii"))
mask.idx <- which(mask>0)

loadMat <- function(sid) {
  message("loading", sid)
  dat <- readRDS(paste0(PATH, sid, "_block_betas.rds"))
}

dset <- lapply(sids, loadMat)
## normalize each image
XEBlocks <- lapply(dset, function(x) {
  idx <- which(x$design$Condition == "Encod")
  X <- x$X[idx,]
  t(apply(X,1, function(vals) { scale(vals) }))
})

XRBlocks <- lapply(dset, function(x) {
  idx <- which(x$design$Condition == "Recall")
  X <- x$X[idx,]
  t(apply(X,1, function(vals) { scale(vals) }))
})


hold_one_out <- function(Yl, Xl, ncomp=10, pfit=NULL) {
  res <- bind_rows(lapply(1:length(Yl), function(i) {
    print(i)
    Xtrain <- Xl[-i]
    Ytrain <- Yl[-i]
    
    Xtest <- Xl[[i]]
    Ytest <- Yl[[i]]
    
    if (is.null(pfit)) {
      print("refitting")
      pfit <- musubada(Ytrain, Xtrain, ncomp=10, center=TRUE, scale=FALSE, normalization="None")
    }
    
    creps <- code_replications(Ytest)
    p <- lapply(sort(unique(creps)), function(rnum) {
      ##browser()
      supX <- Xtest[creps != rnum,]
      supY <- Ytest[creps != rnum]
      
      yt <- Ytest[creps == rnum]
      
      
      predfun1 <- supplementary_predictor(pfit, supX, supY, ncomp=ncomp, type="cosine")
      predfun2 <- supplementary_predictor(pfit, supX, supY, ncomp=ncomp, type="class")
      
      prob <- predfun1(Xtest[creps == rnum,])
      
      auc <- combinedAUC(prob, yt) + .5
      pwinner <- prob[cbind(1:nrow(prob), as.integer(yt))]
      cls <- predfun2(Xtest[creps==rnum,])
      
      list(prob=pwinner, auc=rep(auc,length(yt)), class=cls)
    })
    
    data_frame(auc=unlist(lapply(p, "[[", "auc")), prob=unlist(lapply(p, "[[", "prob")), 
               repnum=creps, pred=unlist(lapply(p, "[[", "class")), 
               observed=Ytest, sid=sids[i])
    
  }))
  
  res$acc <- as.integer(res$observed == res$pred)
  res$repnum <- as.numeric(res$repnum)
  res
  #res2 <- res %>% group_by(sid, repnum) %>% dplyr::summarise(prob=mean(prob), acc=mean(acc))
  #res3 <- res %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))
  #res4 <- res %>% group_by(repnum) %>% dplyr::summarise(prob=mean(prob), acc=mean(acc))
}


Yvid_encode <- lapply(dset, function(df1) factor(df1$design$Video[df1$design$Condition == "Encod"]))
musu_encode_fit <- musubada(Yvid_encode, XEBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="None")

Yvid_recall <- lapply(dset, function(df1) factor(df1$design$Video[df1$design$Condition == "Recall"]))
musu_recall_fit <- musubada(Yvid_recall, XRBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="None")

encode_pred <- hold_one_out(Yvid_encode, XEBlocks, ncomp=10)
encode_pred_s1 <- encode_pred %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))
encode_pred_s2 <- encode_pred %>% group_by(sid, repnum) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))

recall_embed <- hold_one_out(Yvid_recall, XRBlocks, ncomp=10, pfit=musu_encode_fit)
recall_embed_s1 <- recall_embed %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))

## recall in recall
recall_pred <- hold_one_out(Yvid_recall, XRBlocks, ncomp=10)
recall_pred_s1 <- recall_pred %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))


## projection of recall on encoding fit
auc=unlist(lapply(1:length(Yvid_recall), function(i) {
  tmp <- scale(XRBlocks[[i]], center=TRUE, scale=FALSE)
  p1 <- predict(musu_encode_fit, tmp, type="cosine", ncomp=10, table_index=i)
  pwinner <- p1[cbind(1:nrow(p1), as.integer(Yvid_recall[[i]]))]
  mean(pwinner)
  #sum(p1 == Yvid[[i]])/length(Yvid[[i]])
  #combinedAUC(p1, Yvid_recall[[i]]) + .5
}))


df1 <- rbind(data.frame(cosine=encode_pred_s1$prob, classifier="encode", sid=sids),
             data.frame(cosine=recall_embed_s1$prob, classifier="embedded", sid=sids),
             data.frame(cosine=auc, classifier="cross-decode", sid=sids))
df1$group <- ifelse(df1$sid == 5001, "NC", "control")
df1$group[df1$sid == 5002] <- "HC"
df1$group = factor(df1$group)

pdf(file=paste0(PATH, "between_subject.pdf"))
theme_set(theme_gray(base_size = 18))
ggplot(df1, aes(x = factor(classifier), fill = factor(group), y = cosine)) +
  geom_dotplot(binaxis = "y", stackdir = "center") + labs(x= "Classifier Type", y="Accuracy", fill="Group")
dev.off()

write.table(df1, paste0(PATH, "/between_subject_results.txt"))

pcfit <- musu_encode_fit$pca_fit
enc_scores <- musu_encode_fit$scores


sf <- pcfit$d/sum(pcfit$d)
for (i in 1:length(XRBlocks)) {
  print(i)
  
  X <- XRBlocks[[i]]
  fscores <- do.call(rbind, replicate(nrow(X)/11, enc_scores, simplify=FALSE))
  cmat <- cor(X, fscores)
  cmat_scaled <- atanh(cmat) %*% diag(sf)
  cvals <- apply(cmat_scaled, 1, function(vals) mean(vals^2))
  # for (j in 1:5) {
  #   bv <- BrainVolume(cmat[,j], space(mask), indices=mask.idx)
  #   writeVolume(bv, paste0(PATH, "/", sids[i], "_", j, "_cor_loadings.nii"))
  # }
  
  bv <- BrainVolume(cvals, space(mask), indices=mask.idx)
  writeVolume(bv, paste0(PATH, "/", sids[i], "_recall_cor_scaled.nii"))
}


sf <- pcfit$d/sum(pcfit$d)
for (i in 1:length(XRBlocks)) {
  print(i)
  
  X <- XEBlocks[[i]]
  fscores <- do.call(rbind, replicate(nrow(X)/11, enc_scores, simplify=FALSE))
  cmat <- cor(X, fscores)
  cmat_scaled <- atanh(cmat) %*% diag(sf)
  cvals <- apply(cmat_scaled, 1, function(vals) mean(vals^2))
  # for (j in 1:5) {
  #   bv <- BrainVolume(cmat[,j], space(mask), indices=mask.idx)
  #   writeVolume(bv, paste0(PATH, "/", sids[i], "_", j, "_cor_loadings.nii"))
  # }
  
  bv <- BrainVolume(cvals, space(mask), indices=mask.idx)
  writeVolume(bv, paste0(PATH, "/", sids[i], "_encode_cor_scaled.nii"))
}
