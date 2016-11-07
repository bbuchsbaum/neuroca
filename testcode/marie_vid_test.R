library(tibble)
library(dplyr)
library(neuroim)
library(vegan)
library(ggrepel)
library(assertthat)
library(devtools)
library(energy)

devtools::load_all()

PATH <- "~/Dropbox/Brad/neuropls/video_marie/rds/"
group_des <- read.table(paste0(PATH, "/group_design.txt"), header=TRUE)
group_des$patient <- ifelse(group_des$sid %in% c(5001,5002), "patient", "control")
group_des <- subset(group_des, group != "old")
group_des$sid <- factor(group_des$sid)

clusnames <- paste0(PATH, levels(group_des$sid), "_nw_nlm_trial_data_run_all_blockavg_clus_256.nii")
clusvols <- lapply(clusnames, loadVolume)

loadMat <- function(sid) {
  message("loading", sid)
  dat <- readRDS(paste0(PATH, sid, "_nlm_trial_data_all_K256.rds"))
  X <- dat$X
  colnames(X) <- paste0("ROI_", 1:ncol(X))
  des <- tibble::as_data_frame(dat$des)
  des <- add_column(des, sid=rep(sid, nrow(des)))
  Xdes <- bind_cols(des,as_tibble(X))
  Xdes$condition <- ifelse(Xdes$Rating == "Video", "encode", "recall")
  roivars <- 10:ncol(Xdes)
  out <- Xdes %>% group_by(sid, condition, run, Video, time) %>% summarise_each(funs(mean), vars=roivars)
  
}

fill_clus <- function(clusvol, vals, K=256) {
  out <- fill(clusvol, cbind(1:K, vals))
}

fill_and_average <- function(cvols, lmat, K=256) {
  Reduce("+", lapply(1:length(cvols), function(i) {
    fill_clus(cvols[[i]], scale(lmat[,i]), K)
  }))/length(cvols)
}

plot_contributions <- function(contributions) {
  res <- lapply(1:ncol(contributions), function(i) {
    data.frame(sids=sids, comp=paste0("PC_", i), contrubution=contributions[,i])
  })
}


plot_scores <- function(cons, partial, labels, comps=c(1,2), sids, refcons=NULL) {
  if (!is.null(refcons)) {
    proc <- procrustes(refcons, cons, scale=FALSE)
    cons <- cons %*% proc$rotation
    partial <- lapply(partial, function(p) p %*% proc$rotation)
  }
  vnames <- paste0("PC", comps)
  cmat <- cons[,comps]
  df1 <- data.frame(cmat)
  names(df1) <- vnames
  df1$video <- labels
  
  omat <- do.call(rbind, lapply(partial, function(x) x[, comps]))
  row.names(omat) <- NULL
  df2 <- as.data.frame(omat)
  names(df2) <- vnames
  df2$sids <- rep(sids, each=nrow(cmat))
  df2$patient <- df2$sids
  df2$patient[df2$patient == 5001] <- "NC"
  df2$patient[df2$patient == 5002] <- "HC"
  df2$video <- rep(labels, length(partial))
  df3 <- subset(df2, patient %in% c("HC", "NC"))
  df2 <- subset(df2, !(patient %in% c("HC", "NC")))
  
  library(ggrepel)
 
  p1 <- ggplot(df2) +
    geom_point(aes_string(vnames[1], vnames[2], color = "video"), shape=16, size=5, data=df1) + 
    #stat_ellipse(aes(PC1, PC2, color=video)) +
    geom_point(aes_string(vnames[1], vnames[2], color = "video"), shape=16, size=1, alpha=.5, data=df2) +
    geom_point(aes_string(vnames[1], vnames[2], color = "video"), shape=18, size=2, data=df3) +
    geom_text_repel(aes_string(vnames[1], vnames[2], label = "video"), size=4, data=df1) +
    geom_text_repel(aes_string(vnames[1], vnames[2], label = "patient", color="video"), size=4, data=df3) +
    
    theme_classic(base_size = 24)
  
  print(p1)
  p1
  
}


hold_one_out <- function(Yl, Xl, ncomp=10, pfit=NULL) {
  res <- bind_rows(lapply(1:length(Yl), function(i) {
    print(i)
    Xtrain <- Xl[-i]
    Ytrain <- Yl[-i]
    
    Xtest <- Xl[[i]]
    Ytest <- Yl[[i]]
    
    if (is.null(pfit)) {
      print("refitting")
      pfit <- musubada(Ytrain, Xtrain, ncomp=10, center=TRUE, scale=TRUE, normalization="None")
    }
    
    creps <- code_replications(Ytest)
    p <- lapply(sort(unique(creps)), function(rnum) {
      ##browser()
      supX <- Xtest[creps != rnum,]
      supY <- Ytest[creps != rnum]
      
      yt <- Ytest[creps == rnum]
      
      
      predfun1 <- supplementary_predictor(pfit, supX, supY, ncomp=ncomp, type="prob")
      predfun2 <- supplementary_predictor(pfit, supX, supY, ncomp=ncomp, type="class")
      
      prob <- predfun1(Xtest[creps == rnum,])
      
      auc <- combinedAUC(prob, yt) + .5
      pwinner <- prob[cbind(1:nrow(prob), as.integer(yt))]
      cls <- predfun2(Xtest[creps==rnum,])
      
      list(prob=pwinner, auc=rep(auc,length(yt)), class=cls)
    })
    
    data_frame(auc=unlist(lapply(p, "[[", "auc")), prob=unlist(lapply(p, "[[", "prob")), 
               repnum=creps, pred=unlist(lapply(p, "[[", "class")), 
               observed=Ytest, sid=names(Yl)[i])
    
  }))
  
  res$acc <- as.integer(res$observed == res$pred)
  res$repnum <- as.numeric(res$repnum)
  res
  #res2 <- res %>% group_by(sid, repnum) %>% dplyr::summarise(prob=mean(prob), acc=mean(acc))
  #res3 <- res %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))
  #res4 <- res %>% group_by(repnum) %>% dplyr::summarise(prob=mean(prob), acc=mean(acc))
  
}


### set up data
sids <- levels(group_des$sid)
fulldes <- bind_rows(lapply(sids, loadMat))
fulldes$sid <- factor(fulldes$sid)
sids <- levels(fulldes$sid)

#### temporal analysis of encoding
encode_des <- filter(fulldes, condition == "encode")
roicols <- 6:ncol(encode_des)
XBlocks <- plyr::dlply(encode_des, "sid", function(x) as.matrix(x[, roicols]))

## normalize each image
XBlocks <- lapply(XBlocks, function(x) {
  t(apply(x,1, function(vals) vals/sum(vals^2)))
})

Ytime <- plyr::dlply(encode_des, "sid", function(x) ordered(x$time))
mfit_encode_time <- musubada(Ytime, XBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="None")
#r <- resample(mfit)
#tmp = project_copy(mfit, r$Y, r$X)
#mfitr <- musubada(r$Y, r$X, ncomp=10, center=TRUE, scale=FALSE, normalization="None")


## projection of each dataset on to first and second PCs
XProj <- project_cols(mfit_encode_time, ncomp=2)
XProj1 <- lapply(XProj, function(x) t(x[1,,]))
XProj2 <- lapply(XProj, function(x) t(x[2,,]))
red_encode_des <- filter(encode_des, time == 0)
Yvid_encode <- plyr::dlply(red_encode_des, "sid", function(x) x$Video)



## fit musubada for first and second projection
controls_idx <- 1:19
mfit_red_1 <- musubada(Yvid_encode, XProj1, ncomp=10, center=TRUE, scale=TRUE, normalization="None")
mfit_red_2 <- musubada(Yvid_encode, XProj2, ncomp=10, center=TRUE, scale=FALSE, normalization="None")

encode_pred <- hold_one_out(Yvid_encode, XProj1, ncomp=10)
encode_pred_s1 <- encode_pred %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))


#ncfit_red_1 <- musubada(Yvid[20], XProj1[20], ncomp=10, center=TRUE, scale=FALSE, normalization="None")
#hcfit_red_1 <- musubada(Yvid[21], XProj1[21], ncomp=10, center=TRUE, scale=FALSE, normalization="None")

## projection of each subject for each block
# tmp = lapply(1:length(sids), function(i) {
#   d1 <- filter(red_encode_des, sid == sids[i])
#   blocks <- unique(d1$run)
#   Xcur <- XProj1[[i]]
#   m=do.call(cbind, lapply(blocks, function(j) {
#     idx <- which(d1$run == j)
#     supY <- d1$Video[idx]
#     proj <- predict(mfit_red_1, Xcur[idx,], type="scores", table_index=j)
#     apply(proj,2,function(vals) sum(vals^2))
#     
#   }))
#   
# })


## plot space
p1_2 <- plot_scores(mfit_red_1$scores,mfit_red_1$partial_scores, levels(Yvid[[1]]), comps=c(1,2), sids=sids)
p3_4 <- plot_scores(mfit_red_1$scores,mfit_red_1$partial_scores, levels(Yvid[[1]]), comps=c(3,4), sid=sids)
#p2 <- plot_scores(mfit_red_2$scores,mfit_red_2$partial_scores, levels(Yvid[[1]]), refcons=mfit_red_1$scores)


#### temporal analysis of recall
recall_des <- filter(fulldes, condition == "recall")
XBlocks <- plyr::dlply(recall_des, "sid", function(x) as.matrix(x[, roicols]))
## normalize each image
XBlocks <- lapply(XBlocks, function(x) {
  t(apply(x,1, function(vals) vals/sum(vals^2)))
})

Ytime <- plyr::dlply(recall_des, "sid", function(x) ordered(x$time))
mfit_recall <- musubada(Ytime, XBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="None")

## projection of each dataset on to first and second PCs
XProj_recall <- project_cols(mfit_recall, ncomp=2)
XProj1_recall <- lapply(XProj_recall, function(x) t(x[1,,]))
XProj2_recall <- lapply(XProj_recall, function(x) t(x[2,,]))
red_recall_des <- filter(recall_des, time == 0)
Yvid_recall <- plyr::dlply(red_recall_des, "sid", function(x) x$Video)

mfit_recall_vid_1 <- musubada(Yvid_recall, XProj1_recall, ncomp=10, center=TRUE, scale=FALSE, normalization="None")
mfit_recall_vid_2 <- musubada(Yvid_recall, XProj2_recall, ncomp=10, center=TRUE, scale=FALSE, normalization="None")

## embedding recall in encoding model
recall_embed <- hold_one_out(Yvid_recall, XProj1_recall, ncomp=10, pfit=mfit_red_1)
recall_embed_s1 <- recall_embed %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))

## recall in recall
recall_pred <- hold_one_out(Yvid_recall, XProj1_recall, ncomp=10)
recall_pred_s1 <- recall_pred %>% group_by(sid) %>% dplyr::summarise(auc=mean(auc), prob=mean(prob), acc=mean(acc))


## projection of recall on encoding fit
auc=unlist(lapply(1:length(Yvid_recall), function(i) {
  p1 <- predict(mfit_red_1, XProj1_recall[[i]], type="prob", ncomp=10, table_index=i)
  #sum(p1 == Yvid[[i]])/length(Yvid[[i]])
  combinedAUC(p1, Yvid_recall[[i]]) + .5
}))

## ratio of cross-decode and embedded decode
auc_ratio <- auc/recall_embed_s1$auc

## individual cross-projections
auc2 <- unlist(lapply(1:length(Yvid_recall), function(i) {
  #musubada(Yvid[20], XProj1[20], ncomp=10, center=TRUE, scale=FALSE, normalization="None")
  mf <- musubada(Yvid_encode[i], XProj1[i], ncomp=10, center=TRUE, scale=FALSE, normalization="None")
  p1 <- predict(mf, XProj1_recall[[i]], type="prob", ncomp=10, table_index=1)
  #sum(p1 == Yvid[[i]])/length(Yvid[[i]])
  combinedAUC(p1, Yvid_recall[[i]]) + .5
}))


CrawfordHowell <- function(case, control){
  tval <- (case - mean(control)) / (sd(control)*sqrt((length(control)+1) / length(control)))
  degfree <- length(control)-1
  pval <- 2*(1-pt(abs(tval), df=degfree)) #two-tailed p-value
  result <- data.frame(t = tval, df = degfree, p=pval)
  return(result)
}


p1_2 <- plot_scores(mfit_recall_vid_1$scores,mfit_recall_vid_1$partial_scores, levels(Yvid[[1]]))
p3_4 <- plot_scores(mfit_recall_vid_1$scores,mfit_recall_vid_1$partial_scores, levels(Yvid[[1]]), comps=c(1,3))
lds1 <- do.call(cbind, lapply(1:length(XProj1), function(i) loadings(mfit_red_1, i,1)))
lds2 <- do.call(cbind, lapply(1:length(XProj1), function(i) loadings(mfit_red_1, i,2)))
lds3 <- do.call(cbind, lapply(1:length(XProj1), function(i) loadings(mfit_red_1, i,3)))
writeVolume(fill_and_average(clusvols, lds1), paste0(PATH, "/", "encode_avg_lds_video_time1_pc1.nii"))
writeVolume(fill_and_average(clusvols, lds2), paste0(PATH, "/", "encode_avg_lds_video_time1_pc2.nii"))
writeVolume(fill_and_average(clusvols, lds3), paste0(PATH, "/", "encode_avg_lds_video_time1_pc3.nii"))
lds1 <- do.call(cbind, lapply(1:length(XProj1), function(i) loadings(mfit_red_2, i,1)))
lds2 <- do.call(cbind, lapply(1:length(XProj1), function(i) loadings(mfit_red_2, i,2)))
lds3 <- do.call(cbind, lapply(1:length(XProj1), function(i) loadings(mfit_red_2, i,3)))
writeVolume(fill_and_average(clusvols, lds1), paste0(PATH, "/", "encode_avg_lds_video_time2_pc1.nii"))
writeVolume(fill_and_average(clusvols, lds2), paste0(PATH, "/", "encode_avg_lds_video_time2_pc2.nii"))
writeVolume(fill_and_average(clusvols, lds3), paste0(PATH, "/", "encode_avg_lds_video_time2_pc3.nii"))


