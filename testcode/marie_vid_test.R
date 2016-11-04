library(tibble)
library(dplyr)
library(neuroim)
library(vegan)
library(ggrepel)

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


plot_scores <- function(cons, partial, labels, comps=c(1,2), refcons=NULL) {
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
    geom_point(aes_string(vnames[1], vnames[2], color = "video"), shape=16, size=1, alpha=.5, data=df2) +
    geom_point(aes_string(vnames[1], vnames[2], color = "video"), shape=18, size=2, data=df3) +
    geom_text_repel(aes_string(vnames[1], vnames[2], label = "video"), size=4, data=df1) +
    geom_text_repel(aes_string(vnames[1], vnames[2], label = "patient", color="video"), size=4, data=df3) +
    
    theme_classic(base_size = 24)
  
  print(p1)
  p1
  
}


hold_one_out <- function(Yl, Xl) {
  res <- bind_rows(lapply(1:length(Yl), function(i) {
    print(i)
    Xtrain <- Xl[-i]
    Ytrain <- Yl[-i]
    
    Xtest <- Xl[[i]]
    Ytest <- Yl[[i]]
    mfit <- musubada(Ytrain, Xtrain, ncomp=10, center=TRUE, scale=FALSE, normalization="DCor")
    creps <- code_replications(Ytest)
    p <- lapply(sort(unique(creps)), function(rnum) {
      supX <- Xtest[creps != rnum,]
      supY <- Ytest[creps != rnum]
      
      yt <- Ytest[creps == rnum]
      predfun1 <- supplementary_predictor(mfit, supX, supY, ncomp=7, type="prob")
      predfun2 <- supplementary_predictor(mfit, supX, supY, ncomp=7, type="class")
      d <- predfun1(Xtest[creps == rnum,])
      d <- apply(d,1, function(x) x/max(x))
      md <- d[cbind(1:nrow(d), as.integer(yt))]
      cls <- predfun2(Xtest[creps==rnum,])
      
      list(md=md, class=cls)
    })
    
    data_frame(cosine=unlist(lapply(p, "[[", "md")), repnum=creps, pred=unlist(lapply(p, "[[", "class")), observed=Ytest, sid=names(Yl)[i])
    
  }))
  
  res$acc <- as.integer(res$observed == res$pred)
  res2 <- res %>% group_by(sid, repnum) %>% dplyr::summarise(cosine=mean(cosine), acc=mean(acc))
  res2$repnum = as.numeric(res2$repnum)
  allpred[cbind(1:length(labels), as.integer(labels))]
  
  #combinedAUC(allpred, labels)
}





sids <- levels(group_des$sid)
fulldes <- bind_rows(lapply(sids, loadMat))
fulldes$sid <- factor(fulldes$sid)
sids <- levels(fulldes$sid)

#### temporal analysis of encoding
encode_des <- filter(fulldes, Rating == "Video")
roicols <- 10:ncol(encode_des)
XBlocks <- plyr::dlply(encode_des, "sid", function(x) as.matrix(x[, roicols]))
XBlocks <- lapply(XBlocks, function(x) x/sd(x))
Ytime <- plyr::dlply(encode_des, "sid", function(x) ordered(x$time))
mfit <- musubada(Ytime, XBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="None")


## projection of each dataset on to first and second PCs
XProj <- project_cols(mfit, ncomp=2)
XProj1 <- lapply(XProj, function(x) t(x[1,,]))
XProj2 <- lapply(XProj, function(x) t(x[2,,]))
red_encode_des <- filter(encode_des, time == 0)
Yvid <- plyr::dlply(red_encode_des, "sid", function(x) x$Video)
mfit_red_1 <- musubada(Yvid, XProj1, ncomp=10, center=TRUE, scale=FALSE, normalization="None")
mfit_red_2 <- musubada(Yvid, XProj2, ncomp=10, center=TRUE, scale=FALSE, normalization="None")




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



p1_2 <- plot_scores(mfit_red_1$scores,mfit_red_1$partial_scores, levels(Yvid[[1]]))
p3_4 <- plot_scores(mfit_red_1$scores,mfit_red_1$partial_scores, levels(Yvid[[1]]), comps=c(3,4))

p2 <- plot_scores(mfit_red_2$scores,mfit_red_2$partial_scores, levels(Yvid[[1]]), refcons=mfit_red_1$scores)

#encode_trials <- encode_des %>% group_by(sid, repnum, Video) %>% do(
#  mat=project_cols(mfit, as.matrix(select(., roicols)), table_names=.$sid[1], ncomp=2)
#)


#### temporal analysis of recall
recall_des <- filter(fulldes, Rating != "Video")
XBlocks <- plyr::dlply(recall_des, "sid", function(x) as.matrix(x[, roicols]))
Ytime <- plyr::dlply(recall_des, "sid", function(x) ordered(x$time))
mfit2 <- musubada(Ytime, XBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="MFA")





runsub <- function(s, tsel) {
  
  labels <- (filter(alldat, sid == s, time == 0 & Rating == "Video") %>% select(Video))$Video
  
  runs <- unique(filter(alldat, sid == s)$run)
  allpred <- do.call(rbind, (lapply(runs, function(bnum) {
    print(bnum)
    Xtrain <- lapply(tsel, function(i) as.matrix(filter(alldat, sid == s & time == i & Rating == "Video" & run != bnum) %>% select(roicols)))
    Xtest <- lapply(tsel, function(i) as.matrix(filter(alldat, sid == s & time == i  & Rating == "Video" & run == bnum) %>% select(roicols)))
  
    Ytrain <- (filter(alldat, sid == s, time == 0 & Rating == "Video" & run != bnum) %>% select(Video))$Video
    Ytest <- (filter(alldat, sid == s, time == 0 & Rating == "Video" & run == bnum) %>% select(Video))$Video
  
    mfit <- musubada(Ytrain, Xtrain, ncomp=10, center=TRUE, scale=FALSE, normalization="None")
  
    Xtest <- do.call(cbind, Xtest)
    pred <- predict.musubada_result(mfit, Xtest, table_index=1:mfit$ntables, type="prob")
    
  
  })))
  
  allpred[cbind(1:length(labels), as.integer(labels))]
  
  #combinedAUC(allpred, labels)
}


make_preds <- function(tsel) {

  Xsel <- lapply(Xs, function(x) x[, t1 %in% tsel])

  pall <- sapply(1:length(Xs), function(i) {
  
    X <- Xsel[-i]
    Y <- Ys[-i]
  
    res <- musubada(Y, X, ncomp=10, center=TRUE, scale=FALSE, normalization="None")
    Xsup <- Xsel[[i]]
    run <- blocks[[i]]
    vid <- Ys[[i]]

    message("subject", i)
    preds <- do.call(rbind, lapply(sort(unique(run)), function(j) {
      message("run", j)
      keep <- run != j
      heldout <- run == j
      Xembed <- Xsup[keep,]
      Yembed <- vid[keep]
  
      predfun <- supplementary_predictor(res, Xembed, Yembed, type="class", ncomp=10)
      gmeans <- group_means(vid[heldout], Xsup[heldout,])
    
      p <- as.character(predfun(gmeans))
      data.frame(obs=as.character(row.names(gmeans)), pred=as.character(predfun(gmeans)))
    }))
  
    sum(as.character(preds[,1]) == as.character(preds[,2]))/nrow(preds)
  
  })
}

allpred <- lapply(sort(unique(t1)), make_preds)

pmid <- make_preds(sort(unique(t1))[4:7])
plot(pmid)


