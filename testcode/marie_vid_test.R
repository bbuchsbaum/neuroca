library(tibble)
library(dplyr)
library(neuroim)

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
    fill_clus(cvols[[i]], lmat[,i], K)
  }))/length(cvols)
}


fulldes <- bind_rows(lapply(group_des$sids, loadMat))

#### temporal analysis of encoding
encode_des <- filter(fulldes, Rating == "Video")
roicols <- 10:ncol(encode_des)
XBlocks <- plyr::dlply(encode_des, "sid", function(x) as.matrix(x[, roicols]))
Ytime <- plyr::dlply(encode_des, "sid", function(x) ordered(x$time))
mfit <- musubada(Ytime, XBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="MFA")
avg_lds1 <- fill_and_average(clusvols, do.call(cbind, lapply(1:mfit$ntables, function(i) loadings(mfit, i, 1))))
avg_lds2 <- fill_and_average(clusvols, do.call(cbind, lapply(1:mfit$ntables, function(i) loadings(mfit, i, 1))))
## this yields two components
## next, project data on comp1 and 2. compare NC and HC to group.


#### temporal analysis of recall
recall_des <- filter(fulldes, Rating != "Video")
XBlocks <- plyr::dlply(recall_des, "sid", function(x) as.matrix(x[, roicols]))
Ytime <- plyr::dlply(recall_des, "sid", function(x) ordered(x$time))
mfit2 <- musubada(Ytime, XBlocks, ncomp=10, center=TRUE, scale=FALSE, normalization="MFA")
avg_lds1 <- fill_and_average(clusvols, do.call(cbind, lapply(1:mfit$ntables, function(i) loadings(mfit2, i, 1))))
avg_lds2 <- fill_and_average(clusvols, do.call(cbind, lapply(1:mfit$ntables, function(i) loadings(mfit2, i, 2))))
avg_lds3 <- fill_and_average(clusvols, do.call(cbind, lapply(1:mfit$ntables, function(i) loadings(mfit2, i, 3))))

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


