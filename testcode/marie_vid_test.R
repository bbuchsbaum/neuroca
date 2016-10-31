sids <- c("1001", "1002", "1003", "1007", "1008", "1009", "1010",
           "1011", "1013", "1014", "1015", "1016", "13001", "13002",
           "3003", "3004", "3005", "3006", "3008", "2004", "2005", "2006", "2007",
            "2009", "2010", "2013", "2014", "5001", "5002")

PATH <- "~/Dropbox/Brad/neuropls/video_marie/rds/"

loadMat <- function(sid) {
  message("loading", sid)
  dat <- readRDS(paste0(PATH, sid, "_nlm_trial_data_all_K256.rds"))
  X <- dat$X
  colnames(X) <- paste0("ROI_", 1:ncol(X))
  des <- tibble::as_data_frame(dat$des)
  des <- add_column(des, sid=rep(sid, nrow(des)))
  Xdes <- bind_cols(des,as_tibble(X))
  
}

# fold_matrix <- function(X, des) {
#   Xs <- split_matrix(X, des$Video)
#   dsplit <- split(des, des$Video)
#   
#   ret <- lapply(1:length(dsplit), function(i) {
#     do.call(cbind, split_matrix(Xs[[i]], factor(dsplit[[i]]$time)))
#   })
#   
#   nreps <- sapply(ret,nrow)
#   ret <- do.call(rbind, ret)
#   Y <- rep(levels(des$Video), nreps)
#   
#   run <- unlist(lapply(dsplit, function(d) subset(d, time==0)$run))
#   
#   list(Y=Y, X=ret, run=run, time=rep(sort(unique(des$time)), each=ncol(X)))
# }


alldat <- bind_rows(lapply(sids, loadMat))
roicols <- 10:ncol(alldat)

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


