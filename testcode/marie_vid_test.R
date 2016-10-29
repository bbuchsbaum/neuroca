sids <- c("1001", "1002", "1003", "1007", "1008", "1016", "3008", "5001", "5002")

PATH <- "~/Dropbox/Brad/neuropls/video_marie/rds/"

loadMat <- function(sid) {
  message("loading", sid)
  dat <- readRDS(paste0(PATH, sid, "_nlm_trial_data_all_K256.rds"))
  idx <- which(dat$design$Rating == "Video")
  dat$X <- dat$X[idx,]
  dat$design <- dat$design[idx,]
  dat
  # mlist <- lapply(1:length(fnames), function(i) {
  #   d <- subset(des, run == i)
  #   m <- t(as.matrix(loadVector(fnames[i]))[mask.idx,])
  #   folded <- do.call(rbind, lapply(levels(d$Video), function(vid) {
  #     vidx <- which(d$Video == vid & d$condition == "Encod")
  #     vec <- unlist(lapply(vidx, function(j) m[j,]))
  #   }))
  #   
  # })
  # 
}

fold_matrix <- function(X, des) {
  Xs <- split_matrix(X, des$Video)
  dsplit <- split(des, des$Video)
  
  ret <- lapply(1:length(dsplit), function(i) {
    do.call(cbind, split_matrix(Xs[[i]], factor(dsplit[[i]]$time)))
  })
  
  nreps <- sapply(ret,nrow)
  ret <- do.call(rbind, ret)
  Y <- rep(levels(des$Video), nreps)
  
  list(Y=Y, X=ret, time=rep(sort(unique(des$time)), each=ncol(X)))
}


alldat <- lapply(sids, loadMat)
Xlist <- lapply(alldat, "[[","X")
dlist <- lapply(alldat, function(x) x$design)

ret <- lapply(1:length(Xlist), function(i) fold_matrix(Xlist[[i]], dlist[[i]]))
Xs <- lapply(ret, "[[", "X")
Ys <- lapply(ret, "[[", "Y")
time <- lapply(ret, "[[", "time")



