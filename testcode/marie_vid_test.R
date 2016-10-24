sids <- c("1001", "1007", "1008", "1016", "3008", "5001", "5002")

PATH <- "~/Dropbox/Brad/neuropls/video_marie/"
bold <- loadVolume(paste0(PATH, "/mean_bold.nii"))
mask.idx <- which(bold > 600)

loadMat <- function(sid) {
  message("loading", sid)
  des <- read.table(paste0(PATH, sid, "/trial_data_all_blockavg.txt"), header=TRUE)
  fnames <- list.files(paste0(PATH, sid), "nw.*nii.gz", full.names=TRUE)
  mlist <- lapply(1:length(fnames), function(i) {
    d <- subset(des, run == i)
    m <- t(as.matrix(loadVector(fnames[i]))[mask.idx,])
    folded <- do.call(rbind, lapply(levels(d$Video), function(vid) {
      vidx <- which(d$Video == vid & d$condition == "Encod")
      vec <- unlist(lapply(vidx, function(j) m[j,]))
    }))
    
  })
  
  M <- do.call(rbind, mlist)
  dout <- subset(des, time==0 & condition == "Encod")
  list(des=dout, M=M)
  
}

allMat <- lapply(sids, loadMat)
Xlist <- lapply(allMat, "[[","M")
fulldes <- subset(allMat[[1]]$des, condition == "Encod")
Y <- fulldes$Video

