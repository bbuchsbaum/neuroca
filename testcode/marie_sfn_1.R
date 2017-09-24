library(neuroim)
library(tibble)
library(dplyr)
library(vegan)
library(ggrepel)
library(assertthat)
library(devtools)
library(energy)
library(tidyr)
devtools::load_all()

PATH <- "~/Dropbox/Brad/neuropls/video_marie/SFN/"
sids <- scan(paste0(PATH, "/sids.txt"))

mask <- loadVolume(paste0(PATH, "/clus_mask.nii"))
mask.idx <- which(mask>0)

for (sid in sids) {
  print(sid)
  des <- read.table(paste0(PATH, "/", sid, "_design.txt"), header=TRUE)
  fname <- paste0(PATH, "/", sid, "_nwall_betas.nii.gz")
  vec <- as.matrix(loadVector(fname))[mask.idx,]
  
  #out <- list(X=vec, design=des, sid=sid)
  #saveRDS(out, file=paste0(PATH, "/", sid, "_all_betas.RDS"))
  
  fac <- interaction(des$Video, des$run, des$Condition)
  gmeans <- group_means(fac, t(vec))
  
  df1 <- data.frame(fac=levels(fac))
  df1 <- df1 %>% separate(fac, c("Video", "run", "Condition"))
  out2 <- list(X=gmeans, design=df1, sid=sid)
  saveRDS(out2, file=paste0(PATH, "/", sid, "_block_betas.RDS"))
  
}