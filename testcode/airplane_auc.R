library(devtools)
library(neuroim)

load_all()


path <- "/Users/bbuchsbaum/analysis/airplane_paper/searchAUC/"
mask <- loadVolume(paste0(path, "group_global_mask2.nii"))
mask.idx <- which(mask>0)

sids <- read.table(paste0(path, "/all.sids"), header=TRUE)

alldes <- do.call(rbind, lapply(1:nrow(sids), function(i) {
  s <- sids$sid[i]
  
  encod <- list.files(paste0(path, s, "/Searchlight/searchlight_encode/"), "nw.*nii", full.names=TRUE)
  recall <- list.files(paste0(path, s, "/Searchlight/searchlight_recall/"), "nw.*nii", full.names=TRUE)
  cross <- list.files(paste0(path, s, "/Searchlight/searchlight_cross/"), "nw.*nii", full.names=TRUE)
  
  data.frame(sid=s, file=c(encod, recall,cross), condition=rep(c("encode", "recall", "cross"), each=length(encod)),
             group=sids$group[i], run=rep(1:length(encod),3))
}))

alldes$run <- factor(alldes$run)
des_encode <- droplevels(subset(alldes, condition=="encode" & group != "patient"))
des_recall <- droplevels(subset(alldes, condition=="recall" & group != "patient"))
des_cross <- droplevels(subset(alldes, condition=="cross"  & group != "patient"))

volmat_enc <- do.call(rbind, lapply(as.character(des_encode$file), function(fn) {
  loadVolume(fn)[mask.idx]
}))

des_encode$group_run <- interaction(des_encode$run, des_encode$group)
bres <- bada(des_encode$run, volmat_enc, S=des_encode$sid, center=FALSE)
boot_res <- bootstrap(bres, nboot=100)
zboot1 <- boot_res$zboot_loadings[,1]
zout <- BrainVolume(zboot1, space(mask), indices=mask.idx)
writeVolume(zout, paste0(path, "zboot1_encode_run.nii"))

volmat_recall<- do.call(rbind, lapply(as.character(des_recall$file), function(fn) {
  loadVolume(fn)[mask.idx]
}))

des_recall$group_run <- interaction(des_recall$run, des_recall$group)
bres <- bada(des_recall$run, volmat_recall, S=des_recall$sid)
Xresid <- residualize(~ run + group, volmat_recall, des_recall)
bres2 <- bada(des_recall$group_run, Xresid, S=des_recall$sid)
boot_res <- bootstrap(bres2, nboot=100)

zboot1 <- boot_res$zboot_loadings[,1]
zout <- BrainVolume(zboot1, space(mask), indices=mask.idx)
writeVolume(zout, paste0(path, "zboot1_recall_group_run.nii"))

bres <- bada(des_recall$run, volmat_recall, S=des_recall$sid)
boot_res <- bootstrap(bres, nboot=100)
zboot1 <- boot_res$zboot_loadings[,1]
zout <- BrainVolume(zboot1, space(mask), indices=mask.idx)
writeVolume(zout, paste0(path, "zboot1_recall_both_run.nii"))

zboot2 <- boot_res$zboot_loadings[,2]
zout <- BrainVolume(zboot2, space(mask), indices=mask.idx)
writeVolume(zout, paste0(path, "zboot2_recall_both_run.nii"))


volmat_cross<- do.call(rbind, lapply(as.character(des_cross$file), function(fn) {
  loadVolume(fn)[mask.idx]
}))

bres <- bada(des_cross$run, volmat_cross, S=des_cross$sid)
boot_res <- bootstrap(bres, nboot=100)
zboot1 <- boot_res$zboot_loadings[,1]
zout <- BrainVolume(zboot1, space(mask), indices=mask.idx)
writeVolume(zout, paste0(path, "zboot1_cross_both_run.nii"))




des_recall_old <- droplevels(subset(des_recall, group=="old"))
volmat_recall_old <- volmat_recall[des_recall$group == "old",]
bres_old_run <- bada(des_recall_old$run, volmat_recall_old, S=des_recall_old$sid)
boot_res <- bootstrap(bres_old_run, nboot=100)
zboot1 <- boot_res$zboot_loadings[,1]
zout <- BrainVolume(zboot1, space(mask), indices=mask.idx)
writeVolume(zout, paste0(path, "zboot1_recall_old_run.nii"))


