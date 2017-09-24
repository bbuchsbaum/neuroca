library(neuroim)
library(tibble)
library(dplyr)
library(vegan)
library(ggrepel)
library(assertthat)
library(devtools)
library(energy)
library(sda)
library(parallel)

ncores <- 4

devtools::load_all()

PATH <- "~/Dropbox/Brad/neuropls/video_marie/SFN/"
sids <- scan(paste0(PATH, "/sids.txt"))

mask <- loadVolume(paste0(PATH, "/clus_mask.nii"))
mask.idx <- which(mask>0)

# for (sid in sids) {
#   print(sid)
#   x <- readRDS(paste0(PATH, "/", sid, "_block_betas.RDS"))
#   names(x) <- c("Video", "run", "condition")
#   saveRDS(x, paste0(paste0(PATH, "/", sid, "_block_betas.RDS")))
# }

predict.sda2 <- function (object, Xtest, verbose = TRUE, ...) {
  if (missing(object)) {
    stop("A sda fit object must be supplied.")
  }
  if (missing(Xtest)) {
    stop("A new data to predict must be supplied.")
  }
  if (!is.matrix(Xtest)) 
    stop("Test data must be given as matrix!")
  ntest = nrow(Xtest)
  alpha = object$alpha
  cl.count = length(alpha)
  if (ncol(Xtest) != ncol(object$beta)) 
    stop("Different number of predictors in sda object (", 
         ncol(object$beta), ") and in test data (", ncol(Xtest), 
         ")", sep = "")
  beta = object$beta
  if (verbose) 
    cat("Prediction uses", ncol(beta), "features.\n")
  probs = t(tcrossprod(beta, Xtest) + alpha)
}

combinedAUC <- function(Pred, Obs) {
  mean(sapply(1:ncol(Pred), function(i) {
    lev <- levels(Obs)[i]
    pos <- Obs == lev
    pclass <- Pred[,i]
    pother <- rowMeans(Pred[,-i,drop=FALSE])
    Metrics::auc(as.numeric(pos), pclass - pother)-.5
  }))
}


#des <- readRDS(paste0(PATH, "/", 5001, "_block_betas.RDS"))

fit_sda <- function(des) {
  out <- parallel::mclapply(sort(unique(des$design$run)), function(rnum) {
    heldout <- which(des$design$run == rnum & des$design$Condition == "Encod")
    keep <- des$design$run != rnum & des$design$Condition == "Encod"
  
    X <- t(des$X)
    Xtrain <- X[keep,]
    Xtest <- X[heldout,]
  
    Ytrain <- des$design$Video[keep]
    Ytest <- des$design$Video[heldout]
  
    sda.1 <- sda(Xtrain, Ytrain)
    pred1 <- predict.sda2(sda.1, Xtest) 
    pred1a <- predict(sda.1, Xtest)

    heldout2 <- which(des$design$run == rnum & des$design$Condition == "Recall")
    Ytest2 <- des$design$Video[heldout2]
    Xtest2 <- X[heldout2,]
    pred2 <- predict.sda2(sda.1, Xtest2) 
    pred2a <- predict(sda.1, Xtest2)
    
    keep2 <- des$design$run != rnum & des$design$Condition == "Recall"
    
    Xtrain2 <- X[keep2,]
    Ytrain2 <- des$design$Video[keep2]
    sda.2 <- sda(Xtrain2, Ytrain2)
    pred3 <- predict.sda2(sda.2, Xtest2) 
    pred3a <- predict(sda.2, Xtest2)
  
    list(encode_pred=pred1, cross_pred=pred2, recall_pred=pred3, 
         encode_class=pred1a$class, cross_class=pred2a$class, recall_class=pred3a$class)
  
  }, mc.cores=ncores)
  
  Yrecall <- des$design$Video[des$design$Condition == "Recall"]
  Yencode <- des$design$Video[des$design$Condition == "Encod"]

  encode_mat <- do.call(rbind, lapply(out, "[[", "encode_pred"))
  cross_mat <- do.call(rbind, lapply(out, "[[", "cross_pred"))
  recall_mat <- do.call(rbind, lapply(out, "[[", "recall_pred"))
  
  encode_cls <- unlist(lapply(out, "[[", "encode_class"))
  cross_cls <- unlist(lapply(out, "[[", "cross_class"))
  recall_cls <- unlist(lapply(out, "[[", "recall_class"))
  
  
  auc_encode <- combinedAUC(encode_mat, factor(des$design$Video[des$design$Condition == "Encod"]))
  auc_cross <- combinedAUC(cross_mat, factor(des$design$Video[des$design$Condition == "Recall"]))
  auc_recall <- combinedAUC(recall_mat, factor(des$design$Video[des$design$Condition == "Recall"]))
  
  list(auc_cross=auc_cross, auc_recall=auc_recall, auc_encode=auc_encode, 
       acc_encode=sum(encode_cls==Yencode)/length(Yencode),
       acc_cross=sum(cross_cls==Yrecall)/length(Yrecall),
       acc_recall=sum(recall_cls==Yrecall)/length(Yrecall),
       
       sid=des$sid)
}

res <- lapply(sids, function(sid) {
  print(sid)
  des <- readRDS(paste0(PATH, "/", sid, "_all_betas.RDS"))
  ret <- fit_sda(des)
  print(ret)
  ret
})

df2 <- do.call(rbind, lapply(res, function(x) {
  data.frame(sid=x$sid, 
             auc=c(x$auc_encode, x$auc_recall,x$auc_cross),
             acc=c(x$acc_encode, x$acc_recall, x$acc_cross),
             classifier=c("encoding", "recall", "cross-decode"))
             
  }))


df2 <- as.data.frame(df2)
df2$group <- ifelse(df2$sid == 5001, "NC", "control")
df2$group[df2$sid == 5002] <- "HC"
df2$group = factor(df2$group)
df2$classifier <- factor(df2$classifier, levels=c("encoding", "recall", "cross-decode"))
df2$auc <- df2$auc + .5

pdf(file=paste0(PATH, "within_suject.pdf"))
theme_set(theme_gray(base_size = 18))
ggplot(df2, aes(x = factor(classifier), fill = factor(group), y = acc)) +
  geom_dotplot(binaxis = "y", stackdir = "center") + labs(x= "Classifier Type", y="Accuracy", fill="Group")
dev.off()

write.table(df2, paste0(PATH, "/within_subject_results.txt"))

