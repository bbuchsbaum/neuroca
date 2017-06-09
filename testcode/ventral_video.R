path <- "/Users/bbuchsbaum/analysis/hyper/ventral_surface/"

sids <- scan(paste0(path, "sids"), "")


load_mat <- function(fname) {
  tmp <- readRDS(paste0(path, fname))
  d1 <- tmp$left$data
  d2 <- tmp$right$data
  d3 <- cbind(t(d1), t(d2))
  list(mat=d3, design=tmp$design)
}

  

vdat <- lapply(sids, function(sid) {
  v1 <- load_mat(paste0(sid, "_s5_ventral_video_1.RDS"))
  v2 <- load_mat(paste0(sid, "_s5_ventral_video_2.RDS"))
  v3 <- load_mat(paste0(sid, "_s5_ventral_video_3.RDS"))
  des <- rbind(v1$design, v2$design, v3$design)
  des$run <- factor(rep(1:3, each=364))
  
  mat <- rbind(v1$mat, v2$mat, v3$mat)
  
  keep <- apply(mat,2, function(x) sum(x==0)) == 0
  mat <- mat[,keep]
  
  rmat <- resid(lsfit(model.matrix(~ factor(des$run)), mat, intercept=FALSE))
  
  mat <- t(scale(t(rmat)))
  
  ncomp <- fast_estim_ncomp(scale(mat), 1, 300)$bestcomp
  
  
  list(sid=sid, mat=mat, design=des, ncomp=ncomp, keep=keep)
})

ndat <- lapply(vdat, function(vd) {
  sid <- vd$sid
  v1 <- load_mat(paste0(sid, "_s5_ventral_nback.RDS"))
  v1$design$LabelVersion <- factor(paste0(as.character(v1$design$Label), "_", v1$design$Version))
  keep <- vd$keep
  mat <- v1$mat[, keep]
  
  rmat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
  
  
  
  mat <- t(scale(t(rmat)))
  mat <- group_means(v1$design$LabelVersion, mat)
  
  mat_miss <- matrix(NA, nrow(mat)*2, ncol(mat))
  
  levs <- levels(v1$design$LabelVersion)
  labs <- levels(v1$design$Label)
  ver <- sapply(strsplit(levs, "_"), function(x) {
    x[length(x)]
  })
  
  for (i in 1:length(ver)) {
    if (ver[i] == "1") {
      mat_miss[i*2-1,] <- mat[i,]
    } else {
      mat_miss[i*2,] <- mat[i,]
    }
  }
  
  ncomp <- fast_estim_ncomp(scale(mat), 1, 89)$bestcomp
  
  des <- data.frame(label=paste0(rep(labs, each=2), "_", rep(1:2, 90)), version=rep(1:2, 90))
  
  list(sid=sid, mat=mat_miss, design=des, ncomp=ncomp, keep=keep, len=ncol(mat))
})



idat <- lapply(vdat, function(vd) {
  sid <- vd$sid
  v1 <- load_mat(paste0(sid, "_s5_ventral_imagery.RDS"))
  v1$design$LabelVersion <- factor(paste0(as.character(v1$design$label), "_", v1$design$version))
  keep <- vd$keep
  mat <- v1$mat[, keep]
  
  rmat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
  mat <- t(scale(t(rmat)))
  list(sid=sid, mat=mat, design=v1$design, keep=keep, len=ncol(mat))
  
})




library(softImpute)
Xnback <- do.call(cbind, lapply(ndat, "[[", "mat"))
imp_fit <- softImpute(Xnback, rank.max=25, type="als", trace.it=TRUE)
Xcomplete <- complete(Xnback,imp_fit)
lens <- lapply(ndat, "[[", "len")
groups <- rep(1:length(lens), unlist(lens))

#Xbm <- matrix_to_block_matrix(Xcomplete, groups)
#Xl <- as.list(Xbm)
#Y <- ndat[[1]]$design$label
#mbres <- musubada(Y, Xl, ncomp=20, normalization="MFA", rank_k=120)
#df1 <- data.frame(pc1=mbres$scores[,1], pc2=mbres$scores[,2], pc3=mbres$scores[,3], pc4=mbres$scores[,4], pc5=mbres$scores[,5], labels=row.names(mbres$scores))
#ggplot(df1) + geom_point(aes(pc1,pc2), color="red") + geom_text_repel(aes(pc1,pc2, label = labels)) + theme_classic(base_size = 10)


zerostr <- function(vals) {
  ndigits <- log(max(vals), 10) + 1
  nzeros <- ndigits - as.integer(log(vals,10)) -1
  prefix <- sapply(nzeros, function(nz) paste(rep("0", times=nz), collapse=""))
  paste(prefix, vals, sep="")  
}



Xvideo <- do.call(cbind, lapply(vdat, "[[", "mat"))
Xall <- rbind(Xvideo, Xnback)
imp_fit <- softImpute(Xall, rank.max=25, type="als", trace.it=TRUE)
Xcomplete <- complete(Xall,imp_fit)
lens <- lapply(ndat, "[[", "len")
groups <- rep(1:length(lens), unlist(lens))
Xbm <- matrix_to_block_matrix(Xcomplete, groups)
Xl <- as.list(Xbm)
Y <- factor(c(zerostr(1:1092), as.character(ndat[[1]]$design$label)))
nback_ind <- 1093:(1093+179)

mbres <- musubada(Y, Xl, ncomp=20, normalization="MFA", rank_k=200)
library(ggrepel)
df1 <- data.frame(pc1=mbres$scores[nback_ind,1], pc2=mbres$scores[nback_ind,2], pc3=mbres$scores[nback_ind,3], 
                  pc4=mbres$scores[nback_ind,4], pc5=mbres$scores[nback_ind,5], labels=row.names(mbres$scores)[nback_ind])

df1$group <- factor(rep(1:3, each=60))
ggplot(subset(df1, group==1)) + geom_point(aes(pc1,pc2), color="red") + geom_text_repel(aes(pc1,pc2, label = labels)) + theme_classic(base_size = 10)


### train on video
Xl2 <- as.list(matrix_to_block_matrix(Xvideo, groups))
Y2 <- factor(c(zerostr(1:1092)))
mbres <- musubada(Y2, Xl2, ncomp=200, normalization="MFA", rank_k=200)

## complete nback
imp_fit <- softImpute(Xnback, rank.max=50, type="als", trace.it=TRUE)
Xcomplete <- complete(Xnback,imp_fit)
lens <- lapply(ndat, "[[", "len")
groups <- rep(1:length(lens), unlist(lens))
Xbm <- matrix_to_block_matrix(Xcomplete, groups)
Xl <- as.list(Xbm)

proj <- lapply(1:length(lens), function(i) {
  project(mbres, Xl[[i]], ncomp=200, table_index=i)[[1]]
})

Ynback <- factor(rep(as.character(ndat[[1]]$design$label), length(lens)))
Ynback2 <- factor(gsub("_[12]", "", as.character(Ynback)))

Ysubj <- factor(rep(1:length(lens), each=180))
Xcat <- do.call(rbind, proj)

preds <- unlist(lapply(levels(Ysubj), function(sid) {
  print(sids[as.integer(sid)])
  train_idx <- which(Ysubj != sid)
  test_idx <- which(Ysubj == sid)
  sda.r <- sda.ranking(Xcat[train_idx,], Ynback2[train_idx], ranking.score="avg", verbose=FALSE)
  keep <- which(sda.r[, "score"] > 2)
  keepid <- sda.r[, "idx"]
  keepid <- keepid[keep]
  
  sda.1 <- sda(Xcat[train_idx,keepid], Ynback2[train_idx], verbose=FALSE)
  ##rf.1 <- randomForest(Xcat[train_idx,keepid], Ynback2[train_idx], ntree=1500, mtry=8)
  p <- predict(sda.1, Xcat[test_idx,keepid], verbose=FALSE)
  print( sum(p$class == Ynback2[test_idx])/length(test_idx) * 100)
  p$class

}))

Ximagery <- do.call(cbind, lapply(idat, "[[", "mat"))
Yimagery <- idat[[1]]

Xnback <- do.call(cbind, lapply(ndat, "[[", "mat"))
#imp_fit <- softImpute(Xnback, rank.max=50, type="als", trace.it=TRUE)
imp_fit <- softImpute(Xall, rank.max=50, type="als", trace.it=TRUE)
Xcomplete <- complete(Xall,imp_fit)
lens <- lapply(ndat, "[[", "len")
groups <- rep(1:length(lens), unlist(lens))

Xbm <- matrix_to_block_matrix(Xcomplete, groups)
Xl <- as.list(Xbm)
Y <- factor(c(zerostr(1:1092), as.character(ndat[[1]]$design$label)))
nback_ind <- 1093:(1093+179)
mbres <- musubada(Y, Xl, ncomp=100, normalization="MFA", rank_k=200)


Xibm <- as.list(matrix_to_block_matrix(Ximagery, groups))
time <- idat[[1]]$design$time

pres_cos <- lapply(1:length(Xibm), function(i) {
  print(i)
  pred <- predict(mbres, Xibm[[i]], type="cosine", ncomp=100, table_index=i)
  pred <- pred[,nback_ind]
  lver <- idat[[i]]$design$LabelVersion
  react <- pred[cbind(1:length(lver), as.integer(lver))]
  
  dfx <- data.frame(time = idat[[i]]$design$time,
                    cond = ifelse(idat[[i]]$design$cresp == 1, "same", "different"),
                    #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                    acc = idat[[i]]$design$acc,
                    vivid = idat[[i]]$design$vivid,
                    react=react)
  
  dfx %>% group_by(time,cond) %>% dplyr::summarize(cor_acc=cor(react, acc, method="spearman",use="pairwise.complete.obs")) %>% mutate(sid=sids[i])
  #agg1 <- aggregate(react ~ time + vivid, FUN=mean)
  #agg1$sid <- sids[i]
  #agg1
})

pres_cos = do.call(rbind, pres_cos)




