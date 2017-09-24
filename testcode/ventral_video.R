#path <- "/Users/bbuchsbaum/analysis/hyper/ventral_surface/"
path <- "/Users/bbuchsbaum/analysis/hyper/all_surface/"
use_s5 <- FALSE
impute_nback <- TRUE

sids <- scan(paste0(path, "sids"), "")

zerostr <- function(vals) {
  ndigits <- log(max(vals), 10) + 1
  nzeros <- ndigits - as.integer(log(vals,10)) -1
  prefix <- sapply(nzeros, function(nz) paste(rep("0", times=nz), collapse=""))
  paste(prefix, vals, sep="")  
}


load_mat <- function(fname) {
  tmp <- readRDS(paste0(path, fname))
  d1 <- tmp$left$data
  d2 <- tmp$right$data
  d3 <- cbind(t(d1), t(d2))
  
  roi <- c(tmp$left$ROI, tmp$right$ROI)
  nodes <- c(tmp$left$nodes, tmp$right$nodes)
  
  list(mat=d3, design=tmp$design, roi=roi, nodes=nodes)
}

  
library(outliers)


fix_outliers <- function(mat, tscore=4.2) {
  
  nout <- 1
  iter <- 1
  while (nout > 0) {
    scmat <- abs(apply(mat, 2, outliers::scores))
    whout <- which(scmat > tscore, arr.ind=TRUE)
    nout <- nrow(whout)
    print(nout)
    medcol <- apply(mat,2,median)
    medfill <- medcol[whout[,2]]
  
    mat2 <- mat
    mat2[whout] <- medfill
    mat <- mat2
    
    
  }
  
  mat
}

vdat <- lapply(sids, function(sid) {
  
  if (use_s5) {
    v1 <- load_mat(paste0(sid, "_s5_video_1.RDS"))
    v2 <- load_mat(paste0(sid, "_s5_video_2.RDS"))
    v3 <- load_mat(paste0(sid, "_s5_video_3.RDS"))
  } else {
    v1 <- load_mat(paste0(sid, "_video_1.RDS"))
    v2 <- load_mat(paste0(sid, "_video_2.RDS"))
    v3 <- load_mat(paste0(sid, "_video_3.RDS"))
  }
  
  
  des <- rbind(v1$design, v2$design, v3$design)
  des$run <- factor(rep(1:3, each=364))
  
  mat <- rbind(v1$mat, v2$mat, v3$mat)
  
  keep <- apply(mat,2, function(x) sum(x==0)) == 0
  
  mat <- mat[,keep]
  mat <- fix_outliers(mat)
  
  roi <- v1$roi[keep]
  nodes <- v1$nodes[keep]
  
  rmat <- resid(lsfit(model.matrix(~ factor(des$run)), mat, intercept=FALSE))
  
  mat <- t(scale(t(rmat)))
  mat <- scale(mat)
  
  #ncomp <- fast_estim_ncomp(scale(mat), 1, 300)$bestcomp
  
  list(sid=sid, mat=mat, design=des, ncomp=150, keep=keep, roi=roi, nodes=nodes)
})

ndat <- lapply(vdat, function(vd) {
  sid <- vd$sid
  
  if (use_s5) {
    v1 <- load_mat(paste0(sid, "_s5_nback.RDS"))
  } else {
    v1 <- load_mat(paste0(sid, "_nback.RDS"))
  }
  
  v1$design$LabelVersion <- factor(paste0(as.character(v1$design$Label), "_", v1$design$Version))
  keep <- vd$keep
  mat <- v1$mat[, keep]
  mat <- fix_outliers(mat)
  
  rmat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
  
  mat <- group_means(v1$design$LabelVersion, rmat)
  mat <- t(scale(t(mat)))
  mat <- scale(mat)
  
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
  
  #ncomp <- fast_estim_ncomp(scale(mat), 1, 89)$bestcomp
  
  des <- data.frame(label=paste0(rep(labs, each=2), "_", rep(1:2, 90)), version=rep(1:2, 90))
  list(sid=sid, mat=mat_miss, design=des, ncomp=1, keep=keep, len=ncol(mat), roi=vd$roi, nodes=vd$nodes)
})


get_roi_matrix <- function(datlist, rnum) {
  res <- lapply(datlist, function(x) {
    m <- x$mat
    idx <- which(x$roi %in% rnum | (x$roi %in% (rnum + 1000)))
    m[,idx]
  })
}

get_roi_indices <- function(datlist, rnum) {
  lapply(datlist, function(x) {
    which(x$roi %in% rnum | (x$roi %in% (rnum + 1000)))
  })
  
}

library(softImpute)
rois <- sort(unique(vdat[[1]]$roi))


compute_rv <- function(x1, x2) {
  w1 <- is.na(x1[,1])
  w2 <- !is.na(x2[,1])
  coeffRV(x1[w1&w2,], x2[w2&w1,])$rvstd
}


compute_cor <- function(x1,x2) {
  cx1 <- cor(t(x1))
  cx2 <- cor(t(x2))
  
  ut1=cx1[upper.tri(cx1)]
  ut2=cx2[upper.tri(cx2)]
  
  keep <- !is.na(ut1) & !is.na(ut2)
  if (sum(keep) == 0) {
    NA
  } else {
    cor(ut1[keep], ut2[keep], method="spearman")
  }
}


compute_roi_importance <- function() {
  cres <- lapply(rois, function(rnum) {
    print(rnum)
    
    b2 <- block_matrix(get_roi_matrix(ndat, rnum))
    
    clist <- lapply(1:nblocks(b2), function(i) {
    
      unlist(lapply(1:nblocks(b2), function(j) {
        compute_cor(get_block(b2, i), get_block(b2, j))
      }))
    })
    
    cmat <- do.call(cbind, clist)
    diag(cmat) <- 0
    print(mean(cmat, na.rm=TRUE))
    cmat
  })
  
}

roi_cmats <- compute_roi_importance()
roi_imp <- sapply(roi_cmats, function(x) mean(x, na.rm=TRUE))
roi_imp[roi_imp < 0] <- 0

roi_keep <- which(roi_imp > .02)
best_rois <- rois[roi_keep]
roi_wts <- roi_imp[roi_keep]
roi_wts <- (roi_wts^2)/sum(roi_wts)

project_musu <- function(mres, bind, roinum, roi_weight) {
  nback_ind <- 1093:(1093+179)
  Xim <- lapply(1:length(idat), function(i) {
    idat[[i]]$mat[, bind[[i]]]
  })
  
  pres_cos <- lapply(1:length(Xim), function(i) {
    print(i)
    pred <- predict(mres, Xim[[i]], type="cosine", ncomp=mres$ncomp, table_index=i)
    pred <- pred[,nback_ind]
    lver <- idat[[i]]$design$LabelVersion
    react <- pred[cbind(1:length(lver), as.integer(lver))]
    
    dfx <- data.frame(time = idat[[i]]$design$time,
                      cond = ifelse(idat[[i]]$design$cresp == 1, "same", "different"),
                      #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                      acc = idat[[i]]$design$acc,
                      zvivid = scale(idat[[i]]$design$vivid),
                      zconf = scale(idat[[i]]$design$confidence),
                      react=react,
                      sid=sids[i],
                      weight=roi_weight,
                      roinum=roinum)
    
  })
  
  pres_cos = do.call(rbind, pres_cos)
}


get_musu_preds <- function(rnum) {
  print(rnum)
  b1 <- block_matrix(get_roi_matrix(vdat, rnum))
  b2 <- block_matrix(get_roi_matrix(ndat, rnum))
  
  b3 <- rbind(b1,b2)
  
  blens <- block_lengths(b3)
  bind <- get_roi_indices(ndat, rnum)
  
  rmax <- round(sqrt(median(blens)))
  message("rank: ", rmax)
  imp_fit <- softImpute(b3, rank.max=rmax, type="als", thresh=1e-04, trace.it=TRUE)
  Xcomplete <- complete(b3,imp_fit)
  
  groups <- rep(1:length(blens), unlist(blens))
  Xbm <- matrix_to_block_matrix(Xcomplete, groups)
  Xl <- as.list(Xbm)
  Y <- factor(c(zerostr(1:nrow(b1)), as.character(ndat[[1]]$design$label)))
  
  start <- nrow(b1) + 1
  nback_ind <- start:(start+nrow(b2)-1)
  
  mres <- musu_bada(Y, Xl, ncomp=rmax, normalization="None")
  ret <- project_musu(mres, bind, rnum, roi_wts[j])
  rm(mres)
  ret
  
}


musu_preds <- lapply(1:length(best_rois), function(j) {
    rnum <- best_rois[j]
    print(rnum)
    b1 <- block_matrix(get_roi_matrix(vdat, rnum))
    b2 <- block_matrix(get_roi_matrix(ndat, rnum))
    
    b3 <- rbind(b1,b2)
  
    blens <- block_lengths(b3)
    bind <- get_roi_indices(ndat, rnum)
    
    rmax <- round(sqrt(median(blens)))
    message("rank: ", rmax)
    imp_fit <- softImpute(b3, rank.max=rmax, type="als", thresh=1e-04, trace.it=TRUE)
    Xcomplete <- complete(b3,imp_fit)
    
    groups <- rep(1:length(blens), unlist(blens))
    Xbm <- matrix_to_block_matrix(Xcomplete, groups)
    Xl <- as.list(Xbm)
    Y <- factor(c(zerostr(1:nrow(b1)), as.character(ndat[[1]]$design$label)))
    
    start <- nrow(b1) + 1
    nback_ind <- start:(start+nrow(b2)-1)
    
    mres <- musu_bada(Y, Xl, ncomp=rmax, normalization="None")
    ret <- project_musu(mres, bind, rnum, roi_wts[j])
    rm(mres)
    ret
  
})




idat <- lapply(vdat, function(vd) {
  sid <- vd$sid
  if (use_s5) {
    v1 <- load_mat(paste0(sid, "_s5_imagery.RDS"))
  } else {
    v1 <- load_mat(paste0(sid, "_imagery.RDS"))
  }
  
  v1$design$LabelVersion <- factor(paste0(as.character(v1$design$label), "_", v1$design$version))
  keep <- vd$keep
  mat <- v1$mat[, keep]
 
  rmat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
  rmat <- fix_outliers(rmat)
  mat <- t(scale(t(rmat)))
  mat <- scale(mat)
  list(sid=sid, mat=mat, design=v1$design, keep=keep, len=ncol(mat), roi=vd$roi, nodes=vd$nodes)
  
})



### impute missing "version" from n-back task.
library(softImpute)



library(dplyr)

pres_cos <- do.call(rbind, musu_preds)
react_sum <- Reduce("+", lapply(musu_preds, function(p) p$react * p$weight))
pres_cos_sum <- musu_preds[[1]]
pres_cos_sum <- pres_cos_sum %>% mutate(react=react_sum)


reac_by_time <- pres_cos %>% group_by(time, roinum) %>% summarize(mean_react=mean(react))
reac_by_time_sum <- pres_cos_sum %>% group_by(time, sid) %>% summarise(mean_react=mean(react))
qplot(time, mean_react, data=reac_by_time_sum, geom=c("point", "line"))

reac_by_acc <- pres_cos %>% group_by(time, acc) %>% summarize(mean_react=mean(react))
reac_by_acc_sum <- pres_cos_sum %>% group_by(time, acc) %>% summarize(mean_react=mean(react))
qplot(time, mean_react, colour=factor(acc), data=reac_by_acc_sum, geom=c("point", "line"))


reac_by_viv <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% group_by(time, qvivid) %>% summarize(mean_react=mean(react))
reac_by_viv_sum <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% group_by(time, qvivid) %>% summarize(mean_react=mean(react))

reac_by_viv_cond <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% group_by(time, qvivid, cond) %>% summarize(mean_react=mean(react))
reac_by_viv_cond_sum <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% group_by(time, qvivid, cond) %>% summarize(wt_react_sum=sum(wt_react))


reac_by_viv_acc <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% group_by(time, qvivid, acc) %>% summarize(mean_react=mean(react))
reac_by_viv_acc_sum <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% group_by(time, qvivid, acc) %>% summarize(mean_react=mean(react))


reac_by_cond <- pres_cos %>%  group_by(time, cond) %>% summarize(mean_react=mean(react))
reac_by_cond_sum <- pres_cos_sum %>%  group_by(time, cond) %>% summarize(mean_react=mean(react))

reac_by_cond_acc <- pres_cos %>%  group_by(time, cond,acc) %>% summarize(mean_react=mean(react))
reac_by_cond_acc_sum <- pres_cos_sum %>%  group_by(time, cond, acc) %>% summarize(mean_react=mean(react))
qplot(time, mean_react, colour=factor(acc), facets = . ~ cond, data=reac_by_cond_acc_sum, geom=c("point", "line"))




reac_by_acc <- pres_cos %>%  group_by(time, acc) %>% summarize(mean_react=mean(react))

reac_by_acc_viv <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% filter(!is.na(qvivid)) %>% 
                                                                             group_by(time, acc, qvivid) %>% summarize(mean_react=mean(react))

reac_by_acc_viv_sum <- pres_cos %>% mutate(qvivid=ntile(zvivid, 3)) %>% filter(!is.na(qvivid)) %>% 
  group_by(time, acc, qvivid) %>% summarize(wt_react_sum=sum(wt_react))




reac_by_cond_conf <- pres_cos %>% mutate(qconf=ntile(zconf, 4)) %>% filter(!is.na(qconf)) %>% 
  group_by(time, cond, qconf) %>% summarize(mean_react=mean(react))

reac_by_conf <- pres_cos %>% mutate(qconf=ntile(zconf, 4)) %>% filter(!is.na(qconf)) %>% 
  group_by(time, qconf) %>% summarize(mean_react=mean(react))





################################################################
################################################################
################################################################
################################################################
################################################################
# Xnback <- do.call(cbind, lapply(ndat, "[[", "mat"))
# imp_fit <- softImpute(Xnback, rank.max=25, type="als", trace.it=TRUE)
# Xcomplete <- complete(Xnback,imp_fit)
# lens <- lapply(ndat, "[[", "len")
# groups <- rep(1:length(lens), unlist(lens))

#Xbm <- matrix_to_block_matrix(Xcomplete, groups)
#Xl <- as.list(Xbm)
#Y <- ndat[[1]]$design$label
#mbres <- musubada(Y, Xl, ncomp=20, normalization="MFA", rank_k=120)
#df1 <- data.frame(pc1=mbres$scores[,1], pc2=mbres$scores[,2], pc3=mbres$scores[,3], pc4=mbres$scores[,4], pc5=mbres$scores[,5], labels=row.names(mbres$scores))
#ggplot(df1) + geom_point(aes(pc1,pc2), color="red") + geom_text_repel(aes(pc1,pc2, label = labels)) + theme_classic(base_size = 10)



## musu_bada on full data set
mbres <- musu_bada(Y, Xl, ncomp=20, normalization="MFA", rank_k=200)
library(ggrepel)
df1 <- data.frame(pc1=mbres$scores[nback_ind,1], pc2=mbres$scores[nback_ind,2], pc3=mbres$scores[nback_ind,3], 
                  pc4=mbres$scores[nback_ind,4], pc5=mbres$scores[nback_ind,5], labels=row.names(mbres$scores)[nback_ind])

df1$group <- factor(rep(1:3, each=60))
ggplot(subset(df1, group==3), aes(pc1,pc2)) + 
  geom_text_repel(aes(pc1,pc2, label = labels)) + 
  theme_classic(base_size = 10)


### train on video
Xl2 <- as.list(matrix_to_block_matrix(Xvideo, groups))
Y2 <- factor(c(zerostr(1:1092)))
mbres <- musu_bada(Y2, Xl2, ncomp=200, normalization="MFA", rank_k=200)

## complete nback
imp_fit <- softImpute(Xnback, rank.max=50, type="als", trace.it=TRUE)
Xcomplete <- complete(Xnback,imp_fit)
lens <- lapply(ndat, "[[", "len")
groups <- rep(1:length(lens), unlist(lens))
Xbm <- matrix_to_block_matrix(Xcomplete, groups)
Xl <- as.list(Xbm)

## project nback data on to video space
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


## combine imagery data
Ximagery <- do.call(cbind, lapply(idat, "[[", "mat"))

## get imagery labels
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
mbres <- musu_bada(Y, Xl, ncomp=200, normalization="MFA", rank_k=200)
