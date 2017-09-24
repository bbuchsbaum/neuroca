#path <- "/Users/bbuchsbaum/analysis/hyper/ventral_surface/"
args <- commandArgs(trailingOnly=TRUE)
ROI <- args[[1]]
use_s5 <- as.logical(args[[2]])

path <- "/Users/bbuchsbaum/analysis/hyper/all_surface/"
#use_s5 <- TRUE
#impute_nback <- TRUE
use_video <- TRUE


library(softImpute)
library(lme4)
library(broom)


#library(sjPlot)
#library(sjmisc)


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


load_vdat <- function() {
  lapply(sids, function(sid) {
  
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
    
    roi <- v1$roi[keep]
    nodes <- v1$nodes[keep]
  
    rmat <- resid(lsfit(model.matrix(~ factor(des$run)), mat, intercept=FALSE))
    rmat <- fix_outliers(rmat)
    
    sdavg <- mean(apply(rmat,1,sd))
    mat <- rmat/sdavg
    mat <- scale(mat, scale=FALSE)
  
    #ncomp <- fast_estim_ncomp(scale(mat), 1, 300)$bestcomp
  
    list(sid=sid, mat=mat, design=des, ncomp=150, keep=keep, roi=roi, nodes=nodes)
  })
}



load_ndat <- function(vdat=NULL) {
  lapply(1:length(sids), function(i) {
    sid <- sids[i]

    if (use_s5) {
      v1 <- load_mat(paste0(sid, "_s5_nback.RDS"))
    } else {
      v1 <- load_mat(paste0(sid, "_nback.RDS"))
    }
  
    v1$design$LabelVersion <- factor(paste0(as.character(v1$design$Label), "_", v1$design$Version))
    mat <- v1$mat
    
    if (use_video) {
      vd <- vdat[[i]]
      keep <- vd$keep
      roi <- vd$roi
      nodes <- vd$nodes
    }
    else {
      keep <- apply(mat,2, function(x) sum(x==0)) == 0
      roi <- v1$roi[keep]
      nodes <- v1$nodes[keep]
    }
  
    mat <- mat[, keep]
    
    rmat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
    mat <- fix_outliers(rmat)
    
    mat <- group_means(v1$design$LabelVersion, mat)
    sdavg <- mean(apply(mat,1,sd))
    mat <- mat/sdavg
    mat <- scale(mat, scale=FALSE)
   
  
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
    
    des <- data.frame(label=paste0(rep(labs, each=2), "_", rep(1:2, 90)), version=rep(1:2, 90))
    list(sid=sid, mat=mat_miss, design=des, ncomp=1, keep=keep, len=ncol(mat), roi=roi, nodes=nodes)
  })
}

load_idat <- function(rnums) {
  lapply(ndat, function(vd) {
    sid <- vd$sid
    print(sid)
    if (use_s5) {
      v1 <- load_mat(paste0(sid, "_s5_imagery.RDS"))
    } else {
      v1 <- load_mat(paste0(sid, "_imagery.RDS"))
    }
    
    v1$design$LabelVersion <- factor(paste0(as.character(v1$design$label), "_", v1$design$version))
    keep <- vd$keep 
    mat <- v1$mat[, keep]
    
    mat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
    mat <- fix_outliers(mat)
    
    spmat <- split_matrix(mat, ordered(v1$design$time))
    
    rois <- vd$roi
    keep2 <- rois %in% c(rnums, rnums+1000)
    
    spmat <- lapply(spmat, function(m) {
      sdavg <- mean(apply(m,1,sd))
      m <- m/sdavg
      m <- scale(m, scale=FALSE)
      m[,keep2]
     
    })
    
    
    list(sid=sid, spmat=spmat, design=v1$design, keep=keep, len=ncol(mat), roi=rois[keep2], nodes=vd$nodes[keep2])
    
  })
}




if (use_video) {
  vdat <- load_vdat()
  ndat <- load_ndat(vdat)
} else {
  ndat <- load_ndat()
}


  
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


#rois <- sort(unique(ndat[[1]]$roi))


# compute_rv <- function(x1, x2) {
#   w1 <- is.na(x1[,1])
#   w2 <- !is.na(x2[,1])
#   coeffRV(x1[w1&w2,], x2[w2&w1,])$rvstd
# }
# 


compute_avg_cor <- function(x, xl) {
  cl1 <- lapply(xl, function(x1) cor(t(x1)))
  cl2 <- cor(t(x))
  
  cavg <- matrix(rowMeans(do.call(cbind, lapply(cl1, as.vector)), na.rm=TRUE), 180,180)
   
  ut1=cavg[upper.tri(cavg)]
  ut2=cl2[upper.tri(cl2)]
  
  keep <- !is.na(ut1) & !is.na(ut2)
  cor(ut1[keep], ut2[keep], method="spearman")
 
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
    acor <- unlist(lapply(1:nblocks(b2), function(i) {
      
      x <- get_block(b2, i)
      
      oi <- seq(1, nblocks(b2))[-i]
      xl <- lapply(oi, function(i) get_block(b2, i))
      
      compute_avg_cor(x,xl)
    }))
    
    ret <- mean(acor)
    print(ret)
    ret
    
  })
  
}

#roi_imp <- unlist(compute_roi_importance())



project_fscores <- function(mres, idat, tind, scoremat, labs) {
  tvals <- sort(unique(idat[[1]]$design$time))
  
  Xim <- lapply(idat, function(x) x$spmat[[tind]])
  
  pres_cos <- lapply(1:length(Xim), function(i) {
    print(i)
    pred <- predict(mres, Xim[[i]], type="scores", ncomp=mres$ncomp, table_index=i)
    
    sim <- proxy::simil(pred, scoremat, method="cosine")
    
    des <- subset(idat[[i]]$design, time==tvals[tind])
    lver <- des$LabelVersion
    ind <- match(as.character(lver), labs)
    
    react <- sim[cbind(1:length(lver), ind)]
    
    other_version <- ifelse(des$version == 1, 2, 1)
    other_label <- paste0(des$label, "_", other_version)
    
    ind_other <- match(as.character(other_label), labs)
    
    react_other <- sim[cbind(1:length(lver), ind_other)]
    
    dfx1 <- data.frame(time = tvals[tind],
                       cond = ifelse(des$cresp == 1, "same", "different"),
                       #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                       acc = des$acc,
                       zvivid = scale(des$vivid),
                       zconf = scale(des$confidence),
                       react=react,
                       sid=sids[i],
                       label=des$label,
                       label_version=des$LabelVersion,
                       type="react")
    
    dfx2 <- data.frame(time = tvals[tind],
                       cond = ifelse(des$cresp == 1, "same", "different"),
                       #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                       acc = des$acc,
                       zvivid = scale(des$vivid),
                       zconf = scale(des$confidence),
                       react=react_other,
                       sid=sids[i],
                       label=des$label,
                       label_version=des$LabelVersion,
                       type="react_other")
    
    dfx3 <- data.frame(time = tvals[tind],
                       cond = ifelse(des$cresp == 1, "same", "different"),
                       #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                       acc = des$acc,
                       zvivid = scale(des$vivid),
                       zconf = scale(des$confidence),
                       react=react-react_other,
                       sid=sids[i],
                       label=des$label,
                       label_version=des$LabelVersion,
                       type="react_diff")
    
    rbind(dfx1,dfx2,dfx3)
  })
  
  do.call(rbind, pres_cos)
}

project_musu <- function(mres, idat, tind) {
  nback_ind <- 1093:(1093+179)
  
  tvals <- sort(unique(idat[[1]]$design$time))
  
  Xim <- lapply(idat, function(x) x$spmat[[tind]])
  
  pres_cos <- lapply(1:length(Xim), function(i) {
    print(i)
    pred <- predict(mres, Xim[[i]], type="cosine", ncomp=mres$ncomp, table_index=i)
    #browser()
    if (use_video) {
      pred <- pred[,nback_ind]
    }
    
    des <- subset(idat[[i]]$design, time==tvals[tind])
    lver <- des$LabelVersion
    
    if (use_video) {
      ind <- match(as.character(lver), mres$Y[[1]][nback_ind])
    } else {
      ind <- match(as.character(lver), mres$Y[[1]])
    }
    
    react <- pred[cbind(1:length(lver), ind)]
    
    other_version <- ifelse(des$version == 1, 2, 1)
    other_label <- paste0(des$label, "_", other_version)
    
    if (use_video) {
      ind_other <- match(as.character(other_label), mres$Y[[1]][nback_ind])
    } else {
      ind_other <- match(as.character(other_label), mres$Y[[1]])
    }
    
   
    react_other <- pred[cbind(1:length(lver), ind_other)]
   
    
    dfx1 <- data.frame(time = tvals[tind],
                      cond = ifelse(des$cresp == 1, "same", "different"),
                      #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                      acc = des$acc,
                      zvivid = scale(des$vivid),
                      zconf = scale(des$confidence),
                      react=react,
                      sid=sids[i],
                      label=des$label,
                      label_version=des$LabelVersion,
                      type="react")
    
    dfx2 <- data.frame(time = tvals[tind],
                       cond = ifelse(des$cresp == 1, "same", "different"),
                       #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                       acc = des$acc,
                       zvivid = scale(des$vivid),
                       zconf = scale(des$confidence),
                       react=react_other,
                       sid=sids[i],
                       label=des$label,
                       label_version=des$LabelVersion,
                       type="react_other")
    
    dfx3 <- data.frame(time = tvals[tind],
                       cond = ifelse(des$cresp == 1, "same", "different"),
                       #acc = ifelse(idat[[i]]$design$acc == 1, "correct", "incorrect"),
                       acc = des$acc,
                       zvivid = scale(des$vivid),
                       zconf = scale(des$confidence),
                       react=react-react_other,
                       sid=sids[i],
                       label=des$label,
                       label_version=des$LabelVersion,
                       type="react_diff")
    
    rbind(dfx1,dfx2,dfx3)
    
  })
  
  do.call(rbind, pres_cos)
}


get_nback_fscores_video <- function(vdat, ndat, rnum) {
  b1 <- block_matrix(get_roi_matrix(vdat, rnum))
  blens <- block_lengths(b1)
  bind <- get_roi_indices(vdat, rnum)
  
  rmax <- round(min(sqrt(median(blens)), sqrt(nrow(b1))))
  #rmax <- 200
  message("rank: ", rmax)
  
  groups <- rep(1:length(blens), unlist(blens))
  Xl <- as.list(b1)
  Y <- factor(c(zerostr(1:nrow(b1))))
  mres <- musu_bada(Y, Xl, ncomp=rmax, normalization="None")
  
  fset <- lapply(1:length(ndat), function(i) {
    valid <- which(!is.na(ndat[[i]]$mat[,1]))
    fscores <- predict(mres, ndat[[i]]$mat[valid,bind[[i]]], type="scores", table_index=i)
    fscores <- fscores/sd(fscores[,1])
    out <- matrix(NA, nrow(ndat[[i]]$mat), ncol(fscores))
    out[valid,] <- fscores
    list(id=i, valid=valid, fscores=out)
  })

  fscores <- abind(lapply(fset, "[[", "fscores"), along=3)
  scoremat <- apply(fscores, c(1,2), mean, na.rm=TRUE)
  
  labs <- as.character(ndat[[1]]$design$label)
               
  pres <- do.call(rbind, lapply(1:16, function(ti) {
    print(ti)
    project_fscores(mres, idat, ti, scoremat, labs=labs)
  }))
}

get_musu_preds <- function(rnum) {
  print(rnum)
  
  if (use_video) {
    b1 <- block_matrix(get_roi_matrix(vdat, rnum))
    b2 <- block_matrix(get_roi_matrix(ndat, rnum))
  
    b3 <- rbind(b1,b2)
  } else {
    b3 <- block_matrix(get_roi_matrix(ndat, rnum))
  }
  
  blens <- block_lengths(b3)
  bind <- get_roi_indices(ndat, rnum)
  
  rmax <- round(min(sqrt(median(blens)), sqrt(nrow(b3))))
  #rmax <- 200
  message("rank: ", rmax)
  
  imp_fit <- softImpute(b3, rank.max=rmax, type="als", thresh=1e-04, trace.it=TRUE)
  Xcomplete <- complete(b3,imp_fit)
  
  #pcors <- sapply(1:90, function(i) {
  #   cor(Xcomplete[i*2-1,], Xcomplete[i*2,])
  # })
    
   
  groups <- rep(1:length(blens), unlist(blens))
  Xbm <- matrix_to_block_matrix(Xcomplete, groups)
  
  Xl <- as.list(Xbm)
  
  if (use_video) {
    Y <- factor(c(zerostr(1:nrow(b1)), as.character(ndat[[1]]$design$label)))
    start <- nrow(b1) + 1
    nback_ind <- start:(start+nrow(b2)-1)
  } else {
    Y <- factor(as.character(ndat[[1]]$design$label))
    start <- 1
    nback_ind <- 1:nrow(b3)
  }
  
  
  mres <- musu_bada(Y, Xl, ncomp=rmax, normalization="None")
  
  pres <- do.call(rbind, lapply(1:16, function(ti) {
    print(ti)
    project_musu(mres, idat, ti)
  }))
  
}

#best_rois <- rois[which(roi_imp > roi_thresh)]
#best_rois <- c(11125,11126,11127,11130)
idat <- load_idat(ROI)
#pres_cos <- get_musu_preds(ROI)
pres_cos <- get_nback_fscores_video(vdat, ndat, ROI)
library(fmrireg)
imagery_hrf <- hrf_gaussian %>% gen_hrf_lagged(0) %>% gen_hrf_blocked(6) %>% HRF(name="hrf_imagery") 
probe_hrf <- hrf_gaussian %>% gen_hrf_lagged(10) %>% gen_hrf_blocked(2.3) %>% HRF(name="hrf_probe") 

pres_cos$hrf_imag <- imagery_hrf(pres_cos$time)
pres_cos$hrf_probe <- probe_hrf(pres_cos$time)
pres_cos$ROI <- ROI
# library(optimx)
# 
lres1 <- lmer(react ~ hrf_imag + hrf_probe + (1 | sid) + (1 | label) + (0 + hrf_imag | sid) +
                 (0 + hrf_probe | sid) + (0 + hrf_imag | label) + (0 + hrf_probe | label), data=pres_cos, subset=type=="react")
 
 
lres1a <- lmer(react ~ hrf_imag + hrf_probe + (1 | sid) + (1 | label) + (0 + hrf_imag | sid) +
                 (0 + hrf_probe | sid) + (0 + hrf_imag | label) + (0 + hrf_probe | label), data=pres_cos, subset=type=="react_diff")
 
lres2 <- lmer(react ~ hrf_imag + hrf_probe + cond + hrf_probe:cond + (1 | sid) + (0 + hrf_imag | sid) + (0 + hrf_probe | sid) + 
                (0 + hrf_imag:cond | sid), data=pres_cos, subset=type=="react")
               
lres3 <- lmer(react ~ hrf_imag + hrf_probe + hrf_imag:zvivid  + hrf_probe:zvivid + (1 | sid) + (0 + hrf_imag | sid) + (0 + hrf_probe | sid) + 
                 (0 + hrf_imag:zvivid | sid) + (0 + hrf_probe:zvivid | sid), data=pres_cos, subset=type=="react")

lres3a <- lmer(react ~ hrf_imag + hrf_probe + hrf_imag:zvivid  + hrf_probe:zvivid + (1 | sid) + (0 + hrf_imag | sid) + (0 + hrf_probe | sid) + 
                 (0 + hrf_imag:zvivid | sid) + (0 + hrf_probe:zvivid | sid), data=pres_cos, subset=type=="react")

lres3b <- lmer(react ~ hrf_imag + hrf_probe + hrf_imag:zvivid  + hrf_probe:zvivid + (1 | sid) + (0 + hrf_imag | sid) + (0 + hrf_probe | sid) + 
                 (0 + hrf_imag:zvivid | sid) + (0 + hrf_probe:zvivid | sid), data=subset(pres_cos,acc==1), subset=type=="react")
 
lres4 <- lmer(react ~ hrf_imag + hrf_probe + hrf_imag:acc  + hrf_probe:acc  + (1 | sid) + (0 + hrf_imag | sid) + (0 + hrf_probe | sid) + 
                 (0 + hrf_imag:acc | sid) + (0 + hrf_probe:acc | sid), data=pres_cos, subset=type=="react")


out <- list(lres1, lres1a, lres2, lres3, lres3a, lres3b, lres4)
for (i in 1:length(out)) {
  fit <- out[[i]]
  
  tout <- tidy(fit)
  write.table(tout, file=paste0("ROI_", ROI, "lmfit_", i, ".txt"))
}

write.table(pres_cos, file=paste0("ROI_", ROI, "_react.txt"))





# 
# 
# ### impute missing "version" from n-back task.
# 
# library(dplyr)
# 
# 
# reac_by_time <- pres_cos %>% group_by(time, type) %>% summarize(react=mean(react, na.rm=TRUE))
#                                                        
# qplot(time, react, data=reac_by_time, geom=c("point", "line")) + facet_wrap( ~ type)
# 
# reac_by_acc <- pres_cos %>% group_by(time, acc, type) %>% summarize(react=mean(react, na.rm=TRUE))
#                                                               
# qplot(time, react, colour=factor(acc), data=reac_by_acc, geom=c("point", "line")) + facet_wrap( ~ type)
# 
# 
# reac_by_viv <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3)) %>% ungroup() %>% group_by(time, qvivid, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
#   
# qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv, !is.na(qvivid)), geom=c("point", "line")) + facet_wrap( ~ type)
# 
# 
# reac_by_viv_cond <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3)) %>% ungroup() %>% group_by(time, qvivid, cond, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
# 
# qplot(time, react, colour=factor(qvivid), data=subset(reac_by_viv_cond, !is.na(qvivid)), geom=c("point", "line")) + facet_wrap(cond ~ type)
# 
# 
# reac_by_viv_acc <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3)) %>% ungroup() %>% group_by(time, qvivid, acc, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
# 
# qplot(time, react, colour=factor(acc), data=subset(reac_by_viv_acc, !is.na(qvivid)), geom=c("point", "line")) + facet_wrap(qvivid ~ type)
# 
# 
# reac_by_cond <- pres_cos %>%  group_by(time, cond, type) %>% summarize(react=mean(react, na.rm=TRUE)) 
#                                                               
# qplot(time, react, colour=factor(cond), data=reac_by_cond, geom=c("point", "line")) + facet_wrap(~ type)
# 
# reac_by_cond_acc <- pres_cos %>%  group_by(time, cond,acc,type) %>% summarize(react=mean(react, na.rm=TRUE))
#                                                                       
# 
# qplot(time, react, colour=factor(acc), facets = . ~ cond, data=reac_by_cond_acc, geom=c("point", "line")) + facet_wrap(cond ~ type)
# 
# reac_by_conf <- pres_cos %>%  group_by(sid) %>% mutate(qconf=ntile(zconf, 3)) %>% ungroup() %>% group_by(time, qconf, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
#   
# qplot(time, react, colour=factor(qconf), data=subset(reac_by_conf, !is.na(qconf)), geom=c("point", "line")) + facet_wrap(~ type)
# 
# reac_by_conf_cond <- pres_cos %>% group_by(sid) %>% mutate(qconf=ntile(zconf, 3)) %>% ungroup()  %>% group_by(time, qconf, cond, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
# 
# qplot(time, react, colour=factor(qconf), facets = . ~ cond, data=subset(reac_by_conf_cond, !is.na(qconf)), geom=c("point", "line")) +
#   facet_wrap(cond ~ type)
# 
# reac_by_conf_acc <- pres_cos %>%  group_by(sid) %>% mutate(qconf=ntile(zconf, 3)) %>% ungroup() %>% group_by(time, qconf, acc, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
# 
# qplot(time, react, colour=factor(acc), data=subset(reac_by_conf_acc, !is.na(qconf)), geom=c("point", "line")) +
#   facet_wrap(qconf ~ type)
# 
# reac_by_viv_conf <- pres_cos %>% group_by(sid) %>% mutate(qvivid=ntile(zvivid, 3), qconf=ntile(zconf, 3)) %>% ungroup() %>% 
#   group_by(time, qvivid, qconf, type) %>% 
#   summarize(react=mean(react, na.rm=TRUE))
# 
# qplot(time, react, colour=factor(qconf), data=subset(reac_by_viv_conf, !is.na(qvivid) & !is.na(qconf) & type=="react"), geom=c("point", "line")) + 
#   facet_wrap(~ qvivid)

load_vdat_herve <- function(roi_keep) {
  lapply(sids, function(sid) {
    v1 <- load_mat(paste0(sid, "_video_1.RDS"))
    v2 <- load_mat(paste0(sid, "_video_2.RDS"))
    v3 <- load_mat(paste0(sid, "_video_3.RDS"))
    
    
    des <- rbind(v1$design, v2$design, v3$design)
    des$run <- factor(rep(1:3, each=364))
    
    mat <- rbind(v1$mat, v2$mat, v3$mat)
    
    keep <- apply(mat,2, function(x) sum(x==0)) == 0
    
    mat <- mat[,keep]
    
    roi <- v1$roi[keep]
    nodes <- v1$nodes[keep]
    
    rmat <- resid(lsfit(model.matrix(~ factor(des$run)), mat, intercept=FALSE))
    rmat <- fix_outliers(rmat)
    
    keep_roi <- roi %in% roi_keep
    samples <- as.integer(seq(15, nrow(rmat), length.out=100))
    
    list(sid=sid, mat=rmat[samples,keep_roi], design=des[samples,], roi=roi[keep_roi], nodes=nodes[keep_roi], keep=keep, keep_roi=keep_roi)
  })
}

load_ndat_herve <- function(vdat2) {
  ret <- lapply(1:length(sids), function(i) {
    sid <- sids[i]
    v1 <- load_mat(paste0(sid, "_nback.RDS"))
    
    v1$design$LabelVersion <- factor(paste0(as.character(v1$design$Label), "_", v1$design$Version))
    
    mat <- v1$mat
    
    vd <- vdat2[[i]]
    keep <- vd$keep
    roi <- vd$roi
    nodes <- vd$nodes
    
    mat <- mat[, keep]
    rmat <- resid(lsfit(model.matrix(~ factor(v1$design$run)), mat, intercept=FALSE))
    mat <- fix_outliers(rmat)
    
    mat <- group_means(v1$design$LabelVersion, mat)
    mat <- mat[,vd$keep_roi]
    
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
    
    row.names(mat_miss) <- paste0(rep(labs, each=2), "_", rep(1:2, 90))
    mat_miss
  })
  
  names(ret) <- sids
  ret
  
}


