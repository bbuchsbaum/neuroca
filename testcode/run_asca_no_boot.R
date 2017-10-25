library(Matrix)
library(tidyr)
sids <- seq(1001, 1018)


combine_mats <- function(l1, l2) {
  res <- lapply(1:length(l1), function(i) {
    m1 <- l1[[i]]$mat
    m2 <- l2[[i]]$mat
    m3 <- rbind(m1,m2)
    
    d1 <- l1[[i]]$design[, c("block", "trial", "syllable", "consonant", "vowel", "place", "sid", "speaker")]
    d2 <- l2[[i]]$design[, c("block", "trial", "syllable", "consonant", "vowel", "place")]
    
    d2$sid <- d1$sid[1]
    d1$task <- "speech"
    d2$task <- "mouthing"
    d2$speaker <- "internal"
    
    out <- rbind(d1,d2)
    out$task <- factor(out$task)
    out$block2 <- rep(rep(1:4, each=6), length.out=nrow(out))
    
    mred <- neuroca:::collapse(~ task + vowel + consonant + block2, m3, out)
    mred2 <- residualize(~ task, mred$X, mred$design)
    list(mat=mred2, design=mred$design, roi=l1[[i]]$roi)
    
  })
  
}

load_task_mat <- function(fname, sid, rnums="*") {
  tmp <- readRDS(fname)
  d1 <- tmp$lh$dat
  d2 <- tmp$rh$dat
  
  if (rnums[1] != "*") {
    keep_lh <- tmp$mouthing$lh$roi %in% rnums
    keep_rh <- tmp$mouthing$rh$roi %in% rnums
  } else {
    keep_lh <- 1:nrow(d1)
    keep_rh <- 1:nrow(d2)
  }
  
  d3 <- cbind(t(d1[keep_lh,]), t(d2[keep_rh,]))
  des <- tmp$design
  des$ftime <- factor(des$time)
  d3 <- t(scale(t(d3)))
  d3 <- neuropls:::residualize(~ftime, d3, des)
  
  
  des$sid <- sid
  list(mat=d3, design=des, roi=c(tmp$lh$roi, tmp$rh$roi))
  
}


load_mouthing_mat <- function(fname, sid, rnums="*") {
  tmp <- readRDS(fname)
  d1 <- tmp$mouthing$lh$dat
  d2 <- tmp$mouthing$rh$dat
  
  if (rnums[1] != "*") {
    keep_lh <- tmp$mouthing$lh$roi %in% rnums
    keep_rh <- tmp$mouthing$rh$roi %in% rnums
  } else {
    keep_lh <- 1:nrow(d1)
    keep_rh <- 1:nrow(d2)
  }
  
  d3 <- cbind(t(d1[keep_lh,]), t(d2[keep_rh,]))
  d3x <- resid(lsfit(model.matrix( ~ factor(tmp$mouthing$design$block)-1), d3, intercept=FALSE))
  
  d3 <- t(scale(t(d3x)))
  tmp$mouthing$design$sid <- sid
  des <- tmp$mouthing$design
  
  if (length(table(des$block)) == 4) {
    idx <- which(des$block == 4)
    med <- median(idx)
    des$block[idx[idx < med]] <- 4
    des$block[idx[idx >= med]] <- 5
  }
  
  list(mat=d3, design=des, roi=c(tmp$mouthing$lh$roi, tmp$mouthing$rh$roi))
}

load_speech_mat <- function(fname, sid, rnums="*") {
  tmp <- readRDS(fname)
  d1 <- tmp$speech$lh$dat
  d2 <- tmp$speech$rh$dat
  
  if (rnums[1] != "*") {
    keep_lh <- tmp$speech$lh$roi %in% rnums
    keep_rh <- tmp$speech$rh$roi %in% rnums
  } else {
    keep_lh <- 1:nrow(d1)
    keep_rh <- 1:nrow(d2)
  }
 
  d3 <- cbind(t(d1[keep_lh,]), t(d2[keep_rh,]))
  d3x <- resid(lsfit(model.matrix( ~ factor(tmp$speech$design$block)-1), d3, intercept=FALSE))
  d3 <- t(scale(t(d3x)))
  tmp$speech$design$sid <- sid
  list(mat=d3, design=tmp$speech$design, roi=c(tmp$speech$lh$roi, tmp$speech$rh$roi))
}

collapse_delay_vowel<- function(talist) {
  lapply(1:length(talist), function(i) {
    X <- talist[[i]]$mat
    des <- talist[[i]]$design
    des$half <- ifelse(des$run < 6, 1, 2)
    keep <- des$time %in% seq(0.0, 13.5, by=1.5)
    X0 <- X[keep,]
    des0 <- des[keep,]
    des0$run <- factor(des0$run)
    tmp <- neuropls:::collapse(~ cued.vowel + run, X0, des0)
    list(mat=tmp$X, design=tmp$design, roi=talist[[i]]$roi)
  })
}


fnames <- paste0("~/analysis/dual_aud/asca/", sids, "_betas_lang_sp_and_mo.RDS")
#fnames2 <- paste0("../", sids, "c/trial_data/betas_task.RDS")

splist <- lapply(1:length(sids), function(i) load_speech_mat(fnames[i], sids[i]))
molist <- lapply(1:length(sids), function(i) load_mouthing_mat(fnames[i], sids[i]))
bolist <- combine_mats(splist, molist)
#talist <- lapply(1:length(sids), function(i) load_task_mat(fnames2[i], sids[i]))

#talist_vowel <- collapse_delay_vowel(talist)


deslist_sp <- lapply(splist, "[[", "design")
deslist_mo <- lapply(molist, "[[", "design")
deslist_bo <- lapply(bolist, "[[", "design")
#deslist_ta <- lapply(talist, "[[", "design")

rnums <- sort(unique(splist[[1]]$roi))

run_asca <- function(dlist, rnum,nc=2, deslist, fac = "syllable") {
  
  Xret <- lapply(dlist, function(x) {
    idx <- x$roi %in% rnum
    m <- x$mat[,idx]
    list(mat=m, indices=which(idx))
  })
  
  Xlist <- lapply(Xret, "[[", "mat")
  folds <- lapply(deslist, "[[", "block")
  Yl <- lapply(deslist, "[[", fac)
  
  res <- mubada(Xlist, Yl, ncomp=nc, center=TRUE, scale=TRUE, normalization="MFA")
  attr(res, "indices") <- lapply(Xret, "[[", "indices")
  res
  
}

