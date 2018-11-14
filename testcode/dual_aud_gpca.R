library(neuroca)
library(fmrireg)
library(neighborweights)
library(dplyr)
library(tidyr)

dset <- readRDS("~/Dropbox/analysis/dual_aud/roi_test_data.rds")
dset$speech$design$speech <- "speech"
dset$speech$design$onset <- dset$speech$design$onset1/1000
sframe_speech <- sampling_frame(blocklens=rep(197,6), TR=1.5)
ev_speech <- event_model(onset ~ hrf(vowel, prefix="sp", basis="spmg1") + 
                           hrf(consonant, prefix="sp", basis="spmg1") + 
                           hrf(gender, prefix="sp", basis="spmg1"),
                           block= ~ block, durations=4, 
                         data=dset$speech$design, sampling_frame=sframe_speech)



dset$mouthing$design$mouthing <- "mouthing"
dset$mouthing$design$onset <- dset$mouthing$design$rep1/1000
sframe_mouthing <- sampling_frame(blocklens=rep(213,5), TR=1.5)

sframe_task <- sampling_frame(blocklens=rep(185,12), TR=1.5)
des_task <- dset$task$design %>% mutate(speech="speech", delay="delay", mouthing="mouthing",
                                        gender=cued.gender,
                                        onset1 = onset1/1000, onset2=onset2/1000, 
                                        delay.onset=delay.onset/1000, recall.onset=recall.onset/1000)

des_task <- des_task %>% mutate(vowel=cued.vowel, consonant=cued.consonant, syllable=cued.syll)
des_task <- gather(des_task, key="phase", value="onset", onset1, onset2, delay.onset, recall.onset) %>% arrange(block, onset)
des_task <- des_task %>% mutate(durations=case_when(phase=="onset1" ~ 1,
                                                    phase=="onset2" ~ 1,
                                                    phase=="delay.onset" ~ 8,
                                                    phase=="recall.onset" ~ 2))

otherind <- which( (des_task$phase == "onset1" & des_task$cued.order == 2) | (des_task$phase == "onset2" & des_task$cued.order == 1))
des_task$vowel[otherind] <- des_task$other.vowel[otherind]
des_task$consonant[otherind] <- des_task$other.cons[otherind]
des_task$syllable[otherind] <- des_task$other.syll[otherind]

ev_task <- event_model(onset ~ hrf(vowel, prefix="sp", subset=phase %in% c("onset1", "onset2")) + 
                               hrf(consonant, prefix="sp", subset=phase %in% c("onset1", "onset2")) + 
                               #hrf(speech, prefix="sp", subset=phase %in% c("onset1", "onset2")) +
                               hrf(vowel, prefix="mem", durations=8, subset=phase == "delay.onset") +
                               hrf(consonant, prefix="mem",durations=8, subset=phase == "delay.onset") +
                               #hrf(delay, durations=8, subset=phase == "delay.onset") +
                               #hrf(gender, durations=4, prefix="mem", subset=phase == "delay.onset") +
                               hrf(vowel, prefix="rec", durations=4, subset=phase == "recall.onset") +
                               hrf(consonant, prefix="rec", durations=4, subset=phase == "recall.onset"),
                               #hrf(mouthing, prefix="mg", durations=2, subset=phase == "recall.onset"),
                               #hrf(vowel, durations=durations) +
                               #hrf(consonant, durations=durations),
                      block= ~ block, 
                      data=des_task, sampling_frame=sframe_task)

ev_delay <- event_model(onset ~ hrf(syllable, prefix="sp", durations=8,
                                    subset=phase %in% c("delay.onset")),
                          #hrf(vowel, durations=durations) +
                          #hrf(consonant, durations=durations),
                          block= ~ block, 
                          data=des_task, sampling_frame=sframe_task)



do_plot <- function(comp1, comp2, idx) {
  sc <- res$u
  #sc <- scores(res)
  cd1 = as.vector(cor(sc[,comp1], dmat_all))
  cd2 = as.vector(cor(sc[,comp2], dmat_all))
  
  minx = min(cd1) - .4*sd(cd1)
  maxx = max(cd1) + .4*sd(cd1)
  miny = min(cd2) - .4*sd(cd2)
  maxy = max(cd2) + .4*sd(cd2)
  
  labs <- colnames(dmat)
  labs <- sapply(strsplit(labs, "\\["), "[[", 2)
  
  plot(cd1[idx], cd2[idx], type="n")
  text(cd1[idx], cd2[idx], labels=labs[idx])
}
  


run_model2 <- function(train_dset, test_dset, train_model, test_model,
                       k=10, k2=10, ncomp=5, sigma=1, exclude_zero=TRUE, data_knn=TRUE, ds=FALSE) {
  
  dmat <- rbind(design_matrix(train_model))
  rsums <- rowSums(dmat)
  
  bids <- train_model$sampling_frame$blockids
  
  if (exclude_zero) {
    tridx <- (1:length(bids))[rsums !=0]
  } else {
    tridx <- (1:length(bids))
  }
  
  dm <- as.matrix(dmat[tridx,])
  
  #Xmat <- rbind(dset[[1]]$data)
  Xmat <- train_dset$data
  Xtest <- test_dset$data
  
  Xtr <- Xmat[tridx,]
  
  if (sigma > 0) {
    nbd2 <- neighborweights::edge_weights(dm, k=k, 
                                          neighbor_mode="knn", weight_mode="heat", 
                                          sigma=sigma)
    if (ds) {
      Q <- neighborweights::make_doubly_stochastic(nbd2) 
    } else {
      Q <- nbd2
    }
    if (data_knn) {
      nbd_dat <- neighborweights::edge_weights(Xtr, k=k, 
                                               neighbor_mode="knn", weight_mode="normalized", 
                                               sigma=1)
      
      idx_overlap <- which(Q > 0 & nbd_dat > 0)
      idx_enemy <- which(Q == 0 & nbd_dat > 0)
      Q2 <- Q
      Q2[idx_overlap] <- Q[idx_overlap] * (1 + nbd_dat[idx_overlap])
      Q2[idx_enemy] <- nbd_dat[idx_enemy] * (1 - nbd_dat[idx_enemy])
      Q <- Q2
      diag(Q) <- 2
    }
    
   
    Q <- Q/RSpectra::eigs_sym(Q,k=1)$values
    res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=Q, ncomp=ncomp, center=TRUE, scale=TRUE)
  } else {
    res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=rep(1,nrow(Xtr)), ncomp=ncomp, center=TRUE)
  }
  
  dmtest <- as.matrix(design_matrix(test_model))
  proj <- project(res1, Xtest)
  p <- scorepred(proj, scores(res1), "cosine")
  pred <- do.call(rbind, lapply(1:nrow(p), function(i) {
    ord <- order(p[i,], decreasing=TRUE)[1:k2]
    p <- colMeans(dm[ord,])
  }))
  
  library(weights)
 
  cscore = sapply(1:nrow(pred), function(i) cor(pred[i,],dmtest[i,], method="kendall"))
  wts <- rowSums(dmtest)
  tres <- weights::wtd.t.test(cscore, weight=wts/sum(wts))
  list(tstat=tres$coefficients[1], cscore=cscore, wts=wts)
}

run_model_cv <- function(heldout, train_dset, train_model,
                       k=10, k2=10,ncomp=5, sigma=1, 
                       exclude_zero=TRUE, 
                       data_knn=TRUE,
                       ds=FALSE) {
  
  dmat <- rbind(design_matrix(train_model))
  rsums <- rowSums(dmat)
  
  bids <- train_model$sampling_frame$blockids
  
  if (exclude_zero) {
    tridx <- (1:length(bids))[rsums !=0 & bids != heldout]
  } else {
    tridx <- (1:length(bids))
  }
  
  testidx <- seq(1,length(bids))[bids == heldout]
  dm <- as.matrix(dmat[tridx,])
  
  #Xmat <- rbind(dset[[1]]$data)
  Xmat <- train_dset$data

  Xtr <- Xmat[tridx,]
  Xtest <- Xmat[testidx,]
  
  if (sigma > 0) {
    
    nbd2 <- neighborweights::edge_weights(dm, k=k, 
                                          neighbor_mode="knn", weight_mode="heat", 
                                          sigma=sigma)
    if (ds) {
      #Q <- neighborweights::make_doubly_stochastic(nbd2) 
      Q <- neighborweights:::normalize_adjacency(nbd2) 
    } else {
      Q <- nbd2
    }
   
    if (data_knn) {
      nbd_dat <- neighborweights::edge_weights(Xtr, k=k, 
                                            neighbor_mode="knn", weight_mode="normalized", 
                                            sigma=1)
      
      idx_overlap <- which(Q > 0 & nbd_dat > 0)
      idx_enemy <- which(Q == 0 & nbd_dat > 0)
      Q2 <- Q
      Q2[idx_overlap] <- Q[idx_overlap] * (1 + nbd_dat[idx_overlap])
      Q2[idx_enemy] <- nbd_dat[idx_enemy] * (1 - nbd_dat[idx_enemy])
      Q <- Q2
      #diag(Q) <- 2
    }
    
    diag(Q) <- 1/rowSums(Q)
    Q <- Q/RSpectra::eigs_sym(Q,k=1)$values
    
    res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=Q, ncomp=ncomp, center=TRUE, scale=TRUE)
  } else {
    res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=rep(1,nrow(Xtr)), ncomp=ncomp, center=TRUE)
  }
  

  dmtest <- as.matrix(dmat[testidx,])
  proj <- project(res1, Xtest)
  p <- scorepred(proj, scores(res1), "cosine")
  
  pred <- do.call(rbind, lapply(1:nrow(p), function(i) {
    ord <- order(p[i,], decreasing=TRUE)[1:k2]
    p <- colMeans(dm[ord,])
  }))
  
  library(weights)
  #browser()
  
  if (length(grep("basis\\[1\\]", colnames(dmtest))) > 0) {
    idx1 <- grep("basis\\[1\\]", colnames(dmtest))
    idx2 <- grep("basis\\[2\\]", colnames(dmtest))
    idx3 <- grep("basis\\[3\\]", colnames(dmtest))
    #browser()
  
    cscore1 = sapply(1:nrow(pred), function(i) cor(pred[i,idx1],dmtest[i,idx1], method="kendall"))
    cscore2 = sapply(1:nrow(pred), function(i) cor(pred[i,idx2],dmtest[i,idx2], method="kendall"))
    cscore3 = sapply(1:nrow(pred), function(i) cor(pred[i,idx3],dmtest[i,idx3], method="kendall"))
    wts <- rowSums(abs(dmtest))
    cscore <- (cscore1+cscore2+cscore3)/3
  } else {
    wts <- rowSums(abs(dmtest))
    cscore = sapply(1:nrow(pred), function(i) cor(pred[i,],dmtest[i,], method="kendall"))
  }
  
  tres <- weights::wtd.t.test(cscore, weight=wts/sum(wts))
  list(tstat=tres$coefficients[1], cscore=cscore, wts=wts)
}


run_models <- function(d1, d2, m1, m2) {
  grid <- expand.grid(k=seq(10,200,by=80), sigma=c(6), 
                      ncomp=c(4,8,12), data_knn=c(TRUE, FALSE))
  
  gres <- lapply(1:nrow(grid), function(i) {
    g <- grid[i,]
    res <- run_model2(d1,d2, m1, m2, k=g$k, k2=g$k, sigma=g$sigma, ncomp=g$ncomp, data_knn=g$data_knn)
    print(g)
    print(res$tstat)
    res$tstat
  })
  
  
  grid$score <- unlist(gres)
  grid
  
} 

run_models_cv <- function(d1, m1) {
  #grid <- expand.grid(k2=100, k=seq(100, 800, by=80), sigma=c(6), 
  #                    ncomp=c(5), ds=c(FALSE, TRUE))
  
  grid <- expand.grid(k=c(20,100,170), sigma=c(6), 
                      ncomp=c(7, 12,16), ds=c(FALSE), data_knn=c(TRUE, FALSE))
  
  gres <- lapply(1:nrow(grid), function(i) {
    g <- grid[i,]
    print(g)
    res <- lapply(unique(d1$design$block), function(j) {
      run_model_cv(j, d1,m1, k=g$k, k2=g$k, sigma=g$sigma, ncomp=g$ncomp, data_knn=g$data_knn, ds=g$ds)
    })
    cscore = unlist(lapply(res, "[[", "cscore"))
    wts = unlist(lapply(res, "[[", "wts"))
    #browser()
    ret <- sapply(res, "[[", "tstat")
    print(ret)
    weights::wtd.t.test(cscore, weight=wts)$coefficients[1]
   
  })
  
  
  grid$score <- unlist(gres)
  grid
  
} 

ev_speech_syll <- event_model(onset ~ hrf(syllable, prefix="sp",basis="spmg3"), 
                          block= ~ block, durations=1, 
                          data=dset$speech$design, sampling_frame=sframe_speech)

ev_task_speech_syll <- event_model(onset ~ hrf(syllable, prefix="sp", 
                                               subset=phase %in% c("onset1", "onset2")), 
                              block= ~ block, durations=1, 
                              data=des_task, sampling_frame=sframe_task)


res1 <- run_models(dset$speech, dset$task,ev_speech_syll,ev_task_speech_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()


ev_speech_vowel <- event_model(onset ~ hrf(vowel, prefix="sp", basis="spmg1"), 
                              block= ~ block, durations=4, 
                              data=dset$speech$design, sampling_frame=sframe_speech)

ev_task_speech_vowel <- event_model(onset ~ hrf(vowel, prefix="sp", 
                                                subset=phase %in% c("onset1", "onset2"),
                                                basis="spmg1"), 
                                   block= ~ block, durations=1, 
                                   data=des_task, sampling_frame=sframe_task)

ev_task_delay_vowel <- event_model(onset ~ hrf(vowel, prefix="sp", 
                                                durations=8,
                                                subset=phase %in% c("delay.onset"),
                                                basis="spmg1"), 
                                    block= ~ block, 
                                    data=des_task, sampling_frame=sframe_task)

ev_task_delay_syll <- event_model(onset ~ hrf(syllable, prefix="sp", 
                                               durations=8,
                                               subset=phase %in% c("delay.onset"),
                                               basis="spmg1"), 
                                   block= ~ block, 
                                   data=des_task, sampling_frame=sframe_task)

res1 <- run_models(dset$speech, dset$task,ev_speech_vowel,ev_task_speech_vowel)
qplot(k2, score, colour=factor(k), facets = . ~ factor(ncomp), data=res1) + geom_line()


res1 <- run_models(dset$speech, dset$task,ev_speech_vowel,ev_task_delay_vowel)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()




ev_speech_cons <- event_model(onset ~ hrf(consonant, prefix="sp", basis="spmg1"), 
                               block= ~ block, durations=4, 
                               data=dset$speech$design, sampling_frame=sframe_speech)

ev_speech_speaker <- event_model(onset ~ hrf(speaker, prefix="sp", basis="spmg1"), 
                              block= ~ block, durations=4, 
                              data=dset$speech$design, sampling_frame=sframe_speech)

des_task <- des_task %>% mutate(speaker=factor(paste0("S", cued.speaker)))
ev_task_delay_speaker <- event_model(onset ~ hrf(speaker, prefix="sp",
                                                 subset=phase=="delay.onset"), 
                                   block= ~ block, durations=1, 
                                   data=des_task, sampling_frame=sframe_task)


ev_task_speech_cons <- event_model(onset ~ hrf(consonant, prefix="sp", basis="spmg1"), 
                                    block= ~ block, durations=1, 
                                    data=des_task, sampling_frame=sframe_task)

ev_task_rec_syll <- event_model(onset ~ hrf(syllable, prefix="sp", basis="spmg1",
                                            subset=phase == "recall.onset"), 
                                   block= ~ block, durations=1, 
                                   data=des_task, sampling_frame=sframe_task)


res1 <- run_models(dset$speech, dset$task,ev_speech_cons,ev_task_speech_cons)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

ev_mouthing_syll <- event_model(onset ~ hrf(syllable, prefix="sp", basis="spmg3"), 
                              block= ~ block, durations=4, 
                              data=dset$mouthing$design, sampling_frame=sframe_mouthing)

ev_mouthing_vowel <- event_model(onset ~ hrf(vowel, prefix="sp", basis="spmg1"), 
                                block= ~ block, durations=4, 
                                data=dset$mouthing$design, sampling_frame=sframe_mouthing)

ev_mouthing_cons <- event_model(onset ~ hrf(consonant, prefix="sp", basis="spmg1"), 
                                 block= ~ block, durations=4, 
                                 data=dset$mouthing$design, sampling_frame=sframe_mouthing)

res1 <- run_models(dset$speech, dset$mouthing,ev_speech_syll,ev_mouthing_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()


res1 <- run_models(dset$speech, dset$mouthing,ev_speech_vowel,ev_mouthing_vowel)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$speech, dset$task,ev_mouthing_vowel,ev_task_delay_vowel)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$mouthing, dset$task,ev_mouthing_syll,ev_task_delay_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$mouthing, dset$task,ev_mouthing_syll,ev_task_speech_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$speech, dset$task,ev_speech_syll,ev_task_speech_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$speech, dset$mouthing,ev_speech_cons,ev_mouthing_cons)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$mouthing, dset$speech,ev_mouthing_cons,ev_speech_cons)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$mouthing, dset$task,ev_mouthing_syll,ev_task_rec_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$speech, dset$task,ev_speech_syll,ev_task_rec_syll)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()

res1 <- run_models(dset$speech, dset$task,ev_speech_speaker,ev_task_delay_speaker)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=res1) + geom_line()


run_model3 <- function(heldout=1, k=10, ncomp=5, sigma=1, ds=FALSE) {
  print(heldout)
  dmat_syll <- rbind(design_matrix(ev_delay))
  rsums <- rowSums(dmat_syll) 
  
  bids <- sframe_task$blockids
  tridx <- which(bids!=heldout)
  testidx <- which(bids==heldout)
  
  dm <- as.matrix(dmat_syll[tridx,])
  
  Xmat <- rbind(dset[[3]]$data)
  Xtr <- Xmat[tridx,]
  Xtest <- Xmat[testidx,]
  if (sigma > 0) {
    nbd2 <- neighborweights::edge_weights(dm, k=ncomp, 
                                          neighbor_mode="knn", weight_mode="heat", 
                                          sigma=sigma)
    #D <- Diagonal(x=1/sqrt(rowSums(nbd2)))
    #Q <- nbd2 %*% D %*% nbd2
    
    if (ds) {
      Q <- neighborweights::make_doubly_stochastic(nbd2) 
    } else {
      Q <- nbd2
    }
    Q <- Q/RSpectra::eigs_sym(Q,k=1)$values
    #res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=rep(1,nrow(Xtr)), ncomp=10, center=TRUE)
    res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=Q, ncomp=ncomp, center=TRUE)
  } else {
    res1 <- genpca(Xtr, A=rep(1, ncol(Xtr)), M=rep(1,nrow(Xtr)), ncomp=ncomp, center=TRUE)
  }
  
  
  dmtest <- as.matrix(design_matrix(ev_delay))[testidx,]
  proj <- project(res1, Xtest)
  p <- scorepred(proj, scores(res1), "cosine")
  pred <- do.call(rbind, lapply(1:nrow(p), function(i) {
    ord <- order(p[i,], decreasing=TRUE)[1:k]
    p <- colMeans(dm[ord,])
  }))
  
  library(weights)
  cscore = sapply(1:nrow(pred), function(i) cor(pred[i,],dmtest[i,]))
  wts <- rowSums(dmtest)
  tres <- weights::wtd.t.test(cscore, weight=wts/sum(wts))
  tres$coefficients[1]
}


grid <- expand.grid(k=seq(10, 100, by=20), sigma=c(4,8,16), ncomp=c(4,8,12,20))

gres <- lapply(1:nrow(grid), function(i) {
  g <- grid[i,]
  ret <- sapply(1:12, function(j) run_model3(j,k=g$k, sigma=g$sigma, ncomp=g$ncomp))
  print(g)
  #print(ret)
  print(mean(ret))
  mean(ret)
})


grid$score <- unlist(gres)
library(ggplot2)
qplot(ncomp, score, colour=factor(k), facets = . ~ factor(sigma), data=grid)

