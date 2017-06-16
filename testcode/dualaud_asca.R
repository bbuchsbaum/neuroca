library(Matrix)
path <- "/Users/bbuchsbaum/analysis/dual_aud/asca/"

sids <- seq(1001, 1018)


combine_mats <- function(l1, l2) {
  lapply(1:length(l1), function(i) {
    m1 <- l1[[i]]$mat
    m2 <- l2[[i]]$mat
    m3 <- rbind(m1,m2)
    
    d1 <- l1[[i]]$design[, c("block", "trial", "syllable", "consonant", "vowel", "place", "task", "sid")]
    d2 <- l2[[i]]$design[, c("block", "trial", "syllable", "consonant", "vowel", "place", "task", "sid")]
  })
  
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
  d3 <- t(scale(t(d3)))
  tmp$speech$design$sid <- sid
  list(mat=d3, design=tmp$mouthing$design, roi=c(tmp$mouthing$lh$roi, tmp$mouthing$rh$roi))
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
  d3 <- t(scale(t(d3)))
  tmp$speech$design$sid <- sid
  list(mat=d3, design=tmp$speech$design, roi=c(tmp$speech$lh$roi, tmp$speech$rh$roi))
}

fnames <- paste0(path, paste0(sids, "_betas_lang_sp_and_mo.RDS"))

splist <- lapply(1:length(sids), function(i) load_speech_mat(fnames[i], sids[i]))
molist <- lapply(1:length(sids), function(i) load_mouthing_mat(fnames[i], sids[i]))

deslist_sp <- lapply(splist, "[[", "design")
deslist_mo <- lapply(molist, "[[", "design")

##des <- do.call(rbind, deslist)

run_asca <- function(dlist, rnum,nc=2, deslist, form = ~ vowel * consonant) {
  
  Xlist <- lapply(dlist, function(x) {
    idx <- x$roi %in% rnum
    x$mat[,idx]
  })
  
  #form = ~ vowel * consonant * speaker
  #form = ~ vowel * consonant
  folds <- lapply(deslist, "[[", "block")
  res <- musu_asca(Xlist, form, deslist, ncomp=2, svd.method="base")

  ncomp=nc
  
  sapply(1:length(res$results), function(i) {
    print(names(res$results)[i])
    xx <- res$results[[i]]$bada_result
    class(xx) <- "musu_bada"
    n <- length(levels(xx$Y[[1]]))
    chance <- 1/n
    p <- performance(xx, xx$Y, folds=folds, ncomp=ncomp)
    t.test(sapply(p, function(x) sum(x)/length(x) - chance))$statistic
  })
}

rnums <- sort(unique(splist[[1]]$roi))
rnums <- rnums[rnums < 12000]

tres_sp <- do.call(rbind, lapply(rnums, function(rn) {
  print(rn)
  run_asca(splist, c(rn, rn+1000), nc=2, deslist=deslist_sp, form=~ vowel*consonant)
}))

tres_mo <- do.call(rbind, lapply(rnums, function(rn) {
  print(rn)
  run_asca(molist, c(rn, rn+1000), nc=2, deslist=deslist_mo,form=~ vowel*consonant*speaker)
}))

best <- c(11104, 11133, 11134, 11135, 11136, 11141, 11175)

best <- c(best, best+1000)
tbest1 <- run_asca(splist, best, 1, deslist=deslist_sp)
tbest2 <- run_asca(splist, best, 2, deslist=deslist_sp)
tbest3 <- run_asca(splist, best, 3, deslist=deslist_sp)
tbest4 <- run_asca(splist, best, 5, deslist=deslist_sp)
