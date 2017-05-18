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
  v1 <- load_mat(paste0(sid, "_ventral_video_1.RDS"))
  v2 <- load_mat(paste0(sid, "_ventral_video_2.RDS"))
  v3 <- load_mat(paste0(sid, "_ventral_video_3.RDS"))
  
  mat <- rbind(v1$mat, v2$mat, v3$mat)
  
  keep <- apply(mat,2, function(x) sum(x==0)) == 0
  mat <- mat[,keep]
  
  
  ncomp <- fast_estim_ncomp(scale(mat), 1, 300)$bestcomp
  des <- rbind(v1$design, v2$design, v3$design)
  
  list(sid=sid, mat=mat, design=des, ncomp=ncomp)
})