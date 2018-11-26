library(neighborweights)
library(dplyr)
library(Matrix)

base_path <- "/Users/bbuchsbaum/analysis/airplane_paper/gpca"

dat <- readRDS(paste0(base_path, "/", "betas_loc_mvideo.rds"))
sids <- sapply(dat, function(x) x$des$sid[1])
old <- substr(sids, 1,1) == "2"
young <- substr(sids, 1,1) != "2"

master_df1 <- do.call(rbind, lapply(dat, function(x) {
  ret <- x$des %>% mutate(id=row_number()) %>% group_by(Condition, repnum,sid) %>% do({
    ord <- order(.$Video)
    id <- .$id[ord]
    img <- x$mat[id,]
    img <- lapply(1:nrow(img), function(i) {
      vals <- img[i,]
      vals <- vals - mean(vals)
      vals/sqrt(sum(vals^2))
    })
    
    tibble(Video=.$Video[ord], repnum=.$repnum[ord], sid=.$sid[1], rowid=id, img=img)
  })
}))

df_proj <- lapply(1:21, function(i) {
  print(i)
  
  dftrain <- master_df1 %>% filter(repnum != i)
  dftest <- master_df1 %>% filter(repnum == i)
  
  Xtr <- dftrain %>% group_by(Condition, sid) %>% do( {
    X <- do.call(rbind, .$img)
    Y <- .$Video
    tibble(X=list(X), Y=list(Y), Condition=.$Condition[1], repnum=.$repnum[1])
  }) %>% ungroup() %>% mutate(block=row_number())
  
  Xte <- dftest %>% group_by(Condition, sid) %>% do( {
    X <- do.call(rbind, .$img)
    Y <- .$Video

    tibble(X=list(X), Y=list(Y), Condition=.$Condition[1], repnum=.$repnum[1])
  }) %>% ungroup() %>% mutate(block=row_number())
  
  Xl <- Xtr$X
  Yl <- Xtr$Y
  
  mures <- mubada(Yl, Xl, ncomp=11, normalization="None")
  
  Xte %>% rowwise() %>% do({
    xp <- block_project(mures, Xte$X[[.$block]], .$block)
    browser()
    cvals <- sapply(1:ncol(xp), function(i) {
      sum(xp[,i] * mures$scores[,i]) / (sqrt(sum(xp[,i]^2))*sqrt(sum(mures$scores[,i]^2)))
    })
    
    group <- ifelse(substr(.$sid[1],1,1) == "1", "young", "old")
    
    if (substr(.$sid[1],1,1) == "5") {
      group = "patient"
    }
    
    
    rv <- FactoMineR::coeffRV(as.matrix(xp), as.matrix(mures$scores))$rvstd
    tibble(cosvals=c(rv, cvals), dim=0:length(cvals), repnum=.$repnum[1], sid=.$sid[1], condition=.$Condition[1], 
           group=group)
    
  })
    
})
  
df_proj <- do.call(rbind, df_proj)
writeRDS(df_proj, "mubada_proj.rds")

#sc=scorepred(t(scores(mures)), t(xx), type="cosine")
#project_cols(mures, xx)
