library(neuroim)

Xlist <- readRDS("~/Dropbox/Brad/neuropls/speech/sp_matrices.rds")
design <- read.table("~/Dropbox/Brad/neuropls/speech/sp_syllable_design.txt", header=TRUE)
formula = ~ consonant*vowel*speaker