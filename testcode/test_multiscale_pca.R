library(neuroim2)
library(fastcluster)
mask <- read_vol("testdata/global_mask.nii")
vec <- read_vec("testdata/betas_1.nii")

mask.idx <- which(mask>0)
X <- series(vec, mask.idx)

grid <- index_to_coord(neuroim2::space(mask), mask.idx)
hclus <- hclust.vector(grid, method="centroid")
cuts <- c(16,32,64,128,256)