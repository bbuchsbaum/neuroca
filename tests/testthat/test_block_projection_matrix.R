test_that("can construct a block_projection_matrix", {
  
  mat1 <- matrix(rnorm(10*10), 10, 10)
  mat2 <- matrix(rnorm(10*15), 10,15)
  
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(mat2, ncomp=3)
  
  bm <- block_projection_matrix(list(pca1,pca2))
  expect_equal(ncomp(bm), 7)
  expect_equal(ncol(bm), 25)
  expect_equal(nrow(bm), 10)
  expect_equal(nblocks(bm), 2)
  expect_equal(block_lengths(bm), c(10,15))
  
  pm <- project(bm)
  expect_equal(ncol(pm), 7)
  
})