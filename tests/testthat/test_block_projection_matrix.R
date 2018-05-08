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

test_that("can project newdata through a block_projection_matrix", {
  
  mat1 <- matrix(rnorm(10*10), 10, 10)
  mat2 <- matrix(rnorm(10*15), 10,15)
  
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(mat2, ncomp=3)
  
  bm <- block_projection_matrix(list(pca1,pca2))
  mat  <- matrix(rnorm(10*25), 10, 25)
  p <- project(bm, mat)
  expect_equal(dim(p), c(10,7))
})

test_that("can project a newdata sub-block through a block_projection_matrix", {
  
  mat1 <- matrix(rnorm(10*10), 10, 10)
  mat2 <- matrix(rnorm(10*15), 10,15)
  
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(mat2, ncomp=3)
  
  bm <- block_projection_matrix(list(pca1,pca2))
  mat  <- matrix(rnorm(10*10), 10, 10)
  p1 <- project(bm, mat, block_index=1)
  p2 <- project(bm, matrix(rnorm(10*15), 10, 15), block_index=2)
  expect_equal(dim(p1), c(10,4))
  expect_equal(dim(p2), c(10,3))
})