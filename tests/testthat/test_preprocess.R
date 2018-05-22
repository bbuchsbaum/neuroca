test_that("can preprocess a matrix no center, no scale", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- pre_processor(mat1, center=FALSE, scale=FALSE)
  x <- pre_process(pp, mat1)
  x2 <- reverse_pre_process(pp, x)
  expect_equal(x,x2)
  expect_equal(mat1,x)
  expect_equal(x, pp$Xp)
})

test_that("can preprocess a matrix center only", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- pre_processor(mat1, center=TRUE, scale=FALSE)
  x <- pre_process(pp, mat1)
  x2 <- reverse_pre_process(pp, x)
  expect_equal(x2,mat1)
})



test_that("can preprocess a matrix center and scale", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- pre_processor(mat1, center=TRUE, scale=TRUE)
  x <- pre_process(pp, mat1)
  x2 <- reverse_pre_process(pp, x)
  expect_equal(x2,mat1)
})


test_that("can preprocess a matrix with a subind", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- pre_processor(mat1,center=TRUE, scale=FALSE)
  
  res <- pre_process(pp, newdata=mat1[,1:2], subind=1:2)
  
  expect_equal(res, scale(mat1[,1:2], center=TRUE))
})


test_that("can preprocess a block projector", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  mat2 <-  matrix(rnorm(10*10), 10, 10)
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(mat2, ncomp=2)
  
  bm <- block_projector(list(pca1,pca2))
  pp <- pre_processor(bm,center=FALSE, scale=FALSE)
  pdat <- pre_process(pp)
  expect_equal(ncol(pdat), 6)
  expect_equal(project(bm), pdat)
})

test_that("can preprocess a block projector with newdata", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  mat2 <-  matrix(rnorm(10*10), 10, 10)
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(mat2, ncomp=2)
  
  bm <- block_projector(list(pca1,pca2))
  pp <- pre_processor(bm,center=FALSE, scale=FALSE)
  
  mat3 <- cbind(mat1,mat2)
  pdat <- pre_process(pp,mat3)
  
  expect_equal(ncol(pdat), 6)
  expect_equal(project(bm), pdat)
})

test_that("can preprocess a block projector with newdata from a sub-block", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  mat2 <-  matrix(rnorm(10*10), 10, 10)
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(mat2, ncomp=2)
  
  bm <- block_projector(list(pca1,pca2))
  pp <- pre_processor(bm,center=FALSE, scale=FALSE)
  
  mat3 <- cbind(mat2)
  pdat <- pre_process(pp,mat3, block_index=2)
  
  expect_equivalent(project(bm, block_index=2), unclass(pdat))
})






