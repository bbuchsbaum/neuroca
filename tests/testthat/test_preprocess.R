test_that("can preprocess a matrix no center, no scale", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- pass()
  x <- prep(pp, mat1)
  x2 <- x$reverse_transform(x$Xp)
  expect_equal(mat1,x2)
  expect_equal(x$Xp, mat1)
})

test_that("can preprocess a matrix center only", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- center()
  x <- prep(pp, mat1)
  x2 <- x$reverse_transform(x$Xp)
  expect_equal(mat1,x2)
  expect_true(all(mat1 != x$Xp))
})

test_that("can preprocess a matrix with column scaling", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  wts <- 2:16
  pp <- colscale(type="weights", weights=wts)
  x <- prep(pp, mat1)
  x2 <- x$reverse_transform(x$Xp)
  expect_equal(mat1,x2)
  expect_true(all(mat1 != x$Xp))
})

test_that("can reset a prepper with `fresh`", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- center()
  x <- prep(pp, mat1)
  
})



test_that("can preprocess a matrix center and scale", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- standardize()
  x <- prep(pp, mat1)
  x2 <- x$reverse_transform(x$Xp)
  expect_equal(mat1,x2)
  expect_true(all(mat1 != x$Xp))
})

test_that("can compose two pre-processors", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- center() %>% colscale(type="z")
  pp2 <- standardize()
  x <- prep(pp, mat1)
  x0 <- prep(pp2, mat1)
  x2 <- x$reverse_transform(x$Xp)
  expect_equal(mat1,x2)
  expect_equal(x$Xp,x0$Xp)
  expect_true(all(mat1 != x$Xp))
})



test_that("can preprocess a matrix with a colind", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- center()
  
  x <- prep(pp, mat1)
  ret <- x$transform(mat1[,1:2], 1:2)
  
  expect_equal(ret, x$Xp[,1:2])
})

# 
# test_that("can preprocess a block projector", {
#   mat1 <- matrix(rnorm(10*15), 10, 15)
#   mat2 <-  matrix(rnorm(10*10), 10, 10)
#   pca1 <- pca(mat1, ncomp=4)
#   pca2 <- pca(mat2, ncomp=2)
#   
#   bm <- block_projector(list(pca1,pca2))
#   pp <- pre_processor(bm,center=FALSE, scale=FALSE)
#   pdat <- pre_process(pp)
#   expect_equal(ncol(pdat), 6)
#   expect_equal(project(bm), pdat)
# })
# 
# test_that("can preprocess a block projector with newdata", {
#   mat1 <- matrix(rnorm(10*15), 10, 15)
#   mat2 <-  matrix(rnorm(10*10), 10, 10)
#   pca1 <- pca(mat1, ncomp=4)
#   pca2 <- pca(mat2, ncomp=2)
#   
#   bm <- block_projector(list(pca1,pca2))
#   pp <- pre_processor(bm,center=FALSE, scale=FALSE)
#   
#   mat3 <- cbind(mat1,mat2)
#   pdat <- pre_process(pp,mat3)
#   
#   expect_equal(ncol(pdat), 6)
#   expect_equal(project(bm), pdat)
# })
# 
# test_that("can preprocess a block projector with newdata from a sub-block", {
#   mat1 <- matrix(rnorm(10*15), 10, 15)
#   mat2 <-  matrix(rnorm(10*10), 10, 10)
#   pca1 <- pca(mat1, ncomp=4)
#   pca2 <- pca(mat2, ncomp=2)
#   
#   bm <- block_projector(list(pca1,pca2))
#   pp <- pre_processor(bm,center=FALSE, scale=FALSE)
#   
#   mat3 <- cbind(mat2)
#   pdat <- pre_process(pp,mat3, block_index=2)
#   
#   expect_equivalent(project(bm, block_index=2), unclass(pdat))
# })






