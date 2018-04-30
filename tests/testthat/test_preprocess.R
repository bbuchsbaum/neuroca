test_that("can preprocess a matrix no center, no scale", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pp <- pre_processor(mat1, center=FALSE, scale=FALSE)
  x <- pre_process(pp, mat1)
  x2 <- reverse_pre_process(pp, x)
  expect_equal(x,x2)
  expect_equal(mat1,x)
  
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

test_that("can preprocess a pca projector", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  pp <- pre_processor(pca1,scale=FALSE)
  x <- pre_process(pp, mat1)
  expect_equal(x,project(pca1))
})

test_that("can preprocess a pca projector with a subind", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  pp <- pre_processor(pca1,scale=FALSE)
  
  expect_equal(x,project(pca1))
})




