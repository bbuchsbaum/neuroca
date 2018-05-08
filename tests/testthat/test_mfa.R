context("mfa")

block_mat <- block_matrix(lapply(1:3, function(i) matrix(rnorm(10*10), 10, 10)))


test_that("mfa with no scaling is identical to pca", {
  res1 <- mfa(block_mat, center=TRUE, normalization="None")
  res2 <- pca(as.matrix(block_mat), center=TRUE, ncomp=ncomp(res1))
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(sum(diffscores) < 1e-5)
  expect_equal(singular_values(res1), singular_values(res2))
  
  expect_equal(apply(loadings(res1), 2, function(x) sum(x^2)), rep(1, ncomp(res1)))

})

test_that("can project a row vector", {
  res1 <- mfa(block_mat, center=TRUE, ncomp=5, normalization="None")
  x <- rnorm(ncol(block_mat))
  proj <- project(res1, x, comp=1:3)
  expect_equal(ncol(proj), 3)
  expect_equal(nrow(proj), 1)
})

test_that("can project a matrix for a single table", {
  res1 <- mfa(block_mat, center=TRUE, ncomp=5, normalization="None")
  x <- matrix(rnorm(10*10), 10, 10)
  proj <- project(res1, x, comp=1:3, block_index=1)
  expect_equal(ncol(proj), 3)
  expect_equal(nrow(proj), nrow(x))
})

test_that("can run mfa with a block_projection_matrix", {
  bm <- block_projection_matrix(lapply(1:5, function(i) pca(matrix(rnorm(10*10), 10, 10), ncomp=4)))
  
  res1 <- mfa(bm, center=TRUE, ncomp=5, normalization="None")
  
  x <- matrix(rnorm(10*10), 10, 10)
  proj <- project(res1, x, comp=1:5, block_index=1)
  
  x2 <- matrix(rnorm(10*50), 10, 50)
  proj_all  <- project(res1, x2)
  
  expect_equal(ncol(proj), 5)
  expect_equal(nrow(proj), nrow(x))
  
  expect_equal(ncol(proj_all), 5)
})

test_that("can run mfa with a block_projection_matrix with matrices with different number of columns", {
  bm <- block_projection_matrix(lapply(8:12, function(i) pca(matrix(rnorm(10*i), 10, i), ncomp=4)))
  
  res1 <- mfa(bm, center=TRUE, ncomp=5, normalization="None")
  
  x <- matrix(rnorm(10*8), 10, 8)
  proj <- project(res1, x, comp=1:5, block_index=1)
  
  x2 <- matrix(rnorm(10*sum(8:12)), 10, sum(8:12))
  proj_all  <- project(res1, x2)
  
  expect_equal(ncol(proj), 5)
  expect_equal(nrow(proj), nrow(x))
  
  expect_equal(ncol(proj_all), 5)
})




