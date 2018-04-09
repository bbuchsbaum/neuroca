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



