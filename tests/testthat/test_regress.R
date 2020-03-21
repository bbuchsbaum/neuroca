test_that("can carry out regress analysis", {
  basis <- matrix(rnorm(15*10), 15, 10)
  X <- matrix(rnorm(15*20),15,20)
  
  p <- regress(basis, X)
  expect_true(!is.null(p))
})

test_that("can project a regress", {
  basis <- matrix(rnorm(15*10), 15, 10)
  X <- matrix(rnorm(15*20),15,20)
  
  p <- regress(basis, X)
  scores <- project(p,X)
  expect_true(!is.null(scores))
  expect_equal(dim(scores), dim(basis))
})

test_that("can reconstruct a regress", {
  basis <- matrix(rnorm(15*10), 15, 10)
  X <- matrix(rnorm(15*20),15,20)
  
  p <- regress(basis, X)
  Xr <- reconstruct(p)
})

