context("mubada")


test_that("can compute a simple mubada analysis", {
  xlist <- lapply(1:10, function(i) matrix(rnorm(20*10),20,10))
  y <- lapply(1:10, function(i) factor(sample(letters[1:3], 20, replace=TRUE)))
  res1 <- mubada(y, xlist, center=TRUE, normalization="None")
  expect_equal(nrow(scores(res1)),3)
  R <- reconstruct(res1)
  expect_equal(dim(R), c(3, 100))
})

test_that("can compute a mubada analysis on a list of scores", {
  xlist <- lapply(1:10, function(i) scores(pca(matrix(rnorm(20*10),20,10), ncomp=3)))
  y <- lapply(1:10, function(i) factor(sample(letters[1:3], 20, replace=TRUE)))
  res1 <- mubada(y, xlist, center=TRUE, normalization="None")
  expect_equal(nrow(scores(res1)),3)
  R <- reconstruct(res1)
  expect_equal(dim(R), c(3, 30))
})

test_that("can compute a bada analysis on a basic input matrix", {
  X <- matrix(rnorm(100*1000), 100, 1000)
  Y <- factor(rep(letters[1:4], length.out=100))
  S <- factor(rep(1:10, each=10))
  
  bres <- bada(Y, X, S, ncomp=3)
  expect_true(!is.null(bres))
})

test_that("can bootstrap a bada analysis", {
  X <- matrix(rnorm(100*1000), 100, 1000)
  Y <- factor(rep(letters[1:4], length.out=100))
  S <- factor(rep(1:10, each=10))
  
  bres <- bada(Y, X, S, ncomp=3)
  res <- bootstrap(bres, type="projection", nboot=50)
})