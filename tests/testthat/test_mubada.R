context("mubada")


test_that("can compute a simple mubada analysis", {
  xlist <- lapply(1:10, function(i) matrix(rnorm(20*10),20,10))
  y <- lapply(1:10, function(i) factor(sample(letters[1:3], 20, replace=TRUE)))
  res1 <- mubada(y, xlist, center=TRUE, normalization="None")
  expect_equal(nrow(scores(res1)),3)
  R <- reconstruct(res1)
  expect_equal(dim(R), c(3, 100))
})

test_that("can compute a mubada analysis with a projector list", {
  xlist <- lapply(1:10, function(i) pca(matrix(rnorm(20*10),20,10), ncomp=3))
  y <- lapply(1:10, function(i) factor(sample(letters[1:3], 20, replace=TRUE)))
  res1 <- mubada(y, xlist, center=TRUE, normalization="None")
  expect_equal(nrow(scores(res1)),3)
  R <- reconstruct(res1)
  expect_equal(dim(R), c(3, 100))
})