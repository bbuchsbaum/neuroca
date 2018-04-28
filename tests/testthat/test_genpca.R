context("genpca")

mat_10_10 <- matrix(rnorm(10*10), 10, 10)

test_that("pca and genpca have same results with identity matrix for row and column constraints", {
  res1 <- genpca(mat_10_10, center=TRUE)
  res2 <- pca(mat_10_10, center=TRUE, ncomp=ncomp(res1))
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(sum(diffscores) < 1e-5)
  expect_equal(singular_values(res1), singular_values(res2))
  
  expect_equal(apply(loadings(res1), 2, function(x) sum(x^2)), rep(1, ncomp(res1)))
})

test_that("gen_pca with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, center=TRUE)
  res2 <- pca(mat_10_10, center=TRUE, scale=TRUE)
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(abs(sum(diffscores)) < 1e-5)
  expect_equal(singular_values(res1), singular_values(res2))
  
})

test_that("gen_pca with dense column and row constraints works", {
  A <- cov(matrix(rnorm(10*10),10,10))
  M <- cov(matrix(rnorm(10*10),10,10))
  res1 <- genpca(mat_10_10, A=A, M=M, center=TRUE)
  expect_equal(ncomp(res1),length(res1$d))
})

test_that("gen_pca with sparse column and row constraints works", {
  A <- neighborweights::similarity_matrix(mat_10_10, k=3)
  diag(A) <- 1
  M <- neighborweights::similarity_matrix(t(mat_10_10), k=3)
  diag(M) <- 1
  res1 <- genpca(mat_10_10, A=A, M=M, center=TRUE)
})

test_that("can truncate a gen_pca model", {
  res1 <- genpca(mat_10_10, center=TRUE)
  res2 <- truncate(res1,5)
  expect_equal(ncomp(res2), 5)
})

test_that("can project a row vector", {
  A <- cov(matrix(rnorm(10*10),10,10))
  M <- cov(matrix(rnorm(10*10),10,10))
  
  res1 <- genpca(mat_10_10, A=A, M=M, center=TRUE)
  p <- project(res1, mat_10_10[1,])
  expect_equal(dim(p), c(1,ncomp(res1)))
})

test_that("can extract residuals", {
  res1 <- genpca(mat_10_10, center=TRUE)
  resid <- residuals(res1, ncomp=2, mat_10_10)
  expect_equal(sum(res1$d[3:length(res1$d)] ^2), sum(resid^2))
})
