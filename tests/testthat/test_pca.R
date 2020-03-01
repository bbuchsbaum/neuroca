test_that("can run a simple pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1)
  
  expect_equal(length(singular_values(pca1)), ncomp(pca1))
  expect_equal(scores(pca1),project(pca1))
  expect_equal(ncol(pca1), ncol(mat1))
})

test_that("can run a simple nested pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(scores(pca1), ncomp=2)
  
  cp <- compose(pca1, pca2)
  
  expect_equal(ncol(cp), ncol(pca1))
  expect_equal(ncomp(cp), 2)
  expect_equal(dim(reconstruct(cp)), dim(mat1))
})

test_that("can run a triply nested pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(scores(pca1), ncomp=3)
  pca3 <- pca(scores(pca2), ncomp=2)
  cp <- compose_all(pca1,pca2,pca3)
  
  expect_equal(ncomp(pca3), ncomp(cp))
  expect_equal(dim(project(cp)), c(10,2))
})

test_that("can partially project a plain pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  p <- project(pca1, mat1[,1:2], colind=1:2)
  expect_equal(dim(p), c(10, 4))
})

test_that("can truncate a pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=8)
  pca2 <- truncate(pca1, ncomp=5)
  
  expect_true(ncol(scores(pca2)) == 5)
  expect_true(ncol(project(pca2)) == 5)
  expect_equal(dim(residuals(pca2, ncomp=5, xorig=mat1)), c(10,15))
})

test_that("can reconstruct a pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=15)
  recon <- reconstruct(pca1)
  expect_equal(recon, mat1)
})

test_that("can reconstruct a pca on a correlation matrix", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  cmat <- cor(mat1)
  pca1 <- pca(cmat, ncomp=15, preproc=pass())
  recon <- reconstruct(pca1)
  expect_equal(recon, cmat)
})


test_that("can reconstruct a pca and get original data", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=15)
  sc <- scores(pca1)
  ip <- reconstruct(pca1)
  ip2 <- reconstruct(pca1, newdata=sc)
  expect_equal(ip, mat1)
})

test_that("can reconstruct a nested pca and recover original data", {
  mat1 <- matrix(rnorm(10*20), 10, 20)
  pca1 <- pca(mat1, ncomp=10)
  pca2 <- pca(scores(pca1), ncomp=10)
  cp <- compose(pca1,pca2)
  ip <- reconstruct(cp)
  expect_equal(ip, mat1)
})


test_that("can project a single variable onto pca space", {
  mat1 <- matrix(rnorm(10*20), 10, 20)
  pca1 <- pca(mat1, ncomp=10, center=TRUE)
  v1 <- mat1[,1]
  proj <- project_cols(pca1,newdata=v1)
  expect_equal(as.vector(proj), as.vector(pca1$v[1,]))
})

test_that("can project a set of variables onto pca space", {
  mat1 <- matrix(rnorm(10*20), 10, 20)
  mat2 <- matrix(rnorm(10*100), 10, 100)
  pca1 <- pca(mat1, ncomp=10, center=TRUE)
  
  proj <- project_cols(pca1,newdata=mat2)
  expect_equal(dim(proj), c(100, pca1$ncomp))
})


test_that("can run a shrink_pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- shrink_pca(mat1)
  
  
  expect_true(!is.null(scores(pca1)))
  expect_true(!is.null(loadings(pca1)))
  expect_true(!is.null(ncomp(pca1)))
  expect_true(!is.null(project(pca1)))
  expect_true(!is.null(project(pca1, mat1)))
              
})
