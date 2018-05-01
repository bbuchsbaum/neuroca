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
  pca2 <- pca(pca1, ncomp=2)
  
  expect_equal(ncol(pca2), ncol(pca1))
  expect_equal(ncomp(pca2), 2)
  expect_equal(abs(project(pca1, comp=1)), abs(project(pca2, comp=1)))
})

test_that("can run a triply nested pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  pca2 <- pca(pca1, ncomp=3)
  pca3 <- pca(pca2, ncomp=2)
  expect_equal(ncol(pca3), ncol(mat1))
 
})

test_that("can partially project a plain pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=4)
  project(pca1, mat1[,1:2], subind=1:2)
  
})

test_that("can truncate a pca", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  pca1 <- pca(mat1, ncomp=8)
  pca2 <- truncate(pca1, ncomp=5)
  
  expect_true(ncol(scores(pca2)) == 5)
  expect_true(ncol(project(pca2)) == 5)
  expect_equal(dim(residuals(pca2, ncomp=5, xorig=mat1)), c(10,15))
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
