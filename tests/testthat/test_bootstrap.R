test_that("can bootstrap a pca model", {
  mat1 <- matrix(rnorm(100*200), 100, 200)
  pca1 <- pca(mat1, ncomp=10, preproc=center())
  
  bootres <- bootstrap(pca1, nboot=20, k=3)
  expect_true(!is.null(bootres))
})

test_that("can bootstrap a genpca model", {
  mat1 <- matrix(rnorm(100*200), 100, 200)
  M <- runif(100)
  pca1 <- genpca(mat1, M=Matrix::Diagonal(x=M), ncomp=10, preproc=center())
  
  bootres <- bootstrap(pca1, nboot=20, k=3)
  expect_true(!is.null(bootres))
})