

test_that("can project a simple matrix", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  expect_equal(project(mat1), mat1)
})

test_that("can partially project a simple matrix", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  proj <- project(mat1, subind=1:2)
  expect_equivalent(as.matrix(proj[,1:2]), mat1[,1:2])
})

test_that("can partially project a vector", {
  mat1 <- matrix(rnorm(10*15), 10, 15)
  proj <- project(mat1, subind=1)
  expect_equivalent(as.matrix(proj[,1]), mat1[,1])
})
