test_that("input validation rejects invalid arrays", {
  expect_error(mtrclust(matrix(1:6, 2, 3), k = 2), "three-dimensional")
  x <- array(rnorm(24), c(4, 3, 2))
  x[1, 1, 1] <- NA_real_
  expect_error(mtrclust(x, k = 2), "finite")
  x <- array(1, c(5, 3, 2))
  expect_error(mtrclust(x, k = 2), "constant")
})

test_that("k alpha rank and iteration validation are explicit", {
  x <- array(rnorm(5 * 4 * 3), c(5, 4, 3))
  expect_error(mtrclust(x, k = 1), "`k`")
  expect_error(mtrclust(x, k = 6), "`k`")
  expect_error(mtrclust(x, k = 2, alpha = 0.5), "`alpha`")
  expect_error(mtrclust(x, k = 2, rank = 0), "`rank`")
  expect_error(mtrclust(x, k = 2, nstart = 0), "`nstart`")
  expect_error(mtrclust(x, k = 2, max_iter = 0), "`max_iter`")
  expect_error(mtrclust(x, k = 2, seed = NA), "`seed`")
  expect_error(mtrclust(x, k = 2, cutoff_bounds = c(3, 2)), "`cutoff_bounds`")
})
