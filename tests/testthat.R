if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(MTRClust)
  test_check("MTRClust")
} else {
  library(MTRClust)
  dat <- simulate_mtr_panel(n = 24, tt = 12, p = 3, seed = 101)
  fit <- mtrclust(dat$x, k = 3, alpha = 0.05, nstart = 2, max_iter = 4, seed = 101)
  stopifnot(inherits(fit, "mtrclust"))
  stopifnot(length(fit$cluster) == dim(dat$x)[1L])
  stopifnot(is.character(fit$termination))
  stopifnot(is.logical(fit$converged))
  stopifnot(is.logical(fit$rolled_back))
  invisible(summary(fit))
}
