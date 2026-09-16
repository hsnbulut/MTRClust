test_that("package wrapper is equivalent to authority prototype on a small deterministic fit", {
  auth <- new.env(parent = globalenv())
  sys.source("prototype_authority/mtrclust_prototype.R", envir = auth)
  dat <- simulate_mtr_panel(n = 24, tt = 12, p = 3, seed = 505)
  args <- list(k = 3, alpha = 0.05, latent_rank = 2L, spar = 0.55,
               method = "v6_global", nstart = 2L, max_iter = 4L,
               min_cluster = NULL, time_fraction_cut = 0.50,
               global_fraction_cut = 0.50, repair = TRUE,
               rollback_nonconverged = TRUE, selection_rule = "objective",
               cycle_detection = TRUE, cycle_max_period = 6L,
               cycle_tol = 1e-7, seed = 505L, verbose = FALSE,
               return_candidates = TRUE)
  fit_auth <- do.call(auth$mtrclust_v12, c(list(x = dat$x), args))
  fit_pkg <- mtrclust(dat$x, k = 3, alpha = 0.05, rank = 2L, spar = 0.55,
                      method = "global", nstart = 2L, max_iter = 4L,
                      min_cluster = NULL, time_fraction_cut = 0.50,
                      global_fraction_cut = 0.50, repair = TRUE,
                      rollback_nonconverged = TRUE,
                      selection_rule = "objective", cycle_detection = TRUE,
                      cycle_max_period = 6L, cycle_tol = 1e-7,
                      seed = 505L, verbose = FALSE, return_candidates = TRUE)
  expect_equal(fit_pkg$objective, fit_auth$objective, tolerance = 1e-10)
  expect_equal(canonical_partition(fit_pkg$cluster), canonical_partition(fit_auth$cluster))
  expect_equal(fit_pkg$keep, fit_auth$keep)
  expect_equal(fit_pkg$observation_outliers, fit_auth$observation_outliers)
  expect_equal(fit_pkg$termination, fit_auth$termination)
  expect_equal(fit_pkg$converged, fit_auth$converged)
  expect_equal(fit_pkg$rolled_back, fit_auth$rolled_back)
  expect_equal(fit_pkg$iterations, fit_auth$iterations)
  expect_equal(fit_pkg$best_iteration, fit_auth$best_iteration)
  expect_equal(nrow_test(fit_pkg$isolated_cell_outliers), nrow_test(fit_auth$isolated_cell_outliers))
  expect_equal(nrow_test(fit_pkg$sensor_segments), nrow_test(fit_auth$sensor_segments))
  expect_equal(nrow_test(fit_pkg$multivariate_segments), nrow_test(fit_auth$multivariate_segments))
})

test_that("legacy v6_global alias is identical to global", {
  dat <- simulate_mtr_panel(n = 24, tt = 12, p = 3, seed = 606)
  fit_global <- mtrclust(dat$x, k = 3, alpha = 0.05, rank = 2L,
                         method = "global", nstart = 2L, max_iter = 4L,
                         seed = 606L)
  fit_legacy <- mtrclust(dat$x, k = 3, alpha = 0.05, rank = 2L,
                         method = "v6_global", nstart = 2L, max_iter = 4L,
                         seed = 606L)
  expect_equal(fit_global$objective, fit_legacy$objective, tolerance = 1e-10)
  expect_equal(canonical_partition(fit_global$cluster), canonical_partition(fit_legacy$cluster))
  expect_equal(fit_global$keep, fit_legacy$keep)
  expect_equal(fit_global$observation_outliers, fit_legacy$observation_outliers)
  expect_equal(fit_global$termination, fit_legacy$termination)
  expect_equal(fit_global$converged, fit_legacy$converged)
  expect_equal(fit_global$rolled_back, fit_legacy$rolled_back)
  expect_equal(fit_global$iterations, fit_legacy$iterations)
  expect_equal(fit_global$best_iteration, fit_legacy$best_iteration)
})
