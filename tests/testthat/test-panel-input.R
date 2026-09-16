test_that("as_mtr_panel converts complete long data", {
  d <- expand.grid(id = c("a", "b"), time = 1:3, KEEP.OUT.ATTRS = FALSE)
  d$y1 <- seq_len(nrow(d))
  d$y2 <- d$y1 + 10
  arr <- as_mtr_panel(d, id = "id", time = "time", variables = c("y1", "y2"))
  expect_equal(dim(arr), c(2L, 3L, 2L))
  expect_equal(dimnames(arr)[[3]], c("y1", "y2"))
})

test_that("as_mtr_panel rejects duplicate or incomplete panels", {
  d <- expand.grid(id = 1:2, time = 1:2, KEEP.OUT.ATTRS = FALSE)
  d$y <- rnorm(nrow(d))
  expect_error(as_mtr_panel(rbind(d, d[1, ]), "id", "time", "y"), "duplicate")
  expect_error(as_mtr_panel(d[-1, ], "id", "time", "y"), "complete")
})
