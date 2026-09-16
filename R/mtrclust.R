#' Fit MTRClust to a multivariate panel array
#'
#' `mtrclust()` clusters numeric multivariate panel trajectories stored as an
#' `n x T x p` array. It is a user-facing wrapper around the audited MTRClust
#' prototype implementation used in the accompanying manuscript.
#'
#' The algorithm robustly scales variables, initializes a retained set from
#' profile outlyingness, alternates between robust assignment/trimming and
#' smooth low-rank trajectory updates, and records multilevel diagnostic
#' outputs. The returned solution may terminate by convergence, cycle detection,
#' or the iteration limit. Reaching `max_iter` is not convergence; when rollback
#' is enabled, the returned nonconverged solution is restored to the best
#' recorded objective iterate.
#'
#' Missing and non-finite values are not currently supported by the package
#' wrapper. They are rejected explicitly rather than imputed.
#'
#' @param x Numeric three-dimensional array with dimensions `n x T x p`.
#' @param k Number of clusters.
#' @param alpha Observation trimming proportion in `[0, 0.5)`.
#' @param rank Low-rank reconstruction dimension.
#' @param spar Smoothing parameter passed to [stats::smooth.spline()].
#' @param method Residual scaling mode. The default `"global"` uses the
#'   validated global residual scale. `"v6_global"` is accepted only as a
#'   legacy alias for `"global"` and maps to the same internal implementation.
#' @param nstart Number of starts.
#' @param max_iter Maximum number of iterations per start.
#' @param tol_center Relative center-change tolerance for convergence.
#' @param min_cluster Minimum retained units per cluster. If `NULL`, the
#'   authority default is used.
#' @param min_segment Minimum run length for segment diagnostics.
#' @param c_cell Tukey bisquare tuning constant.
#' @param cell_threshold Cell-weight threshold used for raw cell flags.
#' @param time_fraction_cut Unit-time fraction threshold for time diagnostics.
#' @param global_fraction_cut Retained-unit fraction threshold for global
#'   sensor-time diagnostics.
#' @param local_trim Local cell-loss trimming fraction inside each profile loss.
#' @param cutoff_prob Quantile probability for adaptive residual cutoffs.
#' @param cutoff_bounds Lower and upper bounds for adaptive residual cutoffs.
#' @param scale_floor Minimum residual scale used by cluster-scale modes.
#' @param repair Whether to apply minimum retained cluster-size repair.
#' @param rollback_nonconverged Whether to roll nonconverged candidates back to
#'   their best objective iterate.
#' @param trace_history Whether to retain per-iteration trace data.
#' @param selection_rule Candidate-selection rule. The current authority
#'   implementation supports objective-based selection.
#' @param cycle_detection Whether to detect repeated objective cycles.
#' @param cycle_max_period Maximum objective-cycle period to test.
#' @param cycle_tol Relative tolerance for objective-cycle detection.
#' @param seed Integer random seed.
#' @param verbose Whether to print per-iteration progress messages.
#' @param return_candidates Whether to retain all valid internal candidates.
#'
#' @return An S3 object of class `"mtrclust"` containing cluster assignments,
#'   retained indicators, observation outliers, fitted trajectories, objective
#'   history, termination and rollback diagnostics, candidate summaries, and
#'   multilevel outlier diagnostics.
#' @seealso [summary.mtrclust()], [plot.mtrclust()],
#'   [mtrclust_diagnostics()], [as_mtr_panel()]
#' @importFrom stats as.dist cutree hclust IQR kmeans median predict quantile
#'   rnorm smooth.spline
#' @importFrom utils tail
#' @examples
#' dat <- simulate_mtr_panel(n = 30, tt = 15, p = 3, seed = 1)
#' fit <- mtrclust(dat$x, k = 3, alpha = 0.05, rank = 2,
#'                 nstart = 2, max_iter = 20, seed = 1)
#' fit
#' summary(fit)
#' plot(fit, type = "objective")
#' @export
mtrclust <- function(x, k, alpha = 0.08, rank = 2L, spar = 0.55,
                     method = c("global", "cluster_scale", "cluster_adaptive", "v6_global"),
                     nstart = 8L, max_iter = 100L, tol_center = 1e-5,
                     min_cluster = NULL, min_segment = 3L, c_cell = 4.685,
                     cell_threshold = 0.20, time_fraction_cut = 0.50,
                     global_fraction_cut = 0.50, local_trim = 0.10,
                     cutoff_prob = 0.99, cutoff_bounds = c(2.5, 3.5),
                     scale_floor = 0.10, repair = TRUE,
                     rollback_nonconverged = TRUE, trace_history = FALSE,
                     selection_rule = "objective", cycle_detection = TRUE,
                     cycle_max_period = 6L, cycle_tol = 1e-7,
                     seed = 20260827L, verbose = FALSE,
                     return_candidates = FALSE) {
  call <- match.call()
  method <- match.arg(method)
  method_internal <- if (identical(method, "global")) "v6_global" else method
  validate_mtrclust_input(
    x = x, k = k, alpha = alpha, rank = rank, spar = spar,
    nstart = nstart, max_iter = max_iter, tol_center = tol_center,
    min_cluster = min_cluster, min_segment = min_segment,
    cutoff_bounds = cutoff_bounds, seed = seed
  )

  fit <- mtrclust_v12(
    x = x, k = as.integer(k), alpha = alpha, latent_rank = as.integer(rank),
    spar = spar, c_cell = c_cell, cell_threshold = cell_threshold,
    time_fraction_cut = time_fraction_cut,
    global_fraction_cut = global_fraction_cut, local_trim = local_trim,
    nstart = as.integer(nstart), max_iter = as.integer(max_iter),
    tol_center = tol_center, min_cluster = min_cluster,
    min_segment = as.integer(min_segment), method = method_internal,
    cutoff_prob = cutoff_prob, cutoff_bounds = cutoff_bounds,
    scale_floor = scale_floor, repair = repair,
    rollback_nonconverged = rollback_nonconverged,
    trace_history = trace_history, selection_rule = selection_rule,
    cycle_detection = cycle_detection,
    cycle_max_period = as.integer(cycle_max_period),
    cycle_tol = cycle_tol, seed = as.integer(seed), verbose = verbose,
    return_candidates = return_candidates
  )

  fit$call <- call
  fit$arguments <- list(
    k = as.integer(k), alpha = alpha, rank = as.integer(rank), spar = spar,
    method = method, method_internal = method_internal,
    nstart = as.integer(nstart),
    max_iter = as.integer(max_iter), min_cluster = min_cluster,
    seed = as.integer(seed)
  )
  fit$input_dim <- dim(x)
  fit$metadata <- list(
    package = "MTRClust",
    package_version = tryCatch(as.character(utils::packageVersion("MTRClust")),
                               error = function(e) NA_character_),
    r_version = R.version.string
  )
  class(fit) <- c("mtrclust", setdiff(class(fit), "mtrclust"))
  fit
}

validate_mtrclust_input <- function(x, k, alpha, rank, spar, nstart, max_iter,
                                    tol_center, min_cluster, min_segment,
                                    cutoff_bounds, seed) {
  if (!is.array(x) || length(dim(x)) != 3L) {
    stop("`x` must be a numeric three-dimensional array with dimensions n x T x p.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite numeric values; missing-value imputation is not supported.", call. = FALSE)
  }
  dims <- dim(x)
  if (any(dims < 2L)) {
    stop("`x` must have at least two units, two time points, and two variables.", call. = FALSE)
  }
  if (!is_scalar_integerish(k) || k < 2L || k > dims[1]) {
    stop("`k` must be an integer between 2 and the number of units.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha >= 0.5) {
    stop("`alpha` must be a single finite number in [0, 0.5).", call. = FALSE)
  }
  retained <- dims[1] - floor(alpha * dims[1])
  if (retained < k) {
    stop("`alpha` trims too many observations for the requested `k`.", call. = FALSE)
  }
  if (!is_scalar_integerish(rank) || rank < 1L || rank > min(k * dims[2], dims[3])) {
    stop("`rank` must be an integer from 1 to min(k * T, p).", call. = FALSE)
  }
  if (!is.numeric(spar) || length(spar) != 1L || !is.finite(spar)) {
    stop("`spar` must be a single finite numeric smoothing parameter.", call. = FALSE)
  }
  if (!is_scalar_integerish(nstart) || nstart < 1L) {
    stop("`nstart` must be a positive integer.", call. = FALSE)
  }
  if (!is_scalar_integerish(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol_center) || length(tol_center) != 1L || !is.finite(tol_center) || tol_center < 0) {
    stop("`tol_center` must be a non-negative finite number.", call. = FALSE)
  }
  if (!is.null(min_cluster) && (!is_scalar_integerish(min_cluster) || min_cluster < 1L)) {
    stop("`min_cluster` must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is_scalar_integerish(min_segment) || min_segment < 1L) {
    stop("`min_segment` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(cutoff_bounds) || length(cutoff_bounds) != 2L ||
      any(!is.finite(cutoff_bounds)) || cutoff_bounds[1L] >= cutoff_bounds[2L]) {
    stop("`cutoff_bounds` must contain two increasing finite values.", call. = FALSE)
  }
  if (!is_scalar_integerish(seed)) {
    stop("`seed` must be a single integer.", call. = FALSE)
  }
  constant <- vapply(seq_len(dims[3L]), function(j) {
    length(unique(as.vector(x[, , j]))) <= 1L
  }, logical(1))
  if (any(constant)) {
    stop("`x` contains constant variables; remove or combine them before fitting.", call. = FALSE)
  }
  invisible(TRUE)
}

is_scalar_integerish <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && abs(x - round(x)) < sqrt(.Machine$double.eps)
}
