#' Print an MTRClust fit
#' @param x An object returned by [mtrclust()].
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.mtrclust <- function(x, ...) {
  d <- x$input_dim %||% infer_fit_dim(x)
  cat("MTRClust fit\n")
  cat("Data:", d[1L], "units x", d[2L], "time points x", d[3L], "variables\n")
  cat("K:", x$settings$k, " alpha:", x$settings$alpha, "\n", sep = "")
  cat("Retained cluster sizes:", paste(x$retained_cluster_sizes, collapse = ", "), "\n")
  cat("Observation outliers:", length(x$observation_outliers), "\n")
  cat("Objective:", format(signif(x$objective, 6L)), "\n")
  cat("Termination:", x$termination, " | converged:", x$converged,
      " | iterations:", x$iterations, " | rolled back:", x$rolled_back, "\n")
  invisible(x)
}

#' Summarize an MTRClust fit
#' @param object An object returned by [mtrclust()].
#' @param ... Unused.
#' @return An object of class `"summary.mtrclust"`.
#' @export
summary.mtrclust <- function(object, ...) {
  out <- list(
    dimensions = object$input_dim %||% infer_fit_dim(object),
    settings = object$settings,
    cluster_sizes = object$retained_cluster_sizes,
    trimming = data.frame(
      retained = sum(object$keep),
      trimmed = sum(!object$keep),
      trim_rate = mean(!object$keep)
    ),
    termination = data.frame(
      termination = object$termination,
      converged = object$converged,
      cycle_detected = object$cycle_detected,
      rolled_back = object$rolled_back,
      iterations = object$iterations,
      best_iteration = object$best_iteration,
      objective = object$objective
    ),
    diagnostics = mtrclust_diagnostics(object, type = "counts"),
    candidate_summary = object$candidate_summary,
    runtime_seconds = object$runtime_seconds %||% NA_real_
  )
  class(out) <- "summary.mtrclust"
  out
}

#' Print an MTRClust summary
#' @param x An object returned by [summary.mtrclust()].
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.summary.mtrclust <- function(x, ...) {
  cat("MTRClust summary\n")
  cat("Dimensions:", paste(x$dimensions, collapse = " x "), "\n")
  cat("Cluster sizes:", paste(x$cluster_sizes, collapse = ", "), "\n")
  cat("Trimming:", x$trimming$trimmed, "of",
      x$trimming$retained + x$trimming$trimmed,
      sprintf("(%.3f)", x$trimming$trim_rate), "\n")
  cat("Termination:", x$termination$termination,
      "| converged:", x$termination$converged,
      "| rolled back:", x$termination$rolled_back, "\n")
  cat("Diagnostics:\n")
  print(x$diagnostics, row.names = FALSE)
  invisible(x)
}

#' Plot an MTRClust fit
#' @param x An object returned by [mtrclust()].
#' @param type Plot type: `"objective"`, `"sizes"`, `"trajectories"`,
#'   `"outliers"`, or `"diagnostics"`.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
plot.mtrclust <- function(x, type = c("objective", "sizes", "trajectories", "outliers", "diagnostics"), ...) {
  type <- match.arg(type)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  if (type == "objective") {
    y <- x$objective_history
    graphics::plot(seq_along(y), y, type = "b", xlab = "Iteration",
                   ylab = "Objective", main = "MTRClust objective history")
  } else if (type == "sizes") {
    graphics::barplot(x$retained_cluster_sizes, names.arg = seq_along(x$retained_cluster_sizes),
                      xlab = "Cluster", ylab = "Retained units",
                      main = "Retained cluster sizes")
  } else if (type == "trajectories") {
    centers <- x$centers
    graphics::matplot(t(centers[, , 1L]), type = "l", lty = 1, xlab = "Time",
                      ylab = "Scaled center", main = "Cluster trajectories, variable 1")
  } else if (type == "outliers") {
    counts <- c(retained = sum(x$keep), observation_outliers = sum(!x$keep))
    graphics::barplot(counts, ylab = "Units", main = "Observation-level trimming")
  } else {
    d <- mtrclust_diagnostics(x, type = "counts")
    graphics::barplot(d$count, names.arg = d$level, las = 2,
                      ylab = "Count", main = "Multilevel diagnostic counts")
  }
  invisible(x)
}

infer_fit_dim <- function(x) {
  n <- length(x$cluster)
  centers <- x$centers
  c(n, dim(centers)[2L], dim(centers)[3L])
}

`%||%` <- function(x, y) if (is.null(x)) y else x
