#' Extract MTRClust diagnostics
#'
#' Returns compact diagnostic tables from a fitted `"mtrclust"` object.
#'
#' @param fit An object returned by [mtrclust()].
#' @param type `"counts"` for multilevel counts or `"units"` for observation
#'   outlier and retained-unit membership.
#' @return A data frame.
#' @export
mtrclust_diagnostics <- function(fit, type = c("counts", "units")) {
  if (!inherits(fit, "mtrclust")) {
    stop("`fit` must be an object returned by mtrclust().", call. = FALSE)
  }
  type <- match.arg(type)
  if (type == "units") {
    return(data.frame(
      unit = seq_along(fit$cluster),
      cluster = fit$cluster,
      nearest_cluster = fit$nearest_cluster,
      retained = as.logical(fit$keep),
      observation_outlier = !as.logical(fit$keep)
    ))
  }
  raw <- fit$raw_flags %||% list()
  data.frame(
    level = c(
      "observation_outlier",
      "global_sensor_time",
      "temporal_point",
      "multivariate_segment",
      "sensor_segment",
      "isolated_cell",
      "raw_cell_detection",
      "raw_time_detection"
    ),
    count = c(
      length(fit$observation_outliers),
      nrow_or_zero(fit$global_sensor_time_outliers),
      nrow_or_zero(fit$time_point_outliers),
      nrow_or_zero(fit$multivariate_segments),
      nrow_or_zero(fit$sensor_segments),
      nrow_or_zero(fit$isolated_cell_outliers),
      nrow_or_zero(raw$cell_outliers),
      nrow_or_zero(raw$time_point_outliers)
    )
  )
}

nrow_or_zero <- function(x) {
  if (is.null(x)) return(0L)
  d <- dim(x)
  if (is.null(d)) return(length(x))
  nrow(x)
}
