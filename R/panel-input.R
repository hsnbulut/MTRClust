#' Convert long-format data to an MTRClust panel array
#'
#' @param data A data frame in long format.
#' @param id Name of the unit identifier column.
#' @param time Name of the time index column.
#' @param variables Character vector of numeric variable columns.
#' @return A numeric array with dimensions unit by time by variable. Unit and
#'   time labels are stored as dimnames.
#' @examples
#' d <- data.frame(id = rep(1:2, each = 3), time = rep(1:3, 2),
#'                 y1 = rnorm(6), y2 = rnorm(6))
#' as_mtr_panel(d, id = "id", time = "time", variables = c("y1", "y2"))
#' @export
as_mtr_panel <- function(data, id, time, variables) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  required <- c(id, time, variables)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!length(variables)) stop("`variables` must contain at least one column.", call. = FALSE)
  if (!all(vapply(data[variables], is.numeric, logical(1)))) {
    stop("All `variables` columns must be numeric.", call. = FALSE)
  }
  units <- unique(data[[id]])
  times <- sort(unique(data[[time]]))
  out <- array(NA_real_, c(length(units), length(times), length(variables)),
               dimnames = list(as.character(units), as.character(times), variables))
  key <- paste(data[[id]], data[[time]], sep = "\r")
  if (any(duplicated(key))) {
    stop("`data` contains duplicate id-time rows.", call. = FALSE)
  }
  unit_index <- match(data[[id]], units)
  time_index <- match(data[[time]], times)
  for (j in seq_along(variables)) {
    out[cbind(unit_index, time_index, j)] <- data[[variables[j]]]
  }
  if (anyNA(out)) {
    stop("`data` does not form a complete id by time panel.", call. = FALSE)
  }
  out
}
