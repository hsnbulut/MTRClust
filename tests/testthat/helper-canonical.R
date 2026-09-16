canonical_partition <- function(z) {
  z <- as.integer(z)
  out <- integer(length(z))
  if (any(z == 0L)) out[z == 0L] <- 0L
  labs <- setdiff(unique(z), 0L)
  if (!length(labs)) return(out)
  first <- vapply(labs, function(a) min(which(z == a)), integer(1))
  labs <- labs[order(first)]
  for (i in seq_along(labs)) out[z == labs[i]] <- i
  out
}

nrow_test <- function(x) {
  if (is.null(x)) return(0L)
  d <- dim(x)
  if (is.null(d)) length(x) else nrow(x)
}
