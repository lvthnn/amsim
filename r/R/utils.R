assert_probability <- function(x, len = NULL) {
  checkmate::assert_numeric(
    x,
    len = len,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    all.missing = FALSE
  )
}

assert_correlation_matrix <- function(matrix, cross = FALSE) {
  checkmate::assert_matrix(
    matrix,
    mode = "numeric",
    any.missing = FALSE,
    all.missing = FALSE
  )
  checkmate::assert(nrow(matrix) == ncol(matrix))
  checkmate::assert_numeric(matrix, lower = -1, upper = 1)
  if (!cross) checkmate::assert(all(diag(matrix) == diag(diag(nrow(matrix)))))
  checkmate::assert(all(eigen(matrix)$values >= 0))
}

ensure_length <- function(x, len) {
  checkmate::assert_count(len)
  if (length(x) == 1) x <- rep(x, len)
  else checkmate::assert(len(x) == len)
  return(x)
}
