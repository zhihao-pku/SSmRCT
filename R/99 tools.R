#' Title
#'
#' @param A
#' @param B
#'
#' @return
#' @export
#'
#' @examples
combine_matrix <- function(A, B) {
  dim <- nrow(A) * 2
  combined_matrix <- matrix(0, nrow = dim, ncol = dim)
  for (o in 1:nrow(A)) {
    combined_matrix[(2 * o - 1), seq(1, dim, by = 2)] <- A[o, ]
    combined_matrix[(2 * o), seq(1, dim, by = 2)] <- B[o, ]
    combined_matrix[(2 * o - 1), seq(2, dim, by = 2)] <- A[o, ]
    combined_matrix[(2 * o), seq(2, dim, by = 2)] <- B[o, ]
  }
  return(combined_matrix)
}

#' Title
#'
#' @param a
#' @param b
#'
#' @return
#' @export
#'
#' @examples
combine_vector <- function(a, b) {
  combined_vector <- numeric(length(a) + length(b))
  combined_vector[seq(1, length(combined_vector), by = 2)] <- a
  combined_vector[seq(2, length(combined_vector), by = 2)] <- b
  return(combined_vector)
}
