#' Combine vectors or matrices
#'
#' Combine two vectors to one vector or two matrices to one matrix
#'
#' @name combine
#' @param A A vector or matrix
#' @param B A vector or matrix
#'
#' @return A vector or matrix
#'
#' @examples
#' combine(c(1, 2, 3), c(4, 5, 6))
#'
#' combine(
#'   matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE),
#'   matrix(c(5, 6, 7, 8), nrow = 2, byrow = TRUE)
#' )
combine <- function(A, B) {
  if (is.matrix(A)) {
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
  if (is.vector(A)) {
    combined_vector <- numeric(length(A) + length(B))
    combined_vector[seq(1, length(combined_vector), by = 2)] <- A
    combined_vector[seq(2, length(combined_vector), by = 2)] <- B
    return(combined_vector)
  }
}
