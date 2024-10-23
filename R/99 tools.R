#' Title
#'
#' @param A
#' @param B
#'
#' @return
#' @export
#'
#' @examples
#' A <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, ncol = 2, byrow = TRUE)
#' combine(A, A)
#'
#' A <- c(1, 2, 3)
#' B <- c(4, 5, 6)
#' combine(A, B)
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
