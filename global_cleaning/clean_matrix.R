preprocess_matrix <- function(matrix) {
    matrix <- as.matrix(matrix)
    # Remove rows  with NA or Inf
    matrix <- matrix[complete.cases(matrix), , drop = FALSE]
    matrix <- matrix[rowSums(is.infinite(matrix)) == 0, , drop = FALSE]
    # Remove rows with all zeros
    matrix <- matrix[rowSums(matrix) > 0, ]  # Remove rows with all zeros
    # Remove rows with zero variance
    matrix <- matrix[apply(matrix, 1, var) > 0, , drop = FALSE]
    return(matrix)
  }