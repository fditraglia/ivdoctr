collapse_3d_array <- function(myarray){
  out <- aperm(myarray, c(1, 3, 2))
  dim(out) <- c(dim(myarray)[1] * dim(myarray)[3], dim(myarray)[2])
  return(out)
}

#' Convert 3-d array to list of matrixes
#'
#' @param myArray A three-dimensional numeric array.
#' @return A list of numeric matrices.
#' @examples
#' M <- array(c(1, 1, 1, 1, 2, 2, 2, 2), c(2, 2, 2))
#' toList(M)
toList <- function(myArray){
  lapply(seq_len(dim(myArray)[3]), function(i) myArray[,,i,drop = TRUE])
}

