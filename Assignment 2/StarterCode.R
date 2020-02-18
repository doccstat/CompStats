
#' Function that implements K-means algorithm
#'
#' @param X n times p data matrix (samples by features)
#' @param k number of clusters
#' @param S k starting cluster centers, if not supplied the algorithm chooses k random points out of n
#' @param numIter number of iterations, the default is 100
#'
#' @return
#' @export
#'
#' @examples
MyKmeans <- function(X, k, S = NULL, numIter = 100){
 
  # Check whether S is NULL or not. If NULL, initilize based on k random points from X. If not NULL, check for compatibility.
  
  # 
  
}