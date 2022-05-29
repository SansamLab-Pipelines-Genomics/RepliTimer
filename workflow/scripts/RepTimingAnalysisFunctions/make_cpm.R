#' Calculate counts per million
#'
#' This function calculates counts per million from a vector of timing window counts for the RepTiming pipeline.
#'
#' @param A vector of counts for sequencing reads in timing windows.
#' @return A vector with counts per million for each timing window.
#' @examples
#' counts_vector <- rep(10,100000)
#' make_cpm(counts_vector)
#' @export


make_cpm <- function(counts){
  cpmreads <- (counts / sum(counts,na.rm=T))*1000000
  return(cpmreads)
}
