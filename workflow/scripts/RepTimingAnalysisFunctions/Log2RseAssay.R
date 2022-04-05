#' Scale values in rse assay slot.
#'
#' This function takes a Ranged SummarizedExperiment object with a matrix of replication timing values in the assay slot and returns an rse with the log2 values calculated
#'
#' @param rse A ranged summarized experiment object with RT values in the assay slot.
#' @return A ranged summarized experiment object with log2 ratio RT values in the assay slot.
#' @export

Log2RseAssay <- function(rse){
  SummarizedExperiment::assay(rse,withDimnames=FALSE) <- apply(SummarizedExperiment::assay(rse),2,log2)
  return(rse)
}
