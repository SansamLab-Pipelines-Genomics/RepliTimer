#' Separate a SummarizedExperiment object by chromosome and produce a list.
#'
#' This function takes a SummarizedExperiment object as input and returns separate objects based on chromosome.
#'
#' @param RTse is a SummarizedExperiment object.
#' @return Returns a list containing a SummarizedExperiment object for each chromsome.
#' @export

splitSEbyChr <- function(RTse) {
  se_list <- lapply(unique(as.vector(SummarizedExperiment::seqnames(RTse))), function(chrm) {
    RTse[which(as.vector(SummarizedExperiment::seqnames(RTse)) == chrm),]
  })
  names(se_list) <- unique(as.vector(SummarizedExperiment::seqnames(RTse)))
  return(se_list)
}
