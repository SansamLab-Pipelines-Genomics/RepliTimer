#' Separate a SummarizedExperiment object by Replicate_IDs and produce a list.
#'
#' This function takes a SummarizedExperiment object as input and returns separate objects based on Replicate_ID so that replicate samples can be further analyzed.
#'
#' @param RTse is a SummarizedExperiment object, typically an object that has been subset for G1 samples.
#' @return Returns a list containing a SummarizedExperiment object for each Replicate_ID
#' @export

splitSEbyReplicateID <- function(RTse) {
  se_list <- lapply(unique(SummarizedExperiment::colData(RTse)$Replicate_ID), function(Replicate_ID) {
    RTse[, which(SummarizedExperiment::colData(RTse)$Replicate_ID == Replicate_ID)]
  })
  names(se_list) <- unique(SummarizedExperiment::colData(RTse)$Replicate_ID)
  return(se_list)
}
