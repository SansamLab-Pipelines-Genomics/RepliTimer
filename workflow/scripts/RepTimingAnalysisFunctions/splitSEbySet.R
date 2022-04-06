#' Separate a SummarizedExperiment object by Sets and produce a list.
#'
#' This function takes a SummarizedExperiment object as input and returns separate objects based on Set number so that grouped samples can be further analyzed. This function is designed to accessed through the runMedianSteps function.
#'
#' @param RTse is a SummarizedExperiment object, typically an object that has been subset for G1 samples.
#' @return Returns a list containing a SummarizedExperiment object for each Set number.
#' @export

splitSEbySet <- function(RTse) {
  se_list <- lapply(unique(SummarizedExperiment::colData(RTse)$set), function(set) {
    RTse[, which(SummarizedExperiment::colData(RTse)$set == set)]
  })
  names(se_list) <- unique(SummarizedExperiment::colData(RTse)$set)
  return(se_list)
}
