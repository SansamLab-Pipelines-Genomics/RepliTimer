#' Subset a SummarizedExperiment object then calculate medians of cpms for each subset.
#'
#' This function takes a SummarizedExperiment object containing counts per millions in the assays slot and subsets by Set to divide individual experiments. Then the medians for cpms are calculated for each subset and added to the rowData names of each object. This function is designed to be run after subsetting a SummarizedExperiment into G1 or S samples.
#'
#' @param A SummarizedExperiment object with counts per millions in the assays slot.
#' @return A SummarizedExperiment objecs containing medians of cpms in the rowData names slot.
#' @export

runMedianSteps <- function(subsetRTse){
  divideIntoSets <- splitSEbySet(subsetRTse)
  calcTheMedians <- getSubsetSEmedians(divideIntoSets)
  calcTheMedians
}
