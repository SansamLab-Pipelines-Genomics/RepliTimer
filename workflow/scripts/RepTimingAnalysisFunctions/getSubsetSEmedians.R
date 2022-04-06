#' Calculate medians of cpms for a list of SummarizedExperiment objects.
#'
#' This function takes a list of SummarizedExperiment objects created from the splitSEbySet function and calculates the medians of cpms for each object. This function is designed to accessed through the runMedianSteps function.
#'
#' @param A list of SummarizedExperiment objects resulting from the function splitSEbySet
#' @return A list of SummarizedExperiment objects containing medians of cpms in the rowData names slot.
#' @export

getSubsetSEmedians <- function(listfromG1se){
  calcG1medians <- lapply(listfromG1se, CalcMedians)
  calcG1medians
}
