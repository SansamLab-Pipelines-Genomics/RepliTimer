#' Calculates Medians of counts per millions.
#'
#' This function takes a SummarizedExperiment object that has counts per millions in the assays slot and calculates the median values. This function is designed to accessed through the runMedianSteps function.
#'
#' @param A SummarizedExperiment object with counts per millions in the assays slot.
#' @return A SummarizedExperiment object with medians of cpms in the rowData names slot.
#' @export

CalcMedians <- function(oneG1se){
  getCpms <- SummarizedExperiment::assays(oneG1se,withDimnames=FALSE)$cpms
  MedianG1CPMs <- matrixStats::rowMedians(getCpms, na.rm = TRUE)
  SummarizedExperiment::mcols(oneG1se) <- data.frame(MedianG1CPMs)
  oneG1se
}
