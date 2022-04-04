#' Flatten list of RangedSummarizedExperiments
#'
#' This function converts of list of RangedSummarizedExperiments into a single.
#' Only the assay, rowRanges, and colData slots are flattened.
#'
#' @param se_list A list of SummarizedExperiments
#' @return A flattened SummarizedExperiment object with data in the assay, rowRanges, and colData slots.
#' @export

Flatten_SE_list <- function(se_list){
  SummarizedExperiment::SummarizedExperiment(assay=do.call("cbind",lapply(se_list,SummarizedExperiment::assay)),
                       rowRanges=SummarizedExperiment::rowRanges(se_list[[1]]),
                       colData=do.call("rbind",lapply(se_list,SummarizedExperiment::colData)))
}
