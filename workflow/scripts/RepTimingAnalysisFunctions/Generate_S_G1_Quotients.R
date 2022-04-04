#' Divide each S phase cpms by G1 median cpms of same set
#'
#' This function subsets a SummarizedExperiment object based on Type from the colData container. The type will generally be either G1 or S.
#'
#' @param S_RTse_cpms A SummarizedExperiment object with all S-phase cpms in the assay container.
#' @param Median_RTse_cpms_list A list of SummarizedExperiment objects for each set. G1 medians are in the rowData slot.
#' @return A list of SummarizedExperiment objects with S/G1 quotients in the assay slot.
#' @export

Generate_S_G1_Quotients <- function(S_RTse_cpms, Median_RTse_cpms_list){
  Ssets <- splitSEbySet(S_RTse_cpms)
  sets <- intersect(names(Ssets),names(Median_RTse_cpms_list))
  names(sets) <- sets
  quotients_se <- lapply(sets,function(set){
    medianG1cpms <- SummarizedExperiment::rowData(Median_RTse_cpms_list[[set]])$MedianG1CPMs
    quotients <- apply(SummarizedExperiment::assay(Ssets[[set]]),2,function(G1cpms){G1cpms/medianG1cpms})
    Sset <- Ssets[[set]]
    SummarizedExperiment::assay(Sset) <- quotients
    names(SummarizedExperiment::assays(Sset)) <- "quotients"
    Sset
  })
  quotients_se
}
