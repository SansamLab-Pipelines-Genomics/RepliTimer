#' Convert counts to counts per million reads for a SummarizedExperiment object
#'
#' This function takes a SummarizedExperiment object as input and calculates counts per million reads for the counts assay. The counts per million reads will replace the counts in the assay element.
#'
#' @param RTse A SummarizedExperiment object with replication timing windows and sequencing read counts in the "counts" container in "assays".
#' @return A SummarizedExperient object with cpms as the assay element.
#' @export
#' @examples
#' counts <- matrix(rep(10000,100),
#' ncol=10,
#' dimnames = list(1:10,paste("Sample",1:10,sep="_")))
#' Granges <- GenomicRanges::GRanges(seqnames="chr1",IRanges::IRanges(start=1:10,end=(1:10)+1))
#' colData <- data.frame(set=c(1,1,1,1,1,2,2,2,2,2),
#'                       type=rep(c("G1","S"),5))
#' RTse <- SummarizedExperiment::SummarizedExperiment(assays=list("counts"=counts),
#'                                                    rowRanges=Granges,
#'                                                    colData=colData)
#' RTse_cpms <- ConvertAllCountsToCPMs(RTse)
#' SummarizedExperiment::assays(RTse_cpms)$cpms

ConvertAllCountsToCPMs <- function(RTse){
  counts_mx <- SummarizedExperiment::assays(RTse,withDimnames=FALSE)$counts
  allcpms <- apply(counts_mx, 2, make_cpm)
  SummarizedExperiment::assays(RTse) <- list("cpms"=allcpms)
  RTse
}
