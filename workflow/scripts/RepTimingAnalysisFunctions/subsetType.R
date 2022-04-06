#' Subset a SummarizedExperiment object by G1 or S
#'
#' This function subsets a SummarizedExperiment object based on Type from the colData container. The type will generally be either G1 or S.
#'
#' @param RTse A SummarizedExperiment object from the ReplicationTimingData package.
#' @param cycle Specifies the Type to subset from colData. This argument must be in quotes.
#' @return A SummarizedExperiment object containing only the colnames with the specified cell cycle Type.
#' @examples
#' counts <- matrix(rep(10000,100),
#' ncol=10,
#' dimnames = list(1:10,paste("Sample",1:10,sep="_")))
#' Granges <- GenomicRanges::GRanges(seqnames="chr1",IRanges::IRanges(start=1:10,end=(1:10)+1))
#' colData <- data.frame(set=c(1,1,1,1,1,2,2,2,2,2),
#'                       type=rep(c("G1","S"),5))
#' RTdata <- SummarizedExperiment::SummarizedExperiment(assays=list("counts"=counts),
#'                                                    rowRanges=Granges,
#'                                                    colData=colData)
#' onlyG1samples <- subsetType(RTdata,"G1")
#' SummarizedExperiment::colData(onlyG1samples)
#' @export

subsetType <- function(RTse, cycle){
  colData_DF <- SummarizedExperiment::colData(RTse)
  setsToGet <- which(grepl(cycle,colData_DF$type))
  SetType_se <- RTse[,setsToGet]
  SetType_se
}
