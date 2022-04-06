#' Calculate median of quotients for each set of replicates.
#'
#' Calculate median of quotients for each set of replicates.
#'
#' @param A list of SummarizedExperiment object with quotients in assay container and replicate groups marked by "Replicate_ID" in colData.
#' @return A SummarizedExperiment object containing medians of quotients for replicate groups.
#' @export

CalculateReplicateQuotMedians <- function(Quotients_se){
  Quotients_sebyRepID <- splitSEbyReplicateID(Quotients_se)
  RowMeds <- lapply(Quotients_sebyRepID,function(rse){matrixStats::rowMedians(SummarizedExperiment::assay(rse,withDimnames=FALSE))})
  RowMeds <- do.call("cbind",RowMeds)
  row.names(RowMeds) <- row.names(Quotients_se)
  colData <- data.frame(row.names = colnames(RowMeds))
  colData[,names(SummarizedExperiment::colData(Quotients_se))] <- NA
  colData <- lapply(Quotients_sebyRepID,function(rse){
    clDta <- SummarizedExperiment::colData(rse)
    apply(clDta,2,function(clmn){
      if(length(unique(clmn)) > 1){
        return(NA)
      }else{
        return(clmn[1])
      }})
  })
  colData <- do.call("rbind",colData)
  SummarizedExperiment::SummarizedExperiment(assay=RowMeds,
                              rowRanges=SummarizedExperiment::rowRanges(Quotients_se), colData=colData)
}
