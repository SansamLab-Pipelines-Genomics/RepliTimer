#' Export Bedgraph with Specified Sample and Assay.
#'
#' This function takes a RangedSummarizedExperiment object with named data in the assays slot and creates a bedgraph file. An object is not returned. Instead a bedgraph file is saved.
#'
#' @param RTse A ranged SummarizedExperiment object with named data in the assays container.
#' @param SampleName Name of sample. Must match a samplename in the colData rownames. (in quotes)
#' @param Bedgraph_Filename Can include path. (in quotes)
#' @param Assay Assay named in the "assays" slot (in quotes)
#' @export

ExportBedgraphFromRTse <- function(RTse,SampleName,Bedgraph_Filename,Assay){
  assy <- assays(RTse)[[Assay]]
  assy <- assy[,grep(SampleName,colnames(assy))]
  bg_df <- cbind(data.frame(rowRanges(RTse))[,1:3],
                 assy)
  write.table(bg_df,
              file=Bedgraph_Filename,
              quote=F,sep="\t",
              row.names = F,
              col.names = F)
}
