#' Smooth values in rse assay slot.
#'
#' This function takes a Ranged SummarizedExperiment object with a matrix of replication timing values in the assay slot and returns an rse with the values smoothed.
#'
#' @param rse A ranged summarized experiment object with RT values in the assay slot.
#' @return A ranged summarized experiment object with smoothed RT values in the assay slot.
#' @export

SmoothRseAssay <- function(rse){
  #rse <- rse[complete.cases(SummarizedExperiment::assay(rse)),]
  rseByChrm <- splitSEbyChr(rse)
  #RT.ZScores<-split(as.data.frame(Z.score.df),f=Z.score.df$seqnames)
  #matlab csaps smoothing variable (p)
  p=0.0000000000000001
  #p=0.0000000000000001
  spar_value=(1-p)/p #conversion to R variable
  smoothed.list<-lapply(rseByChrm,function(rsechrm){
    #rsechrm <- rsechrm[complete.cases(SummarizedExperiment::assay(rsechrm)),]
    rawRTvalues <- SummarizedExperiment::assay(rsechrm)
    midpoints <- matrixStats::rowMeans2(cbind(SummarizedExperiment::start(rsechrm),
                                              SummarizedExperiment::end(rsechrm)))
    midpoints <- round(midpoints, digits=0)
    smoothed<-apply(rawRTvalues,2,function(raw_values){
      smoothed_values<-c(predict(pspline::smooth.Pspline(
        x = midpoints,
        y = raw_values,
        norder = 2,
        method = 1,
        #spar = spar_value),Zscores$midpoints))
        spar = spar_value),midpoints)) ##
    })
    return(smoothed)
  })
  SummarizedExperiment::assay(rse,withDimnames=FALSE)<-(do.call("rbind",smoothed.list))
  return(rse)
}
