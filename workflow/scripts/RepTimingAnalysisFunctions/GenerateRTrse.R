#' Convert an RTrse with counts to RTrse with quotients and normalized and scaled values.
#'
#' This function takes a RangedSummarizedExperiment object with read counts in the assay slot and returns one with processed RT numbers.
#'
#' @param RTse A ranged SummarizedExperiment object with read counts in the assay container.
#' @param AllowGaps If True, smoothing over gaps longer than 200kb will occur
#' @return Returns an RTrse with quotients ("Quotients"), smoothed ("Smoothed"), scaled ("ZScores") values in assays. For the smoothed median quotients, overall ZScores are calculated and these values are placed in the "mcols" container.
#' @export

GenerateRTrse <- function(RTse,AllowGaps=TRUE){
  # sort the rse
    RTse <- SummarizedExperiment::sort(RTse)
  # convert counts to cpms
    RTse_cpms <- ConvertAllCountsToCPMs(RTse)
  # subset into G1 and S
    onlyG1samples <- subsetType(RTse_cpms,"G1")
    onlySsamples <- subsetType(RTse_cpms,"S")
  # Calculate median cpms for each set
    G1samplesMedians <- runMedianSteps(onlyG1samples)
  # Divide each S phase cpms by G1 median cpms of same set
    Quotients_se_list <- Generate_S_G1_Quotients(onlySsamples, G1samplesMedians)
    Quotients_se <- Flatten_SE_list(Quotients_se_list)
    Quotients_se <- Quotients_se[complete.cases(SummarizedExperiment::assay(Quotients_se,withDimnames=FALSE)),]
  # calculate Smoothed Quotients
    if(isTRUE(AllowGaps)){
      QuotientsSmoothed_se <- SmoothRseAssay(Quotients_se)
    }else{
      QuotientsSmoothed_se <- SmoothRseAssayToGaps(Quotients_se)
    }
  # calculate Smoothed Quotient Z Scores
    QuotientsSmoothedScaled_se <- ZScoreRseAssay(QuotientsSmoothed_se)
  # rebuild rse with all values in assays
    rse <- QuotientsSmoothedScaled_se
    SummarizedExperiment::assay(rse,withDimnames=FALSE) <- NULL
    SummarizedExperiment::assays(rse,withDimnames=FALSE)$ZScores <- SummarizedExperiment::assay(QuotientsSmoothedScaled_se,withDimnames=FALSE)
    SummarizedExperiment::assays(rse,withDimnames=FALSE)$Smoothed <- SummarizedExperiment::assay(QuotientsSmoothed_se,withDimnames=FALSE)
    matches <- GenomicRanges::findOverlaps(SummarizedExperiment::rowRanges(rse),
                                    SummarizedExperiment::rowRanges(Quotients_se),
                                    type="equal")
    Quotients_se <- Quotients_se[subjectHits(matches),]
    SummarizedExperiment::assays(rse,withDimnames=FALSE)$Quotients <- SummarizedExperiment::assay(Quotients_se,withDimnames=FALSE)

  if(is.element('Replicate_ID', names(SummarizedExperiment::colData(Quotients_se)))){
    MedianQuotients_se <- CalculateReplicateQuotMedians(Quotients_se)
    MedianQuotients_se <- SummarizedExperiment::sort(MedianQuotients_se)
    ##
    if(isTRUE(AllowGaps)){
      MedianQuotientsSmoothed_se <- SmoothRseAssay(MedianQuotients_se)
    }else{
      MedianQuotientsSmoothed_se <- SmoothRseAssay(MedianQuotients_se)
    }
    ##
    MedianQuotientsSmoothedScaled_se <- ZScoreRseAssay(MedianQuotientsSmoothed_se)
    rse2 <- MedianQuotientsSmoothedScaled_se
    rse3 <- SummarizedExperiment::SummarizedExperiment(colData = rbind(SummarizedExperiment::colData(rse2),
                                                                       SummarizedExperiment::colData(rse)),
                                                       rowRanges = SummarizedExperiment::rowRanges(rse))
    ZScores <- cbind(SummarizedExperiment::assay(MedianQuotientsSmoothedScaled_se,withDimnames=FALSE),
                     SummarizedExperiment::assay(QuotientsSmoothedScaled_se,withDimnames=FALSE))
    Smoothed <- cbind(SummarizedExperiment::assay(MedianQuotientsSmoothed_se,withDimnames=FALSE),
                      SummarizedExperiment::assay(QuotientsSmoothed_se,withDimnames=FALSE))
    Quotients <- cbind(SummarizedExperiment::assay(MedianQuotients_se,withDimnames=FALSE),
                       SummarizedExperiment::assay(Quotients_se,withDimnames=FALSE))
    #ZScoresAll <- matrix(scale(c(SummarizedExperiment::assay(MedianQuotientsSmoothedScaled_se))),
                               #nrow=nrow(SummarizedExperiment::assay(MedianQuotientsSmoothedScaled_se)))
    SummarizedExperiment::assays(rse3,withDimnames=FALSE)$ZScores <- ZScores
    SummarizedExperiment::assays(rse3,withDimnames=FALSE)$Smoothed <- Smoothed
    SummarizedExperiment::assays(rse3,withDimnames=FALSE)$Quotients <- Quotients
    #SummarizedExperiment::mcols(rse3) <- ZScoresAll
  }
  return(rse3)
}