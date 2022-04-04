#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 4, 2022

#load libraries
library(SummarizedExperiment)
library(GenomicRanges)

#get arguments from command line input
args <- commandArgs(trailingOnly = TRUE)

#set positional arguments to objects
rdsFilename <- args[1]
OutputFilename <- args[2]

# Open rse
RTse <- readRDS(rdsFilename)

# Replace the ranged coordinates to window-centered coordindates.
ranges2 <- SummarizedExperiment::ranges(RTse)
centers <- rowMeans(cbind(GenomicRanges::start(ranges2),GenomicRanges::end(ranges2)))
GenomicRanges::start(ranges2) <- centers
GenomicRanges::end(ranges2) <- centers
SummarizedExperiment::ranges(RTse) <- ranges2

# load functions for calculating quotients and smooting and scaling RT values
sapply(list.files(pattern="[.]R$", path="workflow/scripts/RepTimingAnalysisFunctions/", full.names=TRUE), source)

# calculate quotients and then smooth and scale
Processed_RTse <- GenerateRTrse(RTse,AllowGaps=FALSE)

# save ranged summarized experiment object as .rds file
saveRDS(Processed_RTse,OutputFilename)
