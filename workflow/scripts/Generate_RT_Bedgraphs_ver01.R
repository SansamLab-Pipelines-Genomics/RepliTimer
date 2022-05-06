#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 5, 2022

library(SummarizedExperiment)
library(GenomicRanges)

# load functions
sapply(list.files(pattern="[.]R$", path="workflow/scripts/RepTimingAnalysisFunctions", full.names=TRUE), source)

#get arguments from command line input
args <- commandArgs(trailingOnly = TRUE)

rseFile <- args[1]
assay <- args[2]
assay2 <- paste(trimws(basename(assay)),collapse=" ")

rse <- readRDS(rseFile)

SmplNmes <- row.names(colData(rse))

for (SmplNme in SmplNmes) {
    ExportBedgraphFromRTse(
        RTse=rse,
        SampleName=SmplNme,
        Assay=assay2,
        Bedgraph_Filename=paste(
                         assay,
                         "/",
                         SmplNme,"_",assay2,
                         ".bedgraph",
                         sep="") 
                       )
}
