#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 5, 2022

#get arguments from command line input
args <- commandArgs(trailingOnly = TRUE)

rse <- args[1]
assay <- args[2]
outputFileName <- args[3]
SmplNmes <- row.names(colData(rse))

for (SmplNme in SmplNmes) {
    ExportBedgraphFromRTse(
        RTse=rse,
        SampleName=SmplNme,
        Assay=assay,
        Bedgraph_Filename=paste(
                         assay,
                         "/",
                         SmplNme,"_",assay,
                         ".bedgraph",
                         sep="") 
                       )
}

files2zip <- dir(assay, full.names = TRUE)
zip(zipfile = outputFileName, files = files2zip)
