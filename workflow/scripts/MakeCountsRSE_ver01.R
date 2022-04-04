#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 1, 2022

#get arguments from command line input
args <- commandArgs(trailingOnly = TRUE)

countsTable <- args[1]
colData <- args[2]
RTWindowsBedFile <- args[3]
OutputFilename <- args[4]

# read counts
counts <- read.table("results/merged/test_counts.txt",header=T)

# read coldata
coldata <- read.csv("config/rif1Samples.csv")

# read RT window ranges
ranges.df <- read.table("resources/RTWindows_danRer11_noAlts_ver01.bed")

# convert dataframe to genomicRanges
ranges.gr <- GenomicRanges::makeGRangesFromDataFrame(ranges.df,seqnames.field="V1",start.field="V2",end.field="V3")

# make ranged summarized experiment object
rse <- SummarizedExperiment::SummarizedExperiment(
                    assays=list(counts=counts),
                    rowData=NULL, rowRanges=ranges.gr,
                    colData=coldata,
                    checkDimnames=TRUE)

# save ranged summarized experiment object as .rds file
saveRDS(rse,"results/test.rds")
