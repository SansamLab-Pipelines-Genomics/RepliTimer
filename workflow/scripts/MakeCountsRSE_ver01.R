#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 1, 2022

#load libraries
library(SummarizedExperiment)
library(GenomicRanges)

#get arguments from command line input
args <- commandArgs(trailingOnly = TRUE)

countsTable <- args[1]
colData <- args[2]
RTWindowsBedFile <- args[3]
OutputFilename <- args[4]

# read counts
counts <- read.table(countsTable,header=T)

# read coldata
coldata <- read.csv(colData)

# read RT window ranges
ranges.df <- read.table(RTWindowsBedFile)

# convert dataframe to genomicRanges
ranges.gr <- GenomicRanges::makeGRangesFromDataFrame(ranges.df,seqnames.field="V1",start.field="V2",end.field="V3")

# make ranged summarized experiment object
rse <- SummarizedExperiment::SummarizedExperiment(
                    assays=list(counts=counts),
                    rowData=NULL, rowRanges=ranges.gr,
                    colData=coldata,
                    checkDimnames=TRUE)

# save ranged summarized experiment object as .rds file
saveRDS(rse,OutputFilename)
