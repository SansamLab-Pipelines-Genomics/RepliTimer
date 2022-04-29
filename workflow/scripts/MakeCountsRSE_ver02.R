#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 28, 2022

# load libraries
library(SummarizedExperiment)
library(GenomicRanges)

# get arguments from Snakemake S4 object
countsBedgraphs <- strsplit(snakemake@input[["counts_bedgraphs_list"]], split = " ")
sampleNames <- strsplit(snakemake@params[["sample_list"]], split = ",")
colData <- snakemake@params[["samples_table"]]
RTWindowsBedFile <- snakemake@params[["RT_windows"]]
OutputFilename <- snakemake@output[["rse_counts"]]

# read 4th column of each bedgraph and then cbind to a dataframe
counts <- do.call("cbind",lapply(countsBedgraphs,function(filename){read.table(filename)[,4]}))
names(counts) <- sampleNames

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
