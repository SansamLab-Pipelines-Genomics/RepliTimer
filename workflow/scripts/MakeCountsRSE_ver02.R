#!/usr/bin/env Rscript

# Chris Sansam
# version 01
# April 28, 2022

# load libraries
library(SummarizedExperiment)
library(GenomicRanges)

# get arguments from Snakemake S4 object
countsBedgraphs <- snakemake@params[["counts_bedgraphs_list"]]
print("this is the bedgraphs object before stringsplit")
countsBedgraphs
countsBedgraphs <- strsplit(snakemake@params[["counts_bedgraphs_list"]], split = ",")[[1]]
print("this is the bedgraphs vector after stringsplit")
countsBedgraphs
sampleNames <- strsplit(snakemake@params[["sample_list"]], split = ",")[[1]]
colData <- snakemake@params[["samples_table"]]
RTWindowsBedFile <- snakemake@params[["RT_windows"]]
OutputFilename <- snakemake@output[["rse_counts"]]

# read 4th column of each bedgraph and then cbind to a dataframe
counts <- do.call("cbind",lapply(countsBedgraphs,function(filename){read.table(filename)[,4]}))
names(counts) <- sampleNames

# read coldata
coldata <- read.csv(colData)

# read RT window ranges
ranges.df <- read.table(countsBedgraphs[1])

# convert dataframe to genomicRanges
ranges.gr <- GenomicRanges::makeGRangesFromDataFrame(ranges.df,seqnames.field="V1",start.field="V2",end.field="V3")

print("here is the head of the ranges.gr")
head(ranges.gr)

# make ranged summarized experiment object
rse <- SummarizedExperiment::SummarizedExperiment(
                    assays=list(counts=counts),
                    rowData=NULL, rowRanges=ranges.gr,
                    colData=coldata,
                    checkDimnames=TRUE)
print("Here is some info on the rse made")
rse
print("Here are the dimensions of the assays")
dim(assays(rse)$counts)
# save ranged summarized experiment object as .rds file
saveRDS(rse,OutputFilename)
