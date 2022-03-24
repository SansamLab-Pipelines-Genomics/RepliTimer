#!/usr/bin/env Rscript

# Chris Sansam
# version 04
# March 22, 2022

#get arguments from command line input
args <- commandArgs(trailingOnly = TRUE)

bedFilename <- args[1]
bamFileName <- args[2]
outputFileName <- args[3]

# import bed file
bedFile <- read.table(file=bedFilename,
                      header=F,
                      sep="\t")[,c(1,2,3)]
names(bedFile) <- c("seqnames","start","end")
bedFile <- GenomicRanges::makeGRangesFromDataFrame(bedFile,
                                                   ignore.strand = T,
                                                   seqnames.field=c("seqnames"),
                                                   start.field = "start",
                                                   end.field = "end")
counts <- bamsignals::bamCount(bamFileName,
                   bedFile,
                   verbose=FALSE,
                   paired.end="midpoint")
write.table(data.frame(as.data.frame(bedFile)[,c(1,2,3)],counts),
            file=outputFileName,
            col.names=F,
            row.names=F,
            quote=F,
            sep="\t")
