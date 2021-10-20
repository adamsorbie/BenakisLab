#!/usr/bin/env Rscript
# Author: Adam Sorbie 
# Date: 13/09/21
# Version 0.9.6

library(dada2)
library(optparse)
library(parallel)

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="dada2_qc", 
              help="output folder", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
# Functions

RIGHT <- function(x,n){
  substring(x,nchar(x)-n+1)
}

# check all paths have trailing forward slash 
if (RIGHT(opt$path, 1) != "/") {
  opt$path <- paste0(opt$path, "/")
}
if (RIGHT(opt$out, 1) != "/") {
  opt$out <- paste0(opt$out, "/")
}

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least three arguments must be supplied (path, trimF, trimR)", 
       call.=FALSE)
}

path <- opt$path 
setwd(opt$path)


fnFs <- sort(list.files(pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pattern="_R2_001.fastq.gz", full.names = TRUE))

# get sample names from forward reads 
sample.names <- sapply(strsplit(basename(fnFs), "_S1_L001_R1"), `[`, 1); sample.names

dir.create(opt$out)

qc_F <- plotQualityProfile(fnFs, aggregate=TRUE)
ggsave(paste(opt$out,"quality_profile_f.jpg", sep="/"), qc_F, device = "jpeg", dpi = 200)

qc_R <- plotQualityProfile(fnRs, aggregate=TRUE)
ggsave(paste(opt$out,"quality_profile_r.jpg", sep="/"), qc_R, device = "jpeg", dpi = 200)




