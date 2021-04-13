#!/usr/bin/env Rscript
library(dada2)
library(optparse)

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="dada2_qc", 
              help="output folder", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (path)", call.=FALSE)
}
path <- opt$path 

fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

out_dir <- paste(opt$path, opt$out, sep="/") 
dir.create(out_dir)

jpeg(paste(out_dir,"quality_profile_F.jpg", sep="/"), quality = 100)
plotQualityProfile(fnFs, n=1000000, aggregate = T)
dev.off()

jpeg(paste(out_dir,"quality_profile_R.jpg", sep="/"), quality=100)
plotQualityProfile(fnRs, n=1000000, aggregate = T)
dev.off()