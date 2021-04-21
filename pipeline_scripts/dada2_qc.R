#!/usr/bin/env Rscript
# Author: Adam Sorbie 
# Date: 19/04/20
# Version 0.9.0

library(dada2)
library(optparse)
library(parallel)

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="dada2_qc", 
              help="output folder", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=detectCores(), 
              help="number of threads", metavar="character"),
  make_option(c("-f", "--trimF"), type="integer", default=NULL, 
              help="trim left fwd reads", metavar="character"),
  make_option(c("-r", "--trimR"), type="integer", default=NULL, 
              help="trim left rev reads", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least three arguments must be supplied (path, trimF, trimR)", 
       call.=FALSE)
}

path <- opt$path 
setwd(opt$path)

trimleft <- c(opt$trimF, opt$trimR)

# invoke system command ls and cut to get sample names from fastq filenames
system('ls *R1_001.fastq.gz | cut -f1 -d"-" > samples')

# store sample names as variable
samples <- scan("samples", what = "character")

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

trimFs <- file.path(path, "trimmed_primer", paste0(samples, "_trimmed_primer_R1_001.fastq.gz"))
trimRs <- file.path(path, "trimmed_primer", paste0(samples, "_trimmed_primer_R2_001.fastq.gz"))

primer_rem_out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs, trimLeft = trimleft, 
                                multithread = opt$threads, compress = TRUE, 
                                truncQ = 0, rm.phix = FALSE)
  
dir.create(opt$out)

jpeg(paste(opt$out,"quality_profile_F.jpg", sep="/"), quality = 100)
plotQualityProfile(trimFs, n=1e6, aggregate = T)
dev.off()

jpeg(paste(opt$out,"quality_profile_R.jpg", sep="/"), quality=100)
plotQualityProfile(trimRs, n=1e6, aggregate = T)
dev.off()
