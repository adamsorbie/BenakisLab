#!/usr/bin/env Rscript
# Author: Adam Sorbie 
# Date: 26/04/20
# Version 0.9.5

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
              help="trim left rev reads", metavar="character"),
  make_option(c("-R", "--trim_end"), type="integer", default=0, 
              help="trim right forward and reverse reads", metavar="character") 
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
trimright <- rep(opt$trim_end, 2)

fnFs <- sort(list.files(pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pattern="_R2_001.fastq.gz", full.names = TRUE))

# get sample names from forward reads 
sample.names <- sapply(strsplit(basename(fnFs), "_S1_L001_R1"), `[`, 1); sample.names

dir.create(opt$out)

if (!is.null(trimleft)){ 
  print("Reads untrimmed: trimming left portion of reads, before QC")
  trimFs <- file.path(path, "trimmed_primer", paste0(sample.names, "_trimmed_primer_R1_001.fastq.gz"))
  trimRs <- file.path(path, "trimmed_primer", paste0(sample.names, "_trimmed_primer_R2_001.fastq.gz"))
  primer_rem_out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs, trimLeft = trimleft, 
                                trimRight = trimright, multithread = opt$threads, 
                                compress = TRUE, truncQ = 0, rm.phix = FALSE)
  
  
  jpeg(paste(opt$out,"quality_profile_f.jpg", sep="/"))
  qc_F <- plotQualityProfile(trimFs, aggregate=TRUE)
  print(qc_F)
  dev.off()
  
  jpeg(paste(opt$out,"quality_profile_r.jpg", sep="/"))
  qc_R <- plotQualityProfile(trimRs, aggregate=TRUE)
  print(qc_R)
  dev.off()
  
} else {
  
  print("Reads are pre-trimmed, plotting Quality profile directly")
  
  jpeg(paste(opt$out,"quality_profile_f.jpg", sep="/"))
  qc_F <- plotQualityProfile(fnFs, aggregate=TRUE)
  print(qc_F)
  dev.off()
  
  jpeg(paste(opt$out,"quality_profile_r.jpg", sep="/"))
  qc_R <- plotQualityProfile(fnRs, aggregate=TRUE)
  print(qc_R)
  dev.off()
}



