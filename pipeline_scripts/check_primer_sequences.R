#!/usr/bin/ Rscript

#' Author: Adam Sorbie 
#' Date: 16/04/21
#' Version: 0.3.0
#' Adapted from: DADA2 ITS workflow:https://benjjneb.github.io/dada2/ITS_workflow.html

### LIBRARIES 
library(dada2)
library(Biostrings)
library(ShortRead)
library(optparse)


### CMD OPTIONS

# required options trunclenr trunclenl, 
option_list <- list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", action = "store"),
  make_option(c("-f", "--fwd_primer"), type="character", default=NULL, 
              help="forward primer", action = "store"),
  make_option(c("-r", "--rev_primer"), type="character", default=NULL, 
              help="reverse primer", action = "store"),
  make_option(c("-n", "--n_sample"), type="integer", default=4, 
              help="number of paired reads to sample", action = "store")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$path)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (path)", call.=FALSE)
}

## FUNCTIONS
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

search_primers <- function(R1_filepaths, R2_filepaths, FWD.orients, REV.orients,sample_n) {
  # this function randomly checks a defined number of reads for primers
  subsample <- sample(1:sample_n, length(R1_filepaths), replace = T)
  for (i in subsample){
    hits <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = R1_filepaths[[i]]), 
          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = R2_filepaths[[i]]), 
          REV.ForwardReads = sapply(REV.orients, primerHits, fn = R1_filepaths[[i]]), 
          REV.ReverseReads = sapply(REV.orients, primerHits, fn = R2_filepaths[[i]]))
    print(hits)
  }
}

## MAIN
fnFs <- sort(list.files(opt$path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(opt$path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# get all possible orientations of primer 
FWD.orients <- allOrients(opt$fwd_primer)
REV.orients <- allOrients(opt$rev_primer)

# outpath of N filtered reads
fnFs.filtN <- file.path(opt$path, "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(opt$path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = 6)


search_primers(fnFs.filtN, fnRs.filtN, FWD.orients, REV.orients, 
               sample_n = opt$n_sample)
 
system('rm -rf filtN') 
## TO-do

# double check cmd options are set correctly 
