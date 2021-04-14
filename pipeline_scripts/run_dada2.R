#!/usr/bin/env Rscript

#' Author: Adam Sorbie 
#' Date: 14/04/21
#' Version: 0.6.0


### LIBRARIES 
library(dada2)
library(DECIPHER)
library(optparse)
library(keypress)

### CMD OPTIONS

# required options trunclenr trunclenl, 
option_list <- list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", metavar="character"),
  make_option(c("-f", "--trunc_f"), type="integer", default=NULL, 
              help="length to truncate forward reads to", metavar="character"),
  make_option(c("-r", "--trunc_r"), type="integer", default=NULL, 
              help="length to truncate reverse reads to", metavar="character"),
  make_option(c("-m", "--minlen"), type="integer", default=NULL, 
              help="minimum length of reads", metavar="character"),
  make_option(c("-q", "--int_quality_control"), type="logical", default=TRUE, 
              help="output qc from dada2 run", metavar="character"),
  make_option(c("-c", "--cut_adapt"), type="logical", default=FALSE, 
              help="logical indicating if primers were removed by cutadapt script"
              , metavar="character"),
  make_option(c("-o", "--out"), type="character", default="dada2_qc", 
              help="full path of output folder", metavar="character"),
  make_option(c("-n", "--n_errorsF"), type="integer", default =2,
              help="atomic vector of expected errors for forward reads",
              metavar = "character"), 
  make_option(c("-N", "--n_errorsR"),  type="integer", default =2,
              help="atomic vector of expected errors for reverse reads",
              metavar = "character"),
  make_option(c("-t", "--threads"),  type="integer", default=detectCores(),
              help="Number of threads",
              metavar = "character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$path)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (path)", call.=FALSE)
}

### FUNCTIONS
normalize <- function(asv_tab) {
    rel_tab <- t(100 * t(asv_tab) / colSums(asv_tab))
    return(as.data.frame(rel_tab))
}

filter_abundance <- function(asv_tab) {
  # keeps ASVs which are above 0.25% rel abund in at least one sample 
  rel <- normalize(asv_tab)
  idx <- rownames(rel)
  asv_tab_out <- asv_tab[idx, ]
  return(asv_tab_out)
}


### MAIN
print(paste("ANALYSIS STARTING", Sys.time(), sep=" "))

path <- opt$path 
setwd(path)

dir.create(opt$out)

maxE <- c(opt$n_errorsF, opt$n_errorsR)
trunc_params <- c(opt$trunc_f, opt$trunc_r)

# invoke system command ls and cut to get sample names from fastq filenames
system('ls *R1_001.fastq.gz | cut -f1 -d "_" > samples')
# store sample names as variable
samples <- scan("samples", what = "character")
# forward reads
fnFs <- sort(list.files(pattern="R1_001.fastq.gz", full.names = TRUE))
# reverse reads
fnRs <- sort(list.files(pattern="_R2_001.fastq.gz", full.names = TRUE))


# output path for filtered F reads
filtFs <- file.path("filtered", paste0(samples, "_R1_filt.fastq.gz"))
# output path for filtered R reads
filtRs <- file.path("filtered", paste0(samples, "_R2_filt.fastq.gz"))

# remove filtered reads 
system("rm -rf filtered || true")

if (opt$cut_adapt != TRUE) {
    trimFs <- file.path("trimmed_primer", paste0(samples, "_R1_trimmed_primer.fastq.gz"))
    trimRs <- file.path("trimmed_primer", paste0(samples, "_R2_trimmed_primer.fastq.gz"))
    
    out <- filterAndTrim(trimFs, filtFs, trimRs, filtRs, maxEE=maxE, 
                         truncLen=trunc_params, rm.phix=TRUE, minLen = opt$minlen,
                         truncQ=3, compress=TRUE, multithread=opt$threads)
} else {
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=maxE, 
                       truncLen=trunc_params, rm.phix=TRUE, minLen = opt$minlen,
                       truncQ=3, compress=TRUE, multithread=opt$threads)
}
# learn error rate 
errF <- learnErrors(filtFs, multithread = opt$threads, nbases = 1e8)
jpeg(paste(opt$out, "error_plot_f.jpg", sep="/"))
plotErrors(errF, nominalQ=TRUE)
dev.off()
errR <- learnErrors(filtRs, multithread =opt$threads, nbases = 1e8)
jpeg(paste(opt$out, "error_plot_r.jpg", sep="/"))
plotErrors(errR, nominalQ=TRUE)
dev.off()

# get unique sequences
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- samples
names(derepRs) <- samples

# run dada2 algorithm 
dadaFs <- dada(derepFs, err=errF, multithread = opt$threads)
dadaRs <- dada(derepRs, err=errR, multithread = opt$threads)

# merge paired reads
merged <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# make sequence table - x - abundances, y - seqs
seqtab <- makeSequenceTable(merged)

# either remove or output 
seq_length_distr <- table(nchar(getSequences(seqtab)))

# filter reads by length, removing those below 300/ above 600
seq_lengths <- as.numeric(names(seq_length_distr))

# this is reasonable rule of thumb but  might be good to use quantiles in future
MINLEN <- round(mean(seq_lengths) / 1.25, 0)  
MAXLEN <- round(mean(seq_lengths) * 1.25, 0)
seqlen <- nchar(getSequences(seqtab))
seqtab.len_filt <- seqtab[, seqlen>MINLEN & seqlen<MAXLEN]

# filter chimeras 
seqtab_chim_filt <- removeBimeraDenovo(seqtab.len_filt, method = "consensus", 
                                       multithread=opt$threads, verbose = TRUE)

percent.chim <- 1 - sum(seqtab_chim_filt)/sum(seqtab)
# check number of reads which made it through each step 
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(merged, getN),
                          nonchim=rowSums(seqtab_chim_filt), 
                          nonchim_perc=round(rowSums(seqtab_chim_filt) / sapply(merged, getN)*100, 1),
                          final_perc_reads_retained=round(rowSums(seqtab_chim_filt)/out[,1]*100, 1))

## QC 
qc_list <- list("Sequence_list_distribution" = seq_length_distr, 
                "Stats" = summary_tab)

if (opt$int_quality_control == TRUE){
  
  write.table(qc_list$Stats, paste(opt$out, sep="/"), 
              "study_stats.txt", sep="\t")
  
}

## downloading DECIPHER-formatted SILVA v138 reference
if (file.exists("SILVA_SSU_r138_2019.RData") == FALSE){
  download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", 
                destfile="SILVA_SSU_r138_2019.RData")
}

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")

## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab_chim_filt))

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", 
                   processors=NULL, verbose=TRUE)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_taxa <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa
}))

# generate row and column names
asv_headers <- vector(dim(seqtab_chim_filt) [2], mode = "character")
colnames(asv_taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(asv_taxa) <- gsub(pattern=">", replacement="", x=asv_headers)


## write out files
asv_seqs <- colnames(seqtab_chim_filt)


for (i in 1:dim(seqtab_chim_filt) [2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")  
}


# fasta file 

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(opt$out, "ASV_seqs.fasta", sep = "/"))

# count table 

asv_tab <- t(seqtab_chim_filt)

row.names(asv_tab) <- sub(">", "", asv_headers)

# pre-filter by 0.25% relative abundance (Reitmeier et al 2020, Researchsquare)
asv_tab_filt <- filter_abundance(asv_tab)

write.table(asv_tab_filt, paste(opt$out, "ASV_seqtab.tab", sep = "/"), 
            sep = "\t", quote = F, col.names = NA)

# tax table 

row.names(asv_taxa) <- sub(">", "", asv_headers)
asv_taxa <- as.data.frame(asv_taxa)
asv_taxa$taxonomy <- paste(asv_taxa$Kingdom, asv_taxa$Phylum, asv_taxa$Class, 
                           asv_taxa$Order, asv_taxa$Family, 
                          asv_taxa$Genus, asv_taxa$Species, sep = ";")

asv_tab_tax <- as.data.frame(asv_tab)
asv_tab_tax$taxonomy <- paste(asv_taxa$taxonomy)
write.table(asv_tab_tax, paste(opt$out, "ASV_seqtab_tax.tab", sep = "/"), 
            sep = "\t", quote = F, col.names = NA)


# Phylogenetic tree construction 
setwd(opt$out)

system("muscle -in ASV_seqs.fasta -out aligned.fasta -maxiters 3")
system("FastTree -quiet -nosupport -gtr -nt aligned.fasta > ASV_tree.tre")

print(paste("ANALYSIS COMPLETED", Sys.time(), sep=" "))

# to-do Polishing script 

# necessary 

#' re-check code based on astrobiomike/DADA2 tutorial
#' refine settings to reduce number of ASVs if possible 

# nice additions 

#' enhance taxonomic classification if possible 
#' add third table with best hit formatting 
#' remove any redundant code 
#' tidy up writing out section
#' potentially update muscle to mafft 
