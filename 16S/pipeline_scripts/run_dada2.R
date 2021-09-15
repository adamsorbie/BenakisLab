#!/usr/bin/env Rscript

#' Author: Adam Sorbie 
#' Date: 14/09/21
#' Version: 1.0.2

# to-do Polishing script 

# necessary 

# nice additions 

#' enhance taxonomic classification if possible 
#' add third table with best hit formatting 
#' remove any redundant code 
#' better error handling 

### LIBRARIES 
if(!require("pacman")){
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(dada2, Biostrings, optparse, ggpubr, stringr, tictoc)

### CMD OPTIONS

# required options trunclenr trunclenl, 
option_list <- list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path of read files", action = "store"),
  make_option(c("-f", "--trunc_f"), type="integer", default=NULL, 
              help="length to truncate forward reads to", action = "store"),
  make_option(c("-r", "--trunc_r"), type="integer", default=NULL, 
              help="length to truncate reverse reads to", action = "store"),
  make_option(c("-q", "--int_quality_control"), type="logical", default=TRUE, 
              help="output qc from dada2 run", action = "store"),
  make_option(c("-o", "--out"), type="character", default="dada2_out", 
              help="full path of output folder", action = "store"),
  make_option(c("-n", "--n_errorsF"), type="integer", default =2,
              help="expected errors for forward reads",
              action = "store"), 
  make_option(c("-N", "--n_errorsR"),  type="integer", default =2,
              help=" expected errors for reverse reads",
              action = "store"),
  make_option(c("-t", "--threads"),  type="integer", default=2,
              help="Number of threads",
              action = "store")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check path exists 
if (is.null(opt$path)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (path)", call.=FALSE)
}

# check selected number of threads is ok
if (opt$threads > detectCores()){
  stop("Selected number of threads exceeds available, exiting", call. = FALSE)
}

### FUNCTIONS

# convert to relative abundance
normalize <- function(asv_tab) {
    rel_tab <- t(100 * t(asv_tab) / colSums(asv_tab))
    return(as.data.frame(rel_tab))
}

# strip n chars from right 
RIGHT <- function(x,n){
  substring(x,nchar(x)-n+1)
}

# filter abundance using relative abundance cutoff 
filter_abundance <- function(asv_tab) {
  # keeps ASVs which are above 0.25% rel abund in at least one sample 
  rel <- normalize(asv_tab)
  rel_filt <- rel[rowSums(rel) > 0.25, ]
  idx <- rownames(rel_filt)
  asv_tab_out <- asv_tab[idx, ]
  return(asv_tab_out)
}

checkOS <- function(){
   return(.Platform$OS.type)
}

### MAIN

# check all paths have trailing forward slash 
if (RIGHT(opt$path, 1) != "/") {
  opt$path <- paste0(opt$path, "/")
}
if (RIGHT(opt$out, 1) != "/") {
  opt$out <- paste0(opt$out, "/")
}

# time analysis
tic()
print(paste("ANALYSIS STARTING", Sys.time(), sep=" "))

# set cmd params as variables 
path <- opt$path 
setwd(path)

dir.create(opt$out)

maxE <- c(opt$n_errorsF, opt$n_errorsR)
trunc_params <- c(opt$trunc_f, opt$trunc_r)

# check for illegal characters in fastq filenames and exit if detected 
illegal_chars <- c("-","#", "@", " ")
check_fnames <- str_detect(list.files(pattern = "*.fastq.gz"), 
                           illegal_chars)

if (any(check_fnames == TRUE)){
  print("filenames not in illumina format see: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm, 
        exiting")
  stop()
} else {
  print("filenames ok!")
}

# forward read paths
fnFs <- sort(list.files(pattern="R1_001.fastq.gz", full.names = TRUE))
# reverse read paths
fnRs <- sort(list.files(pattern="R2_001.fastq.gz", full.names = TRUE))

# get sample names from forward reads 
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1); sample.names

# if sample names contain trimmed_primer suffix remove it
if (str_detect(sample.names, "_trimmed_primer") == TRUE){
  sample.names <- gsub(sample.names, pattern = "_trimmed_primer", 
                       replacement = "")
}

# output path for filtered F reads
filtFs <- file.path("filtered", paste0(sample.names, "_filt_R1_001.fastq.gz"))
# output path for filtered R reads
filtRs <- file.path("filtered", paste0(sample.names, "_filt_R2_001.fastq.gz"))

# filter and trim reads 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=maxE, 
                     truncLen=trunc_params, rm.phix=TRUE,
                     truncQ=3, compress=TRUE, multithread=opt$threads)

# learn error rate and plot error profiles 
errF <- learnErrors(filtFs, multithread = opt$threads, nbases = 1e8, 
                    randomize = TRUE)

jpeg(paste(opt$out, "error_plot_f.jpg", sep="/"))
plotErrors(errF, nominalQ=TRUE)
dev.off()

errR <- learnErrors(filtRs, multithread =opt$threads, nbases = 1e8, 
                    randomize = TRUE)
jpeg(paste(opt$out, "error_plot_r.jpg", sep="/"))
plotErrors(errR, nominalQ=TRUE)
dev.off()

# get unique sequences
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# remove filtered reads to reduce disk space 
if (checkOS() == "unix"){
  system("rm -rf filtered || true")
} else if (checkOS() == "windows") {
  system("if exists filtered rmdir /q /s filtered")
}

# run dada2 algorithm 
dadaFs <- dada(derepFs, err=errF, multithread = opt$threads)
dadaRs <- dada(derepRs, err=errR, multithread = opt$threads)

# merge paired reads
merged <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# make sequence table - x - abundances, y - seqs
seqtab <- makeSequenceTable(merged)

# either remove or output 
seq_length_distr <- table(nchar(getSequences(seqtab)))

# filter reads by length
seq_lengths <- as.numeric(names(seq_length_distr))

# remove sequences which are longer or shorter than expected
# 1.25 * mean length is a reasonable rule of thumb 
MINLEN <- round(mean(seq_lengths) / 1.25, 0)  
MAXLEN <- round(mean(seq_lengths) * 1.25, 0)
seqlen <- nchar(getSequences(seqtab))
seqtab.len_filt <- seqtab[, seqlen>MINLEN & seqlen<MAXLEN]

# filter chimeras 
seqtab_chim_filt <- removeBimeraDenovo(seqtab.len_filt, method = "consensus", 
                                       multithread=opt$threads, verbose = TRUE)

percent.nonchim <- sum(seqtab_chim_filt)/sum(seqtab)

# pre-filter by 0.25% relative abundance (Reitmeier et al 2021, ISME communication)
print(paste("pre-filtering dimensions:", dim(t(seqtab_chim_filt))[1], sep=" "))
seqtab_chim_abun_filt <- filter_abundance(t(seqtab_chim_filt))
print(paste("post-filtering dimensions:", dim(seqtab_chim_abun_filt)[1], sep=" ")) 
seqtab_chim_abun_filt <- t(seqtab_chim_abun_filt)

# check number of reads which made it through each step 
getN <- function(x) sum(getUniques(x))
# generate table summarizing steps 
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(merged, getN),
                          nonchim=rowSums(seqtab_chim_filt), 
                          perc_reads_retained=round(rowSums(seqtab_chim_filt)/out[,1]*100, 1), 
                          abund_filt_perc=round(rowSums(seqtab_chim_abun_filt)/rowSums(seqtab_chim_filt)*100, 1))

## QC 
qc_list <- list("Sequence_list_distribution" = seq_length_distr, 
                "Stats" = summary_tab)

# output sequence length distribution and study stats
if (opt$int_quality_control == TRUE){
  
  seq_len_df <- data.frame(seq_length_distr)
  p <- ggpubr::ggbarplot(seq_len_df, x="Var1", y="Freq")
  p <- p + rotate_x_text(angle=45)
  ggsave(paste(opt$out, "sequence_len_dist.png", sep=""), p, device = "png")
  # write out study stats
  write.table(qc_list$Stats, paste(opt$out, "study_stats.txt", sep=""), 
              sep="\t", col.names = NA)
  
}

## Taxonomic classification 

# download silva files and assign taxonomy 

if (!file.exists("silva_nr99_v138.1_train_set.fa.gz")){
  download.file(url="https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1", 
                destfile="silva_nr99_v138.1_train_set.fa.gz")
}

taxa <- assignTaxonomy(seqtab_chim_abun_filt, "silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE, tryRC=TRUE)

if (!file.exists("silva_species_assignment_v138.1.fa.gz")){
  download.file(url="https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1", 
                destfile="silva_species_assignment_v138.1.fa.gz")
}

# addSpecies may improve % classified in some cases
asv_taxa <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz",
                       tryRC=TRUE)

## write output files

# create vector to contain headers
asv_headers <- vector(dim(seqtab_chim_abun_filt) [2], mode = "character")
# generate row and column names 
rownames(asv_taxa) <- gsub(pattern=">", replacement="", x=asv_headers)

# generate fasta headers
for (i in 1:dim(seqtab_chim_abun_filt) [2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")  
}
# fasta file 
asv_seqs <- colnames(seqtab_chim_abun_filt)
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(opt$out, "ASV_seqs.fasta", sep = "/"))

# write feature count table 
asv_tab <- t(seqtab_chim_abun_filt)
row.names(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab, paste(opt$out, "ASV_seqtab.tab", sep = "/"), 
            sep = "\t", quote = F, col.names = NA)


# taxonomy table 
row.names(asv_taxa) <- sub(">", "", asv_headers)
asv_taxa <- as.data.frame(asv_taxa)
asv_taxa$taxonomy <- paste(asv_taxa$Kingdom, asv_taxa$Phylum, asv_taxa$Class, 
                           asv_taxa$Order, asv_taxa$Family, 
                           asv_taxa$Genus, asv_taxa$Species, sep = ";")

# clean up taxa names
asv_taxa$taxonomy <- gsub("NA","", asv_taxa$taxonomy)

# write out feature count table with taxonomy included
asv_tab_tax <- as.data.frame(asv_tab)
asv_tab_tax$taxonomy <- paste(asv_taxa$taxonomy)
write.table(asv_tab_tax, paste(opt$out, "ASV_seqtab_tax.tab", sep = "/"), 
            sep = "\t", quote = F, col.names = NA)

# Phylogenetic tree construction 
setwd(opt$out)

# align seqs with muscle and create tree with FastTree 
if (checkOS() == "unix"){
  system("muscle -in ASV_seqs.fasta -out aligned.fasta -maxiters 3")
  
  system("FastTree -quiet -nosupport -gtr -nt aligned.fasta > ASV_tree.tre")
} else if (checkOS() == "windows"){
  
  #system("muscle.exe -in ASV_seqs.fasta -out aligned.fasta -maxiters 3")
  
  #system("FastTree.exe -quiet -nosupport -gtr -nt aligned.fasta > ASV_tree.tre")
  
  print("phylogenetic tree construction for windows is currently unsupported but
        will be added in the near future")
}
print(paste("ANALYSIS COMPLETED", Sys.time(), sep=" "))

# Send parameters to log file 
sink(paste(opt$out, "parameter_log.txt", sep="/"))
print(paste0("Filepath: ", opt$path))
print(paste0("Forward trunc: ", opt$trunc_f, sep=":"))
print(paste0("Reverse trunc: ", opt$trunc_r, sep=":"))
print(paste0("Internal QC performed: ", opt$int_quality_control))
print(paste0("Output: ", opt$out))
print(paste0("Forward errors: ", opt$n_errorsF))
print(paste0("Reverse errors: ", opt$n_errorsR))
print(paste0("Threads: ", opt$threads))
toc()
closeAllConnections() 
