#' Author Adam Sorbie
#' Date: 2019-08-01
library(tidyverse)
library(phyloseq)

split_taxonomy <- function(otu) {
  # select taxa column and split taxonomy column into separate columns using ; delimiter
   taxonomy_phylo <- dplyr::select(otu, "taxonomy") %>%
    separate("taxonomy", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")
  return(taxonomy_phylo)
}

read_tab_delim <- function(df) {
  # read all tab delimited files using these params
  df_out <- read.table(df, row.names = 1, header = 1, sep="\t", comment.char = "", check.names = F )
  return(df_out)
}

load_phylo <- function(otu, taxa, mapping, tree=NULL) {
  # convert to phyloseq and return list 
  phylo_otu <- otu_table(otu, taxa_are_rows = T)
  
  phylo_tax <- tax_table(as.matrix(taxa))
  
  phylo_map <- sample_data(mapping)
  
  if (exists("tree")){
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_otu, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_otu, phylo_tax, phylo_map))
  }
}

prep_4_phyloseq <- function(otu, mapping, tree=NULL){
  # read files
  otu_taxa <- read_tab_delim(otu)
  metadata <- read_tab_delim(mapping)
  # convert otu into matrix and drop taxonomy col
  count_matrix <- as.matrix(subset(otu_taxa, select=-taxonomy))
  # split taxonomy
  taxonomy_phylo <- split_taxonomy(otu_taxa)
  # make sure taxa table row names match otu 
  rownames(taxonomy_phylo) <- rownames(otu_taxa)
  # use load phylo function to convert all files to phyloseq objs
  if (exists("tree")){
    # if user wants to input a phylogenetic tree this code is called
    out <- load_phylo(count_matrix, taxonomy_phylo, metadata, tree = tree)
  }
  else {
  # without a tree
    out <- load_phylo(count_matrix, taxonomy_phylo, metadata)
  }
  return(out)
}
