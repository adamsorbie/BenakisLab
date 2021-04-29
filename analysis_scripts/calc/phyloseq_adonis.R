library(tidyverse)

phyloseq_adonis <- function(phylo_obj, dist_matrix, factor, ...) {
  metadata_df <- as(sample_data(phylo_obj), "data.frame")
  
  tumor_tissue_ad <- adonis(dist_matrix ~ metadata_df[[factor]],
                            data = metadata_df, ...)
  return(tumor_tissue_ad)
}
