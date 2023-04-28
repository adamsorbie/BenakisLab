#' Metabolomics data-analysis pipeline 
#' Author: Adam Sorbie
#' Date: 05/05/21
#' Version: 0.2.0 

#' This pipeline uses the POMA R package to generate a quick overview of your
#' metabolomics data. Note it is not intended for in depth analysis, only to 
#' get a general idea of differences. 
 
library(POMA)
library(tidyverse)
library(filesstrings)

## FUNCTIONS 

select_samples <- function(df, meta, sampleid_col){
  
  sample_names <- meta %>% pull(.data[[sampleid_col]])
  df_samples <- df %>% 
    select(all_of(sample_names))
  return(df_samples)
}

POMA_import <- function(metab_path, metadata_path, rowname_col, sampleid_col) {
  # do above, chekcing if matrix in correct format, transpose if not 
  # check rownames and sampleID col match 
  # check for NA in meta, replace with 0 
  
  dat <- read_csv(metab_path) %>% 
    column_to_rownames(var=rowname_col)
  
  metadata <- read_csv(metadata_path)
 
  dat_samples <- select_samples(dat, metadata, sampleid_col) %>% 
    t() %>% 
    as.data.frame()
  
  if (metadata[sampleid_col] == rownames(dat_samples)){
    poma <- PomaMSnSetClass(target = metadata, features = dat_samples)
  }
  else {
    print("Error row names and sample ids do not match")
    stop()
  }
  
  return(poma)
  
}


Metabo_report <- function(metab, metadata, metaboliteid_col, sampleIDcol){
  
  poma_obj <- POMA_import(metab, metadata, metaboliteid_col, sampleIDcol)
  PomaEDA(data = poma_obj, 
          imputation = "half_min", 
          normalization = "log_scaling", 
          clean_outliers = F)
  outpath <- paste(find.package("POMA"), "rmd/POMA_EDA_report.pdf", sep = "/")
  file.move(outpath, getwd())
}

