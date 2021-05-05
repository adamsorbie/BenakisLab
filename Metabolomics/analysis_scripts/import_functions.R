#' Functions to handle imputation, and scaling of metabolome data, while handling
#' feature metadata

metabo_list <- function(df, meta, sampleid_col){
  
  sample_names <- meta %>% pull(.data[[sampleid_col]])
  
  features <- df %>% 
    select(all_of(sample_names))
  
  feat_meta <- df %>% 
    select(!all_of(sample_names))
  
  return_list <- list("Features" = features, "Feature_Metadata" = feat_meta)
  return(return_list)
}


half_min_imp <- function(metabo_list) {
  
  # colwise min, fill NA with 0.5 * min 
  half_min <- function(x) replace(x, is.na(x), min(x, na.rm = T) * 0.5)
  df_imp <- replace(metabo_list$Features, TRUE, lapply(metabo_list$Features, half_min))
  
  metabo_list$Features <- df_imp
  
  return(metabo_list)
}

filter_prevalence <- function(metabo_list, prev_trh=0.5) {
  
  feat_filt <- metabo_list$Features[which(rowMeans(!is.na(metabo_list$Features)) > prev_trh), ]
  
  feat_meta_filt <- metabo_list$Feature_Metadata[rownames(feat_filt), ]
  
  return_list <- list("Features" = feat_filt, 
                      "Feature_Metadata" = feat_meta_filt)
  return(return_list)
}