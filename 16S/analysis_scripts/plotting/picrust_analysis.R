library(tidyverse)
library(factoextra)
library(ggfortify)
library(cowplot)

normalize <- function(df, method) {
  if (method == "relative") {
    rel_tab <- t(100 * t(df) / colSums(df))
    return(as.data.frame(rel_tab))
  }
  else {
    min_sum <- min(colSums(df))
    norm_tab <- t(min_sum * t(df) / colSums(df))
    return(as.data.frame(norm_tab))
  }
}

filter_prev_abund <- function(meta, abund, prev) {
  meta_filt <- meta[which(rowSums(meta > abund) 
                          >= prev * ncol(meta)), ]
} 


picrust_pca <- function(picrust_matrix, column_1, column_2, meta_variable) {
  # picrust_matrix - picrust data [matrix]
  # column 1 - 1st numeric column [string]
  # last numeric column [string]
  # metadata variable to colour by [string] 
  #' @description This function generates a PCA plot from picrust data
  #' @param picrust_matrix: matrix matrix of picrust data (samples are rows), 
  #' also containing metadata information   
  #' @param column_1 string 1st numeric column
  #' @param column_2 string last numric column
  #' @param meta_variable string metadata variable to color PCA by
  #' @usage picrust_pca(picrust_matrix, column_1, column_2, meta_variable)
  #' @return PCA plot colored by metavariables
  pca <- prcomp(picrust_matrix[c(col1:col2)], scale. = T)
  scree_plot <- fviz_eig(pca)
  pca_plot <- autoplot(pca, data = ko_tumor_filt_meta,
                       colour = meta_variable) + 
    geom_point(aes(color= {{ meta_variable }} ), size=3, alpha=0.75) + 
    scale_color_manual(values = c("#05A8AA", "#E26D5A")) + 
    theme_cowplot()
  pca_plot$layers[[1]] <- NULL
  return(p)
}