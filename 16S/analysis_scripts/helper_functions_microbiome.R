library(pheatmap)
library(ComplexHeatmap)
library(tidyverse)

get_top_n <- function(ps, n, level = "ASV") {
  if (level != "ASV") {
    ps <- ps %>% tax_fix %>%
      tax_agg(rank = level)
    ps <- ps$ps
  }
  topn <- ps %>%
    transform(transform = "relative") %>%
    psmelt() %>%
    group_by(OTU) %>%
    summarise(Mean_abund = mean(Abundance)) %>%
    slice_max(Mean_abund, n = n) %>%
    pull(OTU)
  return(topn)
}

# update split on and use your functions where possible
plot_heatmap <- function(ps,
                         list_taxa,
                         variable,
                         heatmap.colors,
                         heatmap.type = "complexheatmap",
                         split_on = NULL,
                         ...) {
  ps_in <- prune_taxa(list_taxa, ps) %>%
    format_to_besthit()
  
  ps_z <- transform(ps_in, "Z")
  
  otu.mat <- abundances(ps_z)
  meta <- meta(ps_z)
  
  subset_meta <- subset(meta, select = c(variable))
  
  # Get genus name
  row_ann <- data.frame(tax_table(ps_z)[, 6])
  # ASV names as rownames
  rownames(row_ann) <- rownames(otu.mat)
  
  if (heatmap.type == "pheatmap") {
    p <- pheatmap::pheatmap(
      otu.mat,
      annotation_col = subset_meta,
      annotation_row = row_ann,
      color = heatmap.colors,
      cluster_cols = F,
      border_color = "white",
      show_colnames = F,
      cellwidth = 12,
      cellheight = 12,
      show_rownames = F,
      ...
    )
    return(p)
  }
  else if (heatmap.type == "complexheatmap") {
    unique_taxa <- unique(row_ann[, 1])
    #paired or dark color brewer palette
    if (length(unique_taxa <= 2)) {
      taxa_palette <- brewer.pal(3, "Dark2")[1:length(unique_taxa)]
    }
    else if (length(unique_taxa > 2)) {
      taxa_palette <- brewer.pal(length(unique_taxa), "Dark2")
    }
    
    names(taxa_palette) <- unique_taxa
    
    
    
    row_ha <-
      rowAnnotation(
        df = row_ann,
        col = list(Genus = taxa_palette),
        show_annotation_name = F
      )
    column_ha <- HeatmapAnnotation(
      df = subset_meta,
      col = list(Tissue = c(
        "NT" =  "#30638E",
        "TA" = "#E67F0D",
        "T" = "#D1495B"
      )),
      border = T,
      show_annotation_name = F
    )
    if (!is.null(split_on)) {
      p <-
        Heatmap(
          otu.mat,
          right_annotation = row_ha,
          top_annotation = column_ha,
          col = heatmap.colors,
          cluster_columns = F,
          show_row_names = F,
          show_column_names = F,
          name = "Z-score abundance",
          row_title = NULL,
          column_split = split_on,
          border = "black",
          column_gap = unit(1.5, "mm")
        )
    } else {
      p <-
        Heatmap(
          otu.mat,
          right_annotation = row_ha,
          top_annotation = column_ha,
          col = heatmap.colors,
          cluster_columns = F,
          show_row_names = F,
          show_column_names = F,
          name = "Z-score abundance",
          row_title = NULL,
          border = "black"
        )
    }
    
    
    return(draw(
      p,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    ))
  }
}


create_dummy_tax <- function(asvtab) {
  if ("ASV" %in% colnames(asvtab)) {
    asvtab <- t(asvtab)
  }
  dummy_tax <-
    matrix(data = colnames(asvtab),
           nrow = ncol(asvtab),
           ncol = 6)
  colnames(dummy_tax) <-
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  rownames(dummy_tax) <- colnames(asvtab)
  
  return(dummy_tax)
}
