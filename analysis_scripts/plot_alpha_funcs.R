library(ggplot2)
library(ggpubr)
library(tidyverse)


# this needs to be agnostic to alpha metric 
select_var <- function(df, column, pivot=F){ 
  # list of numeric columns
  alpha_cols <- c("Richness", "Normalized.Richness", "Effective.Richness", "Shannon.Index", 
                  "Shannon.Effective", "Simpson.Index", "Simpson.Effective", "Evenness")

  cols <- c(alpha_cols, column)
  if (pivot == T){
  df_out_long <- select(df, all_of(cols)) %>% 
    pivot_longer(alpha_cols)
  return(df_out_long)
  }
  else {
    df_out <- select(df, all_of(cols))
    return(df_out)
  }
}

select_metric <- function(df, column, alpha_metric){
  combined_df <- select_var(df, column, pivot = T)
  df_metric <- filter(combined_df, combined_df["name"] == alpha_metric)
  return(df_metric)
}

#' plotting function 
plot_alpha <- function(df, variable_col="variable", value_col="value", 
                       fill_var="fill", comparisons_list, xlab, ylab, 
                       p_title, multiple_groups=TRUE){
  
  if (multiple_groups == TRUE) {
    stat_method <- "kruskal.test" 
  } 
  else if (multiple_groups == FALSE) {
    stat_method <- "wilcox.test"
  } 
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(df, aes_string(x=variable_col, y=value_col, 
                                fill = variable_col )) + geom_boxplot(color="black") + 
    labs(x=xlab, y=ylab) + 
    stat_boxplot(geom = "errorbar", width=0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot + stat_compare_means(comparisons = comparisons_list, label = "p.signif", 
                                          test = stat_method, hide.ns = TRUE) + 
    theme_bw() + ggtitle(p_title) + theme(axis.text.x = element_text(size = 10), 
                                          axis.line = element_line(colour = "black"), 
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), 
                                          panel.background = element_blank(),
                                          panel.border = element_blank(), 
                                          plot.title = element_text(hjust = 0.5))
  return(final_plot)
}



make_combined <- function(mapping, alpha_div) {
  meta <- read.table(alpha_div, sep="\t", header = T, 
                     row.names = 1, check.names = F, comment.char = "")
  alpha <- read.table(mapping, sep="\t", header = T, 
                      row.names = 1, check.names = F, comment.char = "")
  combined <- merge(alpha, meta, by = 0)
  row.names(combined) <- combined$Row.names
  combined <- subset(combined, select = -c(Row.names))
  return(combined)
}

filter_combined <- function(combined_df, col, filter_by, not_equal=F) {
  if (not_equal==T){
    combined_out <- combined_df[combined_df[[col]] != filter_by, ]
    return(combined_out)
  }
  else {
    combined_out <- combined_df[combined_df[[col]] == filter_by, ]
    return(combined_out)
  }
}

reorder_factors <- function(combined_df, select_col, reorder_factors=T, factor_levels, alpha_metric) {
  combined_sel <- select_var(combined_df, select_col)
  if (reorder_factors == T){
    combined_sel[[select_col]] <- factor(combined_sel[[select_col]], levels = factor_levels)
  }
  combined_sel_alpha <- select_metric(combined_sel, select_col, alpha_metric)
  return(combined_sel_alpha)
}

ggsave_default_options <- function(filename, plot) {
  ggsave(filename, plot = plot, units = "in", width = 5, height = 5, dpi = 500 )
  
}


