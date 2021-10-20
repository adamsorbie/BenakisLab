library(ggplot2)
library(ggpubr)
library(tidyverse)
library(sigminer)
library(rstatix)
library(cowplot)

# this needs to be agnostic to alpha metric 
select_var <- function(df, column, pivot=F, pd=F){ 
  # list of numeric columns
  if (pd == F){
    alpha_cols <- c("Richness", "Normalized.Richness", "Effective.Richness", "Shannon.Index", 
                  "Shannon.Effective", "Simpson.Index", "Simpson.Effective", "Evenness")
  } 
  else if (pd == T){
    alpha_cols <- c("PD", "SR")
  }
   
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

xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

select_metric <- function(df, column, alpha_metric, pd=F){
  combined_df <- select_var(df, column, pivot = T, pd=pd)
  df_metric <- filter(combined_df, combined_df["name"] == alpha_metric)
  return(df_metric)
}

#' plotting function 
plot_alpha <- function(df, variable_col, value_col, 
                       fill_var="fill", comparisons_list, xlab, ylab, 
                       p_title, multiple_groups=TRUE, col_palette, ...){
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == TRUE) {
    stat_variance <- df %>% 
      kruskal_test(formula)
    stat_test <- df %>%  
      pairwise_wilcox_test(formula, comparisons = comparisons_list,
                           p.adjust.method = "BH") %>% 
      add_significance() %>% 
      add_xy_position(x = variable_col) %>% 
      filter(p.adj < 0.05)
  } 
  else if (multiple_groups == FALSE) {
    stat_test <- df %>% 
      wilcox_test(formula) %>% 
      add_significance() %>% 
      add_xy_position(x = variable_col) %>% 
      filter(p < 0.05)
  }

  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(df, aes_string(x=variable_col, y=value_col, 
                                fill = variable_col )) + geom_boxplot(color="black") + 
    labs(x=xlab, y=ylab) + 
    stat_boxplot(geom = "errorbar", width=0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot + 
    theme_cowplot() + 
    ggtitle(p_title) + 
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          panel.border = element_blank(), 
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = col_palette )

  if (dim(stat_test)[1] == 0){
  plot_out <- final_plot
  }
  else {
  plot_out <- final_plot +
  stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T,
                     step.increase = 0.05, inherit.aes = FALSE, ...)
  }
  
  return(plot_out)
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

reorder_factors <- function(combined_df, select_col, reorder_factors=T, factor_levels, alpha_metric, pd=F) {
  combined_sel <- select_var(combined_df, select_col, pd=pd)
  if (reorder_factors == T){
    combined_sel[[select_col]] <- factor(combined_sel[[select_col]], levels = factor_levels)
  }
  combined_sel_alpha <- select_metric(combined_sel, select_col, alpha_metric, pd=pd)
  return(combined_sel_alpha)
}

ggsave_default_options <- function(filename, plot) {
  ggsave(filename, plot = plot, units = "in", width = 5, height = 5, dpi = 500 )
}


