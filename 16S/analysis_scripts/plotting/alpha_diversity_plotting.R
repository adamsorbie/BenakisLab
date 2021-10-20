select_var <- function(df, column, alpha_cols){
  cols <- c(alpha_cols, column) # should make that modular 
  df_out <- df[cols]
  #cols <- colnames(df_out)
  #df_out_m <- gather(df_out, cols)
  df_out_m <- melt(df_out)
  return(df_out_m)
}

#' plotting function 
plot_alpha <- function(df, variable_col="variable", value_col="value", 
                       fill_var="fill", comparisons_list, xlab, ylab, 
                       p_title, multiple_groups=TRUE, col_palette){
  
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
  final_plot <- plot + stat_compare_means(comparisons = comparisons_list, label = "p.signif", test = stat_method ) + 
    theme_bw() + ggtitle(p_title) + theme(axis.text.x = element_text(size = 10), 
                                          axis.line = element_line(colour = "black"), 
                                          panel.grid.major = element_blank(), 
                                          panel.grid.minor = element_blank(), 
                                          panel.background = element_blank(),
                                          panel.border = element_blank(), 
                                          plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_palette )
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


reorder_factors <- function(combined_df, select_col, factor_levels, alpha_metric) {
  combined_df <- select_var(combined_df, select_col, alpha_metric)
  combined_df[[select_col]] <- factor(combined_df[[select_col]], levels = factor_levels)
  combined_df_alpha <- combined_df[combined_df$variable == alpha_metric, ]
  return(combined_df_alpha)
}

ggsave_default_options <- function(filename, plot) {
  ggsave(filename, plot = plot, units = "in", width = 5, height = 5, dpi = 500 )
  
}