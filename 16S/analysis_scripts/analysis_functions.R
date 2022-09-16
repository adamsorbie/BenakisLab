if(!require("pacman")){
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load(tidyverse, ggpubr, rstatix, picante, phyloseq,
               vegan, ANCOMBC, phangorn, GUniFrac, zoo)

######################### DATA WRANGLING #########################  

# filter dataframe by rownames
filter_rownames <- function(df, filt) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}

meta_to_df <- function(ps) {
  return(as(sample_data(ps), "data.frame"))
} 

ps_to_asvtab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}


#########################    IMPORTING DATA #########################  
split_taxonomy <- function(otu) {
  # select taxa column and split taxonomy column into separate columns using ; delimiter
  taxonomy_phylo <- dplyr::select(otu, "taxonomy") %>%
    separate("taxonomy", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";")
  return(taxonomy_phylo)
}

format_taxonomy <- function(ps) {
  ps_tax <- taxonomy(ps) %>% 
    mutate_all(na_if,"")
  
  ps_tax[] <- t(na.locf(t(ps_tax))) %>% 
    as.data.frame()
  # add unknown to genera which match family
  ps_tax <- ps_tax %>%
    mutate(Genus = case_when(Genus == Family ~ paste0("unknown_", Genus),
                             TRUE ~ Genus)) %>% 
    mutate(Family = case_when(Family == Order ~ paste0("unknown_", Family),
                           TRUE ~ Family)) %>% 
    mutate(Order = case_when(Order == Class ~ paste0("unknown_", Order),
                             TRUE ~ Order)) %>% 
    mutate(Class = case_when(Class == Phylum ~ paste0("unknown_", Class),
                             TRUE ~ Class)) %>% 
    mutate(Phylum = case_when(Phylum == Kingdom ~ paste0("unknown_", Phylum),
                             TRUE ~ Phylum)) %>% 
    mutate(ASV = rownames(ps_tax)) %>% 
    mutate(Highest_classified = paste(rownames(ps_tax), Genus, sep="; "))
  
  rownames(ps_tax) <- ps_tax[, "Highest_classified"]
  taxa_names(ps) <- ps_tax[, "Highest_classified"]
  tax_table(ps) <- tax_table(as.matrix(ps_tax))

  return(ps)
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

import_as_pseq <- function(otu, mapping, tree=NULL){
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

######################### NORMALISATION #########################

# normalise data with minimum sum scaling 
transform <- function(ps, transform="mss") {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq"){
    x <- ps_to_asvtab(ps)
  }
  else {
    print("not a phyloseq object, exiting")
    stop()
  }
  
  if (transform %in% c("mss", "relative")){
    if (transform == "mss") {
      ps_t <- t(min(colSums(x)) * t(x) / colSums(x))
      
    } else if (transform == "relative"){
      ps_t <- t(100 * t(x) / colSums(x))
      
    } 
    otu_table(ps)@.Data <- ps_t
    
    return(ps)
    
  } else {
    print("Not a valid transform, exiting")
    stop()
  }
 
}

######################### PLOTTING AND STATS #########################

# calculate adonis R2 and p-value
phyloseq_adonis <- function(ps, dist_matrix, factor, ...) {
  
  
  meta_df <- meta_to_df(ps) 
  
  if(sum(is.na(meta_df[[factor]])) > 0){
    print("metadata contains NAs, remove these samples with subset_samples
          before continuing")
    return (NULL)
  } else {
    ps_ad <- adonis(unname(dist_matrix) ~ meta_df[[factor]],
                    data = meta_df, ...)
    return(ps_ad)
  }
}


# create statistical formula 
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}

# plot boxplot with stats
plot_boxplot <- function(df, variable_col, value_col, 
                         fill_var="fill", comparisons_list, xlab, ylab, 
                         p_title=NULL, multiple_groups=FALSE, col_palette,
                         group.order=NULL, paired=FALSE, ...){
  # extend color palette with transparent value - required due to way we are 
  # layering plot 
  
  col_palette <- c(col_palette, "transparent")
  
  if (!is.null(group.order)){
    df[, variable_col] <- factor(df[, variable_col], levels = group.order)
  }
  
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == TRUE) {
    if (paired == TRUE){
      stat_variance <- df %>% 
        friedman_test(formula)
      stat_test <- df %>%  
        pairwise_wilcox_test(formula, comparisons = comparisons_list,
                             p.adjust.method = "BH", paired=TRUE) %>% 
        add_significance() %>% 
        add_xy_position(x = variable_col) %>% 
        filter(p.adj < 0.05)
    }
    else {
      stat_variance <- df %>% 
        kruskal_test(formula)
      stat_test <- df %>%  
        pairwise_wilcox_test(formula, comparisons = comparisons_list,
                             p.adjust.method = "BH") %>% 
        add_significance() %>% 
        add_xy_position(x = variable_col) %>% 
        filter(p.adj < 0.05)
    }
  } 
  else if (multiple_groups == FALSE) {
    if (paired == TRUE){
      stat_test <- df %>% 
        wilcox_test(formula, paired=TRUE) %>% 
        add_significance() %>% 
        add_xy_position(x = variable_col) %>% 
        filter(p < 0.05)
    }
    else {
      stat_test <- df %>% 
        wilcox_test(formula) %>% 
        add_significance() %>% 
        add_xy_position(x = variable_col) %>% 
        filter(p < 0.05) 
    }

  }
  
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(df, aes_string(x=variable_col, y=value_col, 
                                fill = variable_col, color=variable_col )) + 
    geom_boxplot(color="black", alpha=0.8,outlier.shape=5, outlier.size=1) + 
    geom_point(size = 1.5, position = position_jitterdodge()) + 
    labs(x=xlab, y=ylab) + 
    stat_boxplot(color="black", geom = "errorbar", width=0.2) 
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot + 
    theme_classic2() + 
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
          plot.title = element_text(hjust = 0.5), 
          legend.position = "None") +
    scale_fill_manual(values = col_palette) + 
    scale_color_manual(values = col_palette) + 
    rotate_x_text(angle = 45)
  
  if (dim(stat_test)[1] == 0){
    plot_out <- final_plot
  }
  else {
    if (multiple_groups == T){
      plot_out <- final_plot +
        stat_pvalue_manual(stat_test, label = "p.adj.signif", hide.ns = T,
                           inherit.aes = FALSE, ...)
    }
    else {
      plot_out <- final_plot +
        stat_pvalue_manual(stat_test, label = "p.signif", hide.ns = T,
                           inherit.aes = FALSE, ...)
    }
  }
  
  return(plot_out)
}

# plot scatter plot with correlation if desired
plot_scatter <- function(df, x, y, point_color, 
                         line_color, fill_color, 
                         xlabel, ylabel, 
                         corr.method=NULL, ...) {
  p <- ggplot(data=df, mapping = aes(x = .data[[ x ]], y = .data[[ y ]] )) +
    geom_point(aes(color=point_color), size=2.5) +
    geom_smooth(method = "lm", color=line_color, fill=fill_color) +
    theme_bw() + 
    theme(legend.position = "None", 
          axis.title.x = element_text(size=14),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12)) + 
    xlab(xlabel) + 
    ylab(ylabel)
  
  if (!is.null(corr.method)) {
    p <- p + stat_cor(method=corr.method, ...)
    return(p)
  }
  else {
    return(p)
  }
}

######################### ALPHA DIVERSITY #########################

# calculate shannon effective
Shannon.E <- function(x) {
  summed <- sum(x)
  shannon.e <- round(exp(-sum(x[x>0]/summed*log(x[x>0]/summed))), digits = 2)
  return(shannon.e)
}

# calculate faiths phylogenetic diversity
Faiths <- function(x, tree) {
  pd <- pd(x, tree, include.root = F)
  return(pd[,1])
}

# calculate richness
Richness <- function(x, detection=0.5) {
  observed <- sum(x>detection)
  return(observed)
}

# calculate all alpha diversity matrices and return dataframe 
calc_alpha <- function(ps, ...) {
  mat_in <- ps_to_asvtab(ps) %>% 
    t()
  
  tree <- phy_tree(ps)
  
  diversity <- setNames(data.frame(matrix(ncol = 3, nrow = nsamples(ps))), 
                        c("Richness", "Shannon.Effective", "Faiths.PD"))
  rownames(diversity) <- rownames(meta_to_df(ps))
  
  diversity$Richness <- apply(mat_in, 1, Richness, ...)
  diversity$Shannon.Effective <- apply(mat_in, 1, Shannon.E)
  diversity$Faiths.PD <- Faiths(mat_in, tree)
  
  return(diversity)
}



######################### BETA DIVERSITY #########################
### wrap these into two functions 

# Calculate Gunifrac in phyloseq objects 
phyloseq_gunifrac <- function(ps, asdist=TRUE) {
  # input phyloseq obj
  # returns dist matrix of GUnifrac distance which can 
  # be used as input for ordinate() function 
  
  #extract otu table and tree
  otu_tab <- t(ps_to_asvtab(ps))
  tree <- phy_tree(ps)
  
  # root tree at midpoint
  rooted_tree <- midpoint(tree)
  
  # calculate gunifrac
  gunifrac <- GUniFrac(otu_tab, rooted_tree, 
                       alpha = c(0.0,0.5,1.0))$unifracs
  # alpha param is weight on abundant lineages, 0.5 has best power
  # so we will extract this 
  if (asdist == T){
    return(stats::as.dist(gunifrac[, , "d_0.5"]))
  }
  else{
    return(gunifrac[, , "d_0.5"])
  }
}


calc_betadiv <- function(ps, dist, ord_method="NMDS") {
  if (ord_method %in% c("NMDS", "MDS", "PCoA")) {
    
    
    if (dist %in% c("unifrac", "wunifrac", "bray")) {
      
      dist_metric <- distance(ps, dist)
      ord <- ordinate(ps, ord_method, dist_metric)
      
      return_list <- list("Distance_Matrix" = dist_metric, 
                          "Ordination" = ord$points)
      return(return_list)
    }
    else if (dist %in% c("gunifrac")) {
      
      gu <- phyloseq_gunifrac(ps)
      ord <- ordinate(ps, ord_method, gu)
      
      return_list <- list("Distance_Matrix" = gu, 
                          "Ordination" = ord$points)
      return(return_list)
    }
  else {
    print("Ordination method not supported, supported methods are:
          NMDS, MDS, PCoA")
    stop()
  }
} 
}

plot_beta_div <- function(ps, ordination, dist_matrix, group_variable, add_ellipse=FALSE, cols) {
  # add support for colors and ellipses
  # significance
  ad <- phyloseq_adonis(ps, dist_matrix, group_variable)
  
  # check sample data dimensions >1 otherwise phyloseq 
  # fails to colour samples due to bug:https://github.com/joey711/phyloseq/issues/541 
  if (dim(sample_data(ps))[2] <2) {
    # add repeat of first column as dummy 
    sample_data(ps)[ , 2] <- sample_data(ps)[ ,1]
  }
  plot <- plot_ordination(ps, ordination , color=group_variable)
  plot$layers[[1]] <- NULL
  
  plot_out <- plot + geom_point(size=3, alpha=0.75) +
    theme_bw() + 
    scale_fill_manual(values=cols) +
    scale_color_manual(values = cols) +
    labs(caption = bquote(Adonis~R^2 ~ .(round(ad$aov.tab$R2[1], 2)) ~~ p-value ~ .(ad$aov.tab$`Pr(>F)`[1])))
  
  return(plot_out)
}



######################### Differential Abundance #########################
add_taxonomy_da <- function(ps, list_da_asvs, da_res) {
  
  taxonomy <- as.data.frame(tax_table(ps))
  
  da_taxonomy <- taxonomy %>%
    select(Highest_classified)%>% 
   filter(Highest_classified %in% list_da_asvs)

  da_tax_out <- da_res %>%
    rownames_to_column(var="Group") %>%
    column_to_rownames("ASV") %>%
    merge(da_taxonomy, by=0)
  
  return(da_tax_out)
}

ancom_da <- function(ps, formula, group, ord=NULL, zero_thresh=0.9) {
  
  if (!is.null(ord)){
    group_levels <- levels(get_variable(ps, group))
    sample_data(ps)[[group]] <- factor(group_levels, levels = ord)
  }
 
  
  ps_f <- format_taxonomy(ps)
  
  res <- ancombc(phyloseq = ps_f, formula = formula, p_adj_method = "BH", 
                 zero_cut = zero_thresh, group = group, struc_zero = TRUE, 
                 neg_lb = FALSE, tol = 1e-5, max_iter = 100, conserve = TRUE, 
                 alpha = 0.05, global = FALSE)
  
  res_df <- data.frame(
    ASV = row.names(res$res$beta),
    beta = unlist(res$res$beta),
    se = unlist(res$res$se),
    W = unlist(res$res$W),
    pval = unlist(res$res$p_val),
    qval = unlist(res$res$q_val),
    da = unlist(res$res$diff_abn))

  res_da <- res_df %>%
    filter(da == T)

  da_asvs <- res_da$ASV
  da_tax <- add_taxonomy_da(ps_f, da_asvs, res_da)
  
  return(da_tax)
}



plot_da <- function(ancom_da, groups, cols) {
  
  ancom_da_plot <- ancom_da %>% 
    mutate(enriched_in = ifelse(beta > 0, groups[2], 
                                groups[1]))
  
  ancom_da_plot_sort <- ancom_da_plot %>% 
    arrange(beta) %>% 
    mutate(Highest_classified=factor(Highest_classified, levels=Highest_classified))
  
  
  p <- ggplot(ancom_da_plot_sort, aes(x = Highest_classified, y = beta, color = enriched_in)) +
    geom_point(size = 5) +
    labs(y = paste("\nLog2 Fold-Change", groups[1], "vs", groups[2], sep=" "), x = "") +
    theme_bw() + 
    theme(axis.text.x = element_text(color = "black", size = 14),
          axis.text.y = element_text(color = "black", size = 14),
          axis.title.y = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = "none") +
    coord_flip() +
    geom_hline(yintercept = 0, linetype="dotted") + 
    scale_color_manual(values=cols)
  
  print(p)
}








