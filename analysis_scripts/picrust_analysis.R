
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