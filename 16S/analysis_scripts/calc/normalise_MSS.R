library(phyloseq)

transform_mss <- function(x) {
  x_untransformed <- x
  if (length(is(x)) == 1 && class(x) == "phyloseq"){
    x <- abundances(x)
  }
  else {
    print("not a phyloseq or data not integers, exiting")
    stop()
  }
  xt <- t(min(colSums(x)) * t(x) / colSums(x))
  otu_table(x_untransformed)@.Data <- xt
  x_out <- x_untransformed
  return(x_out)
}
    

   