transform_mss <- function(x) {
  if (length(is(x)) == 1 && class(x) == "phyloseq"){
    x <- abundances(x)
  }
  else {
    print("not a phyloseq or data not integers, exiting")
    stop()
  }
  xt <- t(min(colSums(x)) * t(x) / colSums(x))
  return(xt)
}
    

   