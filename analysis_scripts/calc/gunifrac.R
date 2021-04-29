# Calculate Gunifrac in phyloseq objects 

library(phangorn)
library(GUniFrac)

phyloseq_gunifrac <- function(phyloseq_obj, asdist=TRUE) {
  # input phyloseq obj
  # returns dist matrix of GUnifrac distance which can 
  # be used as input for ordinate() function 
  
  #extract otu table and tree
  otu_tab <- t(otu_table(phyloseq_obj))
  tree <- phy_tree(phyloseq_obj)
  
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

otu_gunifrac <- function(otu, tree, asdist=TRUE){
  rooted_tree <- midpoint(tree)
  
  gunifrac <- GUniFrac(otu, rooted_tree, 
                       alpha = c(0.0, 0.5, 1.0))$unifracs
  if (asdist == T){
    return(stats::as.dist(gunifrac[, , "d_0.5"]))
  }
  else{
    return(gunifrac[, , "d_0.5"])
  }
}