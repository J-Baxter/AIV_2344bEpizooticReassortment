ml_treefiles <- list.files(path = './2024Jun01/region_alignments_withBLAST/iqtree',
                           pattern = 'treefile',
                           full.names = T)

ml_trees <- lapply(ml_treefiles, read.tree) %>%
  lapply(., function(x) {
    x$tip.label <-  gsub("(.*)(\\d{4}-\\d{2}(?:-\\d{2})?).*", "\\1\\2", x$tip.label) %>%
      sub("\\.(?=[^\\.]*$)", "|", ., perl = TRUE)
    return(x)})



mapply(write.tree,
       ml_trees, 
       ml_treefiles)
