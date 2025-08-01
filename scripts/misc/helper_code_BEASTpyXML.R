txt_files <- c(list.files('./2025Feb26/region_traits',
                        pattern = '.txt',
                        full.names = T),
               list.files('./2025Feb26/reassortant_traits',
                          pattern = '.txt',
                          full.names = T) )

csv <- lapply(txt_files, read_delim) %>%
  lapply(., function(x) x %>% select(tipnames, lat, long) %>%
           rename(taxon = tipnames))

csv_files <- gsub('txt$', 'csv', txt_files) 

mapply(write_csv,
       csv,
       csv_files)


traits <- list.files('./2025Feb26/empirical_tree_distributions',
           pattern = '.csv')

trees <- list.files('./2025Feb26/empirical_tree_distributions',
           pattern = '.trees')

fasta <- list.files('./2025Feb26/empirical_tree_distributions',
           pattern = '.fasta')

stem <- gsub('.*empirical_tree_distributions\\/|\\.fasta$', '', fasta) 

cmd <- c()
for (i in 1:length(traits)){
  cmd[i] <- paste0('python temp.py --fasta ', fasta[i],' --empirical-tree-model --empirical-tree-distribution ',trees[i],' --continuous-phylogeo --continuous-phylogeo-coords ',traits[i],' --chain-length 4000000 --log-every 400 --file-stem ' ,stem[i], '_empiricaltree_traits')
}

write_lines(cmd, './2025Feb26/empirical_tree_distributions/cmd.sh')
