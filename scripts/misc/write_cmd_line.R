treefiles <- list.files(path = './2024Jul12/region_beasttreefile',
                        pattern = 'combined10000')

shell <- paste0('/Applications/bin/treeannotator -burnin 1000000 -heights median ',
       treefiles,
       ' ',
       gsub('.trees$', '_mcc.tree', treefiles))

writeLines(shell,
            './2024Jul12/region_beasttreefile/treeannotator.sh')



xmlfiles <- list.files('./2024Aug18/region_subsampled_beasttraits/', pattern = '.xml')
treefiles <- list.files('./2024Aug18/region_subsampled_alignments/', pattern = '.xml')
s
shell <- paste0('python3 GetEmpiricalXML.py ',
                xmlfiles,
                ' ',
                unique(gsub('_\\d.xml$', '_1000.trees', treefiles)))

writeLines(shell,
           './2024Aug18/region_subsampled_beasttraits/getempirical.sh')