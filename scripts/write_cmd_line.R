treefiles <- list.files(path = './2024Jul12/region_beasttreefile',
                        pattern = 'combined10000')

shell <- paste0('/Applications/bin/treeannotator -burnin 1000000 -heights median ',
       treefiles,
       ' ',
       gsub('.trees$', '_mcc.tree', treefiles))

writeLines(shell,
            './2024Jul12/region_beasttreefile/treeannotator.sh')



treefiles <- list.files('./2024Aug18/reassortant_subsampled_beasttraits/')
s
shell <- paste0('python3 GetEmpiricalXML.py ',
                treefiles,
                ' ',
                gsub('_traits.xml$', '_relaxLn_constant_2000.trees', treefiles))

writeLines(shell,
           './2024Aug18/reassortant_subsampled_beasttraits/getempirical.sh')