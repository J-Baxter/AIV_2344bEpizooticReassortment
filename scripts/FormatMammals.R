# Format Mammals


#grepl('^otter$', source) & grepl('europe', region) ~ 'eurasian otter',
#grepl('badger', source) & grepl('europe', region) ~ 'european badger'
# Needs mammal csv

FormatMammals <- function(x){
  mammals <- c(
    "red fox" = "(red|blue|silver) fox|^fox$|vulpes vulpes",
    "brown rat" = "r{0,1}attus norvegicus",
    "seal sp."= '^seal$|^sea lion$|seal sp\\.',
    'mink sp.'= '^mink$|wild mink|mink sp\\.',
    'harbour seal'= "harbou{0,1}r seal",
    'ferret'= "mustela furo|^ferret$",
    'european polecat'= 'mustela putorius|european polecat',
    'beech marten'= "(stone|beech) marten",
    'bear sp.'= "^bear$|bear sp\\.",
    "cetacea sp."= "dolphin$|porpoise|cetacea sp\\.",
    "fisher"= "pekania pennanti",
    "lynx sp."= "lynx sp\\.|^lynx$",
    "mustelid sp."= "^polecat$",
    'skunk sp.'= '^skunk$|skunk sp\\.',
    'otter sp.'= '^otter$|lutralutra|otter sp\\.',
    'feline sp.'= '^lion$|^cat$|domestic cat|feline|lion|feline sp\\.',
    'canine sp.'= 'canine',
    "arctic fox"= "arctic[ -]fox|vulpes lagopus",
    'human' = '^humans{0,1}$')
  
  for (i in 1:length(mammals)){
    if(any(grepl(mammals[[i]], x, fixed = TRUE))){
      x <- names(mammals)[[i]]
    }
  }
  
  return(x)
}