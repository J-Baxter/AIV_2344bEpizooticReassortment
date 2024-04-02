# Format Mammals
FormatMammal <- function(x){
  mammals <- c(
    "red fox" = "(red|blue|silver) fox|^fox$|vulpes vulpes",
    "brown rat" = "r{0,1}attus norvegicus",
    'harbor seal'= "harbou{0,1}r seal",
    'ferret'= "mustela furo|^ferret$",
    'european polecat'= 'mustela putorius|european polecat',
    'beech marten'= "(stone|beech) marten",
    "fisher"= "pekania pennanti|^fisher$",
    "arctic fox"= "arctic[ -]fox|vulpes lagopus",
    'human' = '^humans{0,1}$',
    'lion' = '^lion$|leo panthera',
    'eared seal sp.' = '^sea lion$',
    'true seal sp.' = '^seal$|seal sp\\.',
    'bear sp.'= "^bear$|bear sp\\.",
    'mustelid sp.'= '^mink$|wild mink|mink sp\\.|polecat|^otter$|badger',
    'common otter' = 'luta {0,1}lutra|common otter',
    "dolphin sp."= "^dolphin$",
    "porpoise sp."= "^porpoise$",
    "lynx"= "lynx sp\\.|^lynx$",
    'skunk sp.'= '^skunk$|skunk sp\\.',
    'feline sp.'= '^cat$|domestic cat|feline',
    'canine sp.'= 'canine'
    )
  
  for (i in 1:length(mammals)){
    if(any(grepl(mammals[[i]], x, fixed = TRUE))){
      x <- names(mammals)[[i]]
    }
  }
  
  return(x)
}

