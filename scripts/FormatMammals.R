################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            Regex functions to format host (mammals)
## Author:             James Baxter
## Date Created:       2025-08-01
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6, # Avoid scientific notation
  digits = 7 # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)


################################### DATA #######################################
# Read and inspect data

################################### MAIN #######################################
# Main analysis or transformation steps
FormatMammal <- function(x) {
  mammals <- c(
    "red fox" = "(red|blue|silver) fox|^fox$|vulpes vulpes",
    "brown rat" = "r{0,1}attus norvegicus",
    "harbour seal" = "harbou{0,1}r seal",
    "ferret" = "mustela furo|^ferret$",
    "european polecat" = "mustela putorius|european polecat",
    "beech marten" = "(stone|beech) marten",
    "fisher" = "pekania pennanti|^fisher$",
    "arctic fox" = "arc{0,1}tic[ -]fox|vulpes lagopus",
    "human" = "^humans{0,1}$|homo sapiens{0,1}",
    "lion" = "^lion$|leo panthera",
    "eared seal sp." = "^sea lion$",
    "brown bear" = "(brown|blue) bear|ursus arctos",
    "american black bear" = "^black bear$|ursus americanus",
    "true seal sp." = "^seal$|seal sp\\.",
    "atlantic grey seal" = "gr[ea]y seal",
    "bear sp." = "^bear$|bear sp\\.",
    "mustelid sp." = "^mink$|wild mink|mink sp\\.|polecat|^otter$|badger|white mink",
    "domestic ferret" = "^ferret$|mustela furo",
    "common otter" = "lutra {0,1}lutra|common otter",
    "aurochs" = "dairy cattle|bos taurus",
    "dolphin sp." = "^dolphin$",
    "common dolphin" = "(short[- ]beaked|long[- ]beaked) common dolphin|^common dolphin$",
    "porpoise sp." = "^porpoise$",
    "lynx" = "lynx sp\\.|^lynx$",
    "skunk sp." = "^skunk$|skunk sp\\.",
    "feline sp." = "^cat$|domestic cat|feline",
    "canid sp." = "canine",
    "murid sp." = "^murine$",
    "mammal sp," = "mammal",
    "cougar" = "mountain lion|cougar|puma concolor",
    "raccoon" = "raccoon",
    "bovid sp." = "bovine"
  )

  for (i in 1:length(mammals)) {
    if (any(grepl(mammals[[i]], x))) {
      x <- names(mammals)[[i]]
    }
  }

  return(x)
}

################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################
