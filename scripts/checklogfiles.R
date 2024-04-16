library(beastio)
library(tidyverse)


logfiles <- list.files('~/Downloads/temp',
                       full.names = T)


imported_logs <- lapply(logfiles, function(x) try(readLog(x)))


probs <- cbind.data.frame('filename' = rep(NA, length(imported_logs)),
                          'complete' = rep(NA, length(imported_logs)),
                          'problems' = rep(NA, length(imported_logs)))

for( i in 1:length(imported_logs)){
  print(i)
  probs[i,1] <- logfiles[[i]]
  probs[i,2] <- nrow(imported_logs[[i]]) > 9000
  
  lowESS <- try(checkESS(imported_logs[[i]]))
  if(length(lowESS)==0){
    probs[i,3] <- NA
  }else{
    lowESS_char <- paste0(lowESS %>% names(), ' (', lowESS, ")")
    probs[i,3] <- paste(lowESS_char, collapse = ', ')
  }
  
}