# Script to combine all Matlab files into one database
# Blair Shevlin, May 2023

library(tidyverse)

path = "C:/Users/blair/Documents/Research/Sinai-UG-ED/ugfit"
out = "C:/Users/blair/Documents/Research/Sinai-UG-ED/"
files = list.files(path = path, pattern = "*.csv",)

df <- data.frame(IDs = character(),
                 BIC = numeric(),
                 group = character(),
                 condition = character(),
                 model = character(),
                 alpha = numeric(),
                 beta = numeric(),
                 epsilon = numeric(),
                 delta = numeric())
for (i in files) {
  tmp = read.csv(paste0(path,"/",i))
  
  # Figure out which group
  if (grepl("BN",i)) {
    tmp$group = "BN"
    
  } else if (grepl("BED",i)) {
    tmp$group = "BED"
    
  } else if (grepl("BNHC",i)) {
    tmp$group = "BNHC"
    
  } else if (grepl("BEDHC",i)) {
    tmp$group = "BEDHC"
  }
  
  if (grepl("IC",i)) {
    tmp$condition = "IC"
  } else if (grepl("NC",i)) {
    tmp$condition = "NC"
  }
  
  # Figure out which model
  if (grepl("UG0",i)){
    tmp <- tmp %>%
      mutate(model = "UG0",
             alpha = params_1,
             beta = params_2,
             epsilon = params_3,
             delta = NA) %>%
      select(!c(params_1,params_2,params_3))
    
  } else if (grepl("UG1",i)) {
    tmp <- tmp %>%
      mutate(model = "UG1",
             alpha = params_1,
             beta = params_2,
             epsilon = params_3,
             delta = params_4) %>%
      select(!c(params_1,params_2,params_3,params_4))
    
  } else if (grepl("UG2",i)) {
    tmp <- tmp %>%
      mutate(model = "UG2",
             alpha = params_1,
             beta = params_2,
             epsilon = params_3,
             delta = params_4) %>%
      select(!c(params_1,params_2,params_3,params_4))
    
  } else if (grepl("UG3",i)) {
    tmp <- tmp %>%
      mutate(model = "UG3",
             alpha = params_1,
             beta = params_2,
             epsilon = params_3,
             delta = params_4) %>%
      select(!c(params_1,params_2,params_3,params_4))
  }
  
  df <- rbind(df,tmp)
}

write.csv(df,paste0(out,"FULL_WEIGHTED_FIT_ED_May24_2023.csv"),row.names = F)











