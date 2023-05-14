# Importing all csv files from a folder, joining them, and exporting it as a single csv
library(dplyr)
library(tidyverse)
library(readr)

osystem <- "windows"
pathtrainfiles <- ifelse(osystem == "ubuntu", "/media/federico/WD300_2/Fisica/Otros/Cesar/training" ,"f:/Fisica/Otros/Cesar/training")

datatrainlist <- list.files(path = pathtrainfiles)
datatrainlist <- file.path(pathtrainfiles, datatrainlist)

datatrainfr <- do.call("rbind",lapply(datatrainlist,FUN=function(files){ read.csv(files)}))

datatrainfr <- datatrainfr %>% relocate(set, .after = Archivo)  #dataframe with all csv's info concatenated

# new column with full path and filename of each train file
datatrainfr <- datatrainfr %>% mutate(fullfile = paste('f:/Fisica/Otros/Cesar/Sets/',set,'/','CR94_1300_P5_',Archivo,'_dist698p6768_diffractogram.txt', sep = ""))

# exporting to csv
write_csv(datatrainfr, path = paste(pathtrainfiles,'/','training_set.csv', sep = ""))
