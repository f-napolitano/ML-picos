# script for separate all files in folder and subfolders into 2 sets, training
# and evaluation sets

library(dplyr)
library(caret)

osystem <- "windows"

filetrainpath <- ifelse(osystem == "ubuntu", 
                        "/media/federico/WD300_2/Fisica/Otros/Cesar/Sets_training/CR94-1200",
                        "f:/Fisica/Otros/Cesar/Sets_training/CR94-1200")

temp_df <- list.files(filetrainpath, pattern = "*.csv", recursive = FALSE)
DataTrainList <- lapply(temp_df, function(x) {
  drop_na(read.csv(paste0(filetrainpath, "/", x), stringsAsFactors = FALSE))
})

rm(temp_df)

seq_folders <- seq(1:length(DataTrainList))

test_index <- lapply(seq_folders, function(x) {
  createDataPartition(DataTrainList[[x]][[2]], times = 1, p = 0.49, list = FALSE)
})

train_set <- lapply(seq_folders, function(x) {
  DataTrainList[[x]][-test_index[[x]],]
})

test_set <- lapply(seq_folders, function(x) {
  DataTrainList[[x]][test_index[[x]],]
})

sum(sapply(train_set, nrow))
sum(sapply(test_set, nrow))

folders <- list.dirs(filetrainpath)[- 1]

train_set <- lapply(seq_folders, function(x) {
    colnames(train_set[[x]]) <- c('file', 'zona1', 'zona2', 'zona3', 'zona4', 'zona5')
    train_set[[x]] %>% mutate(set = folders[x])
})

train_set <- bind_rows(train_set)


test_set <- lapply(seq_folders, function(x) {
  colnames(test_set[[x]]) <- c('file', 'zona1', 'zona2', 'zona3', 'zona4', 'zona5')
  test_set[[x]] %>% mutate(set = folders[x]) #%>% mutate(fullfile = paste0(filetrainpath, test_set[[x]]$Archivo))
})

test_set <- bind_rows(test_set)

train_set <- train_set[-c(290), ]   #delete a random row to make train_set and test_set of the same length


write.csv(train_set, file = paste0(filetrainpath, "CR94-1200_train_set.csv"))
write.csv(test_set, file = paste0(filetrainpath, "CR94-1200_test_set.csv"))