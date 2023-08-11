# Spinoff of ML_picos2.R. Now we're separating the training stuff and the
# application of the already trained code to real data on two separate scripts.
# Here is the later one, the application into real data. It takes the best values
# calculated in ML_picos2.R and written in csv files.

library(tidyverse)
library(ggplot2)
library(data.table)
library(pracma)
library(dplyr)
library(caret)
library(lattice)
library(tibble)

osystem <- "windows"

# ----------- reading of the calculated parameters from training code ----------
filepath <- ifelse(osystem == "ubuntu", 
                   "/media/federico/WD300_2/Fisica/Otros/Cesar/training/", 
                   "f:/Fisica/Otros/Cesar/training/")

best_ispeak <- read.csv(file = paste(filepath, "best_ispeak.csv", sep = ""), header = TRUE)
model_comp <- read.csv(file = paste(filepath, "Maximum-Guess_performance.csv", sep = ""), header = TRUE)
validos_filt <- read.csv(file = paste(filepath, "validos_filt.csv", sep = ""), header = TRUE)
Fscore_beta <- read.csv(file = paste(filepath, "Fscore_beta.csv", sep = ""), header = TRUE)


peak_threshold <-  Fscore_beta$threshold1[which.max(Fscore_beta$maximum1)]


# ----------- importing data to be analized ------------------------------------
filepath <- ifelse(osystem == "ubuntu", 
                   "/media/federico/WD300_2/Fisica/Otros/Cesar/CR94-1200/CR94-1200_smooth/", 
                   "f:/Fisica/Otros/Cesar/CR94-1200/CR94-1200_smooth/")
fileEvalfiles <- read.table(file = paste(filepath, "CR94-1200_folders.txt", sep = ""))
fileEvalfiles <- rename(fileEvalfiles, fullfolder = V1) %>%
  mutate(namefolder = NA)

fileEvalfiles$namefolder <- substr(fileEvalfiles$fullfolder, 
                                   start = ifelse(osystem == "ubuntu", 71, 50),
                                   stop = str_length(fileEvalfiles$fullfolder) - 1)
fileEvalfiles$fullfolder <- paste(fileEvalfiles$fullfolder, "/", sep="")

filenames_full <- fileEvalfiles$fullfolder %>% 
  lapply(function(x) list.files(x, pattern = "*.txt", full.names = TRUE))  #list with all filenames and paths

lst_Eval_resultsM <- vector(mode = 'list', length = length(filenames_full)) #list with all results (empty for now)

# ------------ parameters ------------------------------------------------------
subzone_size <- 3/8    # as a fraction of validos_filt$min (min distance to another reflexion)

# -------------- Math on every dataset -----------------------------------------
folder_number <- seq(1:length(filenames_full))

ii <- 1

sapply(folder_number, function(ii) {
  
dataEvallist <- lapply(filenames_full[[ii]], function(x) read.csv(x))
seq_long_folder <- seq(1:length(dataEvallist))

# extract the "zones" of each diffractogram in folder ii
zona <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y){
    dataEvallist[[x]] %>% filter(V1 > (validos_filt$twothetaM[y] - validos_filt$min[y]/2) & 
                                 V1 < (validos_filt$twothetaM[y] + validos_filt$min[y]/2))
         })
})

print(sprintf("Done zones 1 to %i folder %i (total: %i diffractograms)", 5, ii, max(seq_long_folder)))

#extract the subzones of each diffractogram in folder ii where background will be calculated
subzona <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y){
    zona[[x]][[y]] %>% filter(V1 < (validos_filt$twothetaM[y] - validos_filt$min[y] * subzone_size) | 
                         V1 > (validos_filt$twothetaM[y] + validos_filt$min[y] * subzone_size)) 
      #select(V2)
         })
})


print(sprintf("Done subzones 1 to %i folder %i (total: %i diffractograms)", 5, ii, max(seq_long_folder)))

# compute maximum in each zone, sd in each subzone, for each diffractogram in folder ii
max_M <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y) {
    max(zona[[x]][[y]]$smooth)
  })
})

#sd_M <- lapply(seq_long_folder, function(x) {
#  lapply(seq(1:5), function(y) {
#    sd(subzona[[x]][[y]]$V2)
#  })
#})

sd_M_left <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y) {
    temp <- subzona[[x]][[y]] %>%
      filter(V1 < validos_filt$twothetaM[y]) %>%
      select(smooth)
    sd(temp$smooth)
  })
})

sd_M_right <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y) {
    temp <- subzona[[x]][[y]] %>% 
      filter(V1 > validos_filt$twothetaM[y]) %>%
      select(smooth)
    sd(temp$smooth)
  })
})

rm(temp)
#media_M <- lapply(seq_long_folder, function(x) {
#  lapply(seq(1:5), function(y) {
#    mean(subzona[[x]][[y]]$V2)
#  })
#})

sd_M <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y) {
    (sd_M_left[[x]][[y]] + sd_M_right[[x]][[y]]) * 0.5
  })
})

# calculate background in each zone and value of bkg at martensite peak position

bkground_coeff <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y) {
    fitbkground <- lm(subzona[[x]][[y]]$smooth ~ subzona[[x]][[y]]$V1)
    fitbkground$coef
  })
})

media_M <- lapply(seq_long_folder, function(x) {
  lapply(seq(1:5), function(y) {
    bkground_coeff[[x]][[y]][1] + bkground_coeff[[x]][[y]][2] * validos_filt$twothetaM[y]
  })
})

#for (x in seq_long_folder) {
#  zona[[1]][[5]] %>%  ggplot(aes(V1, smooth)) + geom_point(size = 3, alpha = .9, color = "grey") +
#    geom_line(aes(V1, bkground_coeff[[1]][[5]][1] + bkground_coeff[[1]][[5]][2] * V1), color = "red")
#}


print(sprintf("Done calculating MAX, MEAN and SD 1 to %i folder %i (total: %i diffractograms)", 5, ii, max(seq_long_folder)))
# now lets build the output dataframe with the results
df_Eval_resultsM <- data.frame("full_name" = filenames_full[[ii]])

temp_column <- c(rbind(paste0("zona",1:5,"_max"), 
                       paste0("zona",1:5,"_mean"), 
                       paste0("zona",1:5,"_sd"), 
                       paste0("prominence",1:5), 
                       paste0("peak",1:5)))

df_Eval_resultsM[,temp_column] <- NA
rm(temp_column)

temp <- transpose(as.data.frame(sapply(seq_long_folder, function(x) {
  sapply(1:5, function(y) {
    max_M[[x]][[y]]
  })
})))

df_Eval_resultsM[, (2+((1:5)-1)*5)] <- temp[,(1:5)]

temp <- transpose(as.data.frame(sapply(seq_long_folder, function(x) {
  sapply(1:5, function(y) {
    media_M[[x]][[y]]
  })
})))

df_Eval_resultsM[, (3+((1:5)-1)*5)] <- temp[, (1:5)]

temp <- transpose(as.data.frame(sapply(seq_long_folder, function(x) {
  sapply(1:5, function(y) {
    sd_M[[x]][[y]]
  })
})))

df_Eval_resultsM[, (4+((1:5)-1)*5)] <- temp[, (1:5)]

df_Eval_resultsM <- df_Eval_resultsM %>% mutate(prominence1 = (zona1_max - zona1_mean) / zona1_sd)

df_Eval_resultsM[, (5+((1:5)-1)*5)] <- (df_Eval_resultsM[, (2+((1:5)-1)*5)] - df_Eval_resultsM[, (3+((1:5)-1)*5)]) / df_Eval_resultsM[, (4+((1:5)-1)*5)]

df_Eval_resultsM[, (6+((1:5)-1)*5)] <- ifelse(df_Eval_resultsM[, (5+((1:5)-1)*5)] > peak_threshold, "TRUE", "FALSE")

lst_Eval_resultsM[[ii]] <- df_Eval_resultsM

write.csv(df_Eval_resultsM, paste(fileEvalfiles$fullfolder[ii], fileEvalfiles$namefolder[ii],".csv", sep=""), row.names=TRUE)

print(sprintf("Done writing output file folder %i", ii))

})

