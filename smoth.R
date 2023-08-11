library(tidyverse)
library(ggplot2)
library(dplyr)

filepath <- ifelse(osystem == "ubuntu", 
                   "/media/federico/WD300_2/Fisica/Otros/Cesar/CR94-1200/CR94-1200_evaluation/", 
                   "f:/Fisica/Otros/Cesar/CR94-1200/CR94-1200_evaluation/")
fileEvalfiles <- read.table(file = paste(filepath, "CR94-1200_folders.txt", sep = ""))
fileEvalfiles <- rename(fileEvalfiles, fullfolder = V1)
fileEvalfiles$fullfolder <- paste(fileEvalfiles$fullfolder, "/", sep="")
filenames_full <- fileEvalfiles$fullfolder %>% 
  lapply(function(x) list.files(x, pattern = "*.txt", full.names = TRUE))  #list with all filenames and paths

folder_number <- seq(1:length(filenames_full))

smooth_totalpoints <- 3173.99 / 4724 # tengo que chequar por que no funciona el que esta dentro del sapply

sapply(folder_number, function(ii) {
    
    smooth_file_origin <- lapply(filenames_full[[ii]], function(x) read.table(x))
    
    smooth_file_dest <- lapply(smooth_file_origin, function(x){
      #smooth_totalpoints <- diff(range(x))/nrow(x)
      smooth_fit <- loess(V2 ~ V1, degree = 1, span = 0.005*smooth_totalpoints, data = x)
      x %>% mutate(smooth = smooth_fit$fitted)
    } )
    
    dest_folder <- paste0("f:/Fisica/Otros/Cesar/CR94-1200/CR94-1200_smooth/", ii,"/")
    
    tail(strsplit(filenames_full[[1]][[1]], "/")[[1]], n = 1)
    
    sapply(seq(1:length(smooth_file_dest)), function(x) {
      write.csv(smooth_file_dest[[x]], paste(dest_folder, tail(strsplit(filenames_full[[ii]][[x]], "/")[[1]], n = 1), sep=""), row.names=FALSE, col.names = TRUE)
    })

})


#smooth_file_origin <- read.table("f:/Fisica/Otros/Cesar/Sets_training/CR94-1200/CR94_1200_P1HP_cyc1_20MPa_dyn-09098_dist698p6768_diffractogram.txt")
#
#plot(smooth_file_origin$V1, smooth_file_origin$V2)
#
#smooth_totalpoints <- diff(range(smooth_file_origin))/nrow(smooth_file_origin)
#
#
#smooth_fit <- loess(V2 ~ V1, degree=1, span = 0.005*smooth_totalpoints, data=smooth_file_origin)
#smooth_file_origin <- smooth_file_origin %>% mutate(smooth = smooth_fit$fitted)
#smooth_file_origin %>% ggplot(aes(V1, V2)) +
#  geom_point(size = 3, alpha = .5, color = "grey") +
#  geom_line(aes(V1, smooth), color = "red")
#
#write.table(smooth_file_origin, file = "f:/Fisica/Otros/Cesar/Sets_training/CR94-1200/prueba_smooth.dat", row.names = FALSE, col.names = FALSE)
