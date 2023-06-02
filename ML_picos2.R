# Variation of ML_picos1.R now already knowing if there are peaks in 10% of the whole data sets
# here we're going to take this 10% and separate them in train and evaluation sets (50/50, random)
# and make an algorithm that calculates existance of peaks in train sets and compare it with was manually done

library(tidyverse)
library(ggplot2)
library(data.table)
library(pracma)
library(dplyr)
library(caret)
library(lattice)
library(tibble)

osystem <- "windows"

#----------- Reflection determination as predictor -------------
# In this block we're going to find where the martensite reflection are in a calculated diffractogram
# but looking only those that are reasonable far away fromm those of the austenite and gamma phases
# reflections found in this way are the one that are going to be used as predictors

filepath <- ifelse(osystem == "ubuntu", "/media/federico/WD300_2/Fisica/Otros/Cesar/SimulacionXRD/" ,"f:/Fisica/Otros/Cesar/SimulacionXRD/")
marten <- read.csv(file = paste(filepath, "Martensita.csv", sep = ""), header = TRUE) %>% select(2:4)
aust <- read.csv(file = paste(filepath, "Austenita.csv", sep = ""), header = TRUE) %>% select(2:4)
gammaphase <- read.csv(file = paste(filepath, "Gamma.csv", sep = ""), header = TRUE) %>% select(2:4)

validos <- marten

# calculating distances between nearest (austenite, gamma) reflections to each one of martensite
diff_aust_validos <- sapply(validos$twothetaM, function(x) aust$twothetaA[ which.min(abs(x - aust$twothetaA))] - x)
validos <- validos %>% mutate(diff_aust = diff_aust_validos)
rm(diff_aust_validos)
diff_gamma_validos <- sapply(validos$twothetaM, function(x) gammaphase$twothetaG[ which.min((abs(x - gammaphase$twothetaG)))] - x)
validos <- validos %>% mutate(diff_gamma = diff_gamma_validos)
rm(diff_gamma_validos)

validos <- validos %>% rowwise() %>% mutate(min = min(abs(diff_aust), abs(diff_gamma))) %>% 
  filter(abs(diff_aust) > 0.1 & abs(diff_gamma) > 0.1) %>% 
  mutate(distancia = sqrt(diff_aust**2 + diff_gamma**2))

validos %>% ggplot() + aes(abs(diff_aust), abs(diff_gamma), label = round(twothetaM, 2)) + geom_point() + geom_text(hjust=0, vjust=-0.5)
validos %>% ggplot() + aes(twothetaM, distancia) + geom_point()

distancia_min <- 0.5 # a measure of the min distance required between martensite and other phases peaks to be considered
validos_filt <- validos %>% filter(distancia > distancia_min & abs(diff_gamma) > 0.25 & abs(diff_aust) > 0.25)

#----------- peaks position manual check and fine tunning of their position by hand --------
validos_filt[3,1] <- 8.70


#---------- parameters declaration ---------------------------
ispeak <- 4 # how much noise's sd has to have a peak maximum to be considered as a TRUE peak

# -------------- creating train and test set --------------
filetrainpath <- ifelse(osystem == "ubuntu", "/media/federico/WD300_2/Fisica/Otros/Cesar/training/" ,"f:/Fisica/Otros/Cesar/training/")
selected_sets <- read.csv(file = paste(filetrainpath, "training_set.csv", sep = ""), header = TRUE)

test_index <- createDataPartition(selected_sets$set, times = 1, p = 0.489, list = FALSE)

train_set <- selected_sets %>% slice(-test_index)
test_set <- selected_sets %>% slice(test_index)

# --------------- Importing files from train set and dataframes creation -------------
datatrainlist <- lapply(as.character(train_set$fullfile), function(x)read.table(x))
datatrainfr <- bind_cols(datatrainlist)

Output_train_list <- train_set$Archivo
Output_train_list <- as.data.frame(Output_train_list) %>% mutate(set = train_set$set)
Output_train_list <- rename(Output_train_list, filename = Output_train_list)

diff1 <- datatrainfr

#zona_i definen un rango de datos alrededor de los picos de la martensita definidos en "validos_filt" de forma que no agarre pico de otra fase
zona1 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[1] - validos_filt$min[1] /2) & V1...1 < (validos_filt$twothetaM[1] + validos_filt$min[1] /2))
zona2 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[2] - validos_filt$min[2] /2) & V1...1 < (validos_filt$twothetaM[2] + validos_filt$min[2] /2))
zona3 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[3] - validos_filt$min[3] /2) & V1...1 < (validos_filt$twothetaM[3] + validos_filt$min[3] /2))
zona4 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[4] - validos_filt$min[4] /2) & V1...1 < (validos_filt$twothetaM[4] + validos_filt$min[4] /2))
zona5 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[5] - validos_filt$min[5] /2) & V1...1 < (validos_filt$twothetaM[5] + validos_filt$min[5] /2))

# ------------- calculo de sd y mean ---------------------------
# descarto la region donde deberia estar el pico de la reflexion martensita, me quedo con los bordes para calcular mean y sd
intsd1 <- zona1 %>% filter(V1...1 < (validos_filt$twothetaM[1] - validos_filt$min[1] /4) | V1...1 > (validos_filt$twothetaM[1] + validos_filt$min[1] /4)) 
intsd1 <- intsd1[, c(FALSE, TRUE)] 
intsd2 <- zona2 %>% filter(V1...1 < (validos_filt$twothetaM[2] - validos_filt$min[2] /4) | V1...1 > (validos_filt$twothetaM[2] + validos_filt$min[2] /4))
intsd2 <- intsd2[, c(FALSE, TRUE)]
intsd3 <- zona3 %>% filter(V1...1 < (validos_filt$twothetaM[3] - validos_filt$min[3] /4) | V1...1 > (validos_filt$twothetaM[3] + validos_filt$min[3] /4))
intsd3 <- intsd3[, c(FALSE, TRUE)]
intsd4 <- zona4 %>% filter(V1...1 < (validos_filt$twothetaM[4] - validos_filt$min[4] /4) | V1...1 > (validos_filt$twothetaM[4] + validos_filt$min[4] /4))
intsd4 <- intsd4[, c(FALSE, TRUE)]
intsd5 <- zona5 %>% filter(V1...1 < (validos_filt$twothetaM[5] - validos_filt$min[5] /4) | V1...1 > (validos_filt$twothetaM[5] + validos_filt$min[5] /4))
intsd5 <- intsd5[, c(FALSE, TRUE)]

Output_train_list <- Output_train_list %>% mutate(sd1 = transpose(data.frame(lapply(intsd1, function (x) sd(x)))), 
                                         media1 = transpose(data.frame(lapply(intsd1, function (x) mean(x)))), 
                                         sd2 = transpose(data.frame(lapply(intsd2, function (x) sd(x)))), 
                                         media2 = transpose(data.frame(lapply(intsd2, function (x) mean(x)))), 
                                         sd3 = transpose(data.frame(lapply(intsd3, function (x) sd(x)))), 
                                         media3 = transpose(data.frame(lapply(intsd3, function (x) mean(x)))), 
                                         sd4 = transpose(data.frame(lapply(intsd4, function (x) sd(x)))), 
                                         media4 = transpose(data.frame(lapply(intsd4, function (x) mean(x)))), 
                                         media5 = transpose(data.frame(lapply(intsd5, function (x) mean(x)))), 
                                         sd5 = transpose(data.frame(lapply(intsd5, function (x) sd(x)))))

#--------------- Peak determination ------------------------
intsd1 <- zona1[, c(FALSE, TRUE)]
intsd2 <- zona2[, c(FALSE, TRUE)]
intsd3 <- zona3[, c(FALSE, TRUE)]
intsd4 <- zona4[, c(FALSE, TRUE)]
intsd5 <- zona5[, c(FALSE, TRUE)]

# first method: searching the maximum and compare it with the sd of the boundaries in each zone
Output_train_list <- Output_train_list %>% mutate(Maximo1 = sapply(intsd1, function (x) max(x)), 
                                                  pico1 = diag( sapply(intsd1, function(x) { ifelse(max(x) > media1 + ispeak * sd1, TRUE, FALSE)} )), 
                                                  Maximo2 = sapply(intsd2, function(x) max(x)),
                                                  pico2 = diag( sapply(intsd2, function(x) { ifelse(max(x) > media2 + ispeak * sd2, TRUE, FALSE)} )),
                                                  Maximo3 = sapply(intsd3, function(x) max(x)),
                                                  pico3 = diag( sapply(intsd3, function(x) { ifelse(max(x) > media3 + ispeak * sd3, TRUE, FALSE)} )),
                                                  Maximo4 = sapply(intsd4, function(x) max(x)),
                                                  pico4 = diag( sapply(intsd4, function(x) { ifelse(max(x) > media4 + ispeak * sd4, TRUE, FALSE)} )),
                                                  Maximo5 = sapply(intsd5, function(x) max(x)),
                                                  pico5 = diag( sapply(intsd5, function(x) { ifelse(max(x) > media5 + ispeak * sd5, TRUE, FALSE)} )))

# second method: fitting a gaussian inside each zone, took their heights and compare them with de sd of the boundaries

fitG =
  function(x,y,mu,sig,scale,background){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2]) + background
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }

Output_train_list[, c('fitGmu1', 'fitGsigma1', 'fitGscale1', 'MaxGauss1', 'peakG1',
                      'fitGmu2', 'fitGsigma2', 'fitGscale2', 'MaxGauss2', 'peakG2',
                      'fitGmu3', 'fitGsigma3', 'fitGscale3', 'MaxGauss3', 'peakG3',
                      'fitGmu4', 'fitGsigma4', 'fitGscale4', 'MaxGauss4', 'peakG4',
                      'fitGmu5', 'fitGsigma5', 'fitGscale5', 'MaxGauss5', 'peakG5')
                  ] <- NA

# calculate the gaussian fit in each zone
for(i in 1:round(ncol(diff1)/2)) {
  
  fitGauss <- as.list(fitG(zona1[[2*i-1]], zona1[[2*i]], validos_filt$twothetaM[1], 1, 1, Output_train_list$media1[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_train_list$fitGmu1[i] <- as.numeric(peakpar[1])
  Output_train_list$fitGsigma1[i] <- as.numeric(peakpar[2])
  Output_train_list$fitGscale1[i] <- as.numeric(peakpar[3])
  
#  plot(zona1[[2*i-1]], zona1[[2*i]], main = Output_train_list$Output_train_list[i])
#  lines(zona1[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona1[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media1[[1]][i], col = "red", pch = 20)
  
  
  fitGauss <- as.list(fitG(zona2[[2*i-1]], zona2[[2*i]], validos_filt$twothetaM[2], 1, 1, Output_train_list$media2[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_train_list$fitGmu2[i] <- as.numeric(peakpar[1])
  Output_train_list$fitGsigma2[i] <- as.numeric(peakpar[2])
  Output_train_list$fitGscale2[i] <- as.numeric(peakpar[3])
  
#  plot(zona2[[2*i-1]], zona2[[2*i]], main = Output_train_list$Output_train_list[i])
#  lines(zona2[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona2[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media2[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona3[[2*i-1]], zona3[[2*i]], validos_filt$twothetaM[3], 1, 1, Output_train_list$media3[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_train_list$fitGmu3[i] <- as.numeric(peakpar[1])
  Output_train_list$fitGsigma3[i] <- as.numeric(peakpar[2])
  Output_train_list$fitGscale3[i] <- as.numeric(peakpar[3])
  
#  plot(zona3[[2*i-1]], zona3[[2*i]], main = Output_train_list$Output_train_list[i])
#  lines(zona3[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona3[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media3[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona4[[2*i-1]], zona4[[2*i]], validos_filt$twothetaM[4], 1, 1, Output_train_list$media4[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_train_list$fitGmu4[i] <- as.numeric(peakpar[1])
  Output_train_list$fitGsigma4[i] <- as.numeric(peakpar[2])
  Output_train_list$fitGscale4[i] <- as.numeric(peakpar[3])
  
#  plot(zona4[[2*i-1]], zona4[[2*i]], main = Output_train_list$Output_train_list[i])
#  lines(zona4[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona4[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media4[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona5[[2*i-1]], zona5[[2*i]], validos_filt$twothetaM[5], 1, 1, Output_train_list$media5[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_train_list$fitGmu5[i] <- as.numeric(peakpar[1])
  Output_train_list$fitGsigma5[i] <- as.numeric(peakpar[2])
  Output_train_list$fitGscale5[i] <- as.numeric(peakpar[3])
  
#  plot(zona5[[2*i-1]], zona5[[2*i]], main = Output_train_list$Output_train_list[i])
#  lines(zona5[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona5[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media5[[1]][i], col = "red", pch = 20)

}

# calculate the height of each gaussian fit, compare them with sd and decide if its a peak or not
Output_train_list$MaxGauss1 <- Output_train_list$fitGscale1 * dnorm(Output_train_list$fitGmu1, Output_train_list$fitGmu1, Output_train_list$fitGsigma1) + Output_train_list$media1
Output_train_list$peakG1 <- ifelse(Output_train_list$MaxGauss1 > Output_train_list$media1 + ispeak * Output_train_list$sd1, TRUE, FALSE)
Output_train_list$MaxGauss2 <- Output_train_list$fitGscale2 * dnorm(Output_train_list$fitGmu2, Output_train_list$fitGmu2, Output_train_list$fitGsigma2) + Output_train_list$media2
Output_train_list$peakG2 <- ifelse(Output_train_list$MaxGauss2 > Output_train_list$media2 + ispeak * Output_train_list$sd2, TRUE, FALSE)
Output_train_list$MaxGauss3 <- Output_train_list$fitGscale3 * dnorm(Output_train_list$fitGmu3, Output_train_list$fitGmu3, Output_train_list$fitGsigma3) + Output_train_list$media3
Output_train_list$peakG3 <- ifelse(Output_train_list$MaxGauss3 > Output_train_list$media3 + ispeak * Output_train_list$sd3, TRUE, FALSE)
Output_train_list$MaxGauss4 <- Output_train_list$fitGscale4 * dnorm(Output_train_list$fitGmu4, Output_train_list$fitGmu4, Output_train_list$fitGsigma4) + Output_train_list$media4
Output_train_list$peakG4 <- ifelse(Output_train_list$MaxGauss4 > Output_train_list$media4 + ispeak * Output_train_list$sd4, TRUE, FALSE)
Output_train_list$MaxGauss5 <- Output_train_list$fitGscale5 * dnorm(Output_train_list$fitGmu5, Output_train_list$fitGmu5, Output_train_list$fitGsigma5) + Output_train_list$media5
Output_train_list$peakG5 <- ifelse(Output_train_list$MaxGauss5 > Output_train_list$media5 + ispeak * Output_train_list$sd5, TRUE, FALSE)


# -------------- creating predictors dataframe ---------------------
# notice that predictors were already calculated just before but then algorithm directly put TRUE/FALSE instead of a number
# evaluate to make it there with the number and in the new dataframe here already calculate the TRUE/FALSE thing

# PredPicoM# and PredPicoG# will be how much the "peak height" will be greater than sd in each region (for both methods)
Output_train_list <- add_column(Output_train_list, PredPicoM1 = NA, .after = "Maximo1")
Output_train_list <- add_column(Output_train_list, PredPicoM2 = NA, .after = "Maximo2")
Output_train_list <- add_column(Output_train_list, PredPicoM3 = NA, .after = "Maximo3")
Output_train_list <- add_column(Output_train_list, PredPicoM4 = NA, .after = "Maximo4")
Output_train_list <- add_column(Output_train_list, PredPicoM5 = NA, .after = "Maximo5")

Output_train_list <- add_column(Output_train_list, PredPicoG1 = NA, .after = "MaxGauss1")
Output_train_list <- add_column(Output_train_list, PredPicoG2 = NA, .after = "MaxGauss2")
Output_train_list <- add_column(Output_train_list, PredPicoG3 = NA, .after = "MaxGauss3")
Output_train_list <- add_column(Output_train_list, PredPicoG4 = NA, .after = "MaxGauss4")
Output_train_list <- add_column(Output_train_list, PredPicoG5 = NA, .after = "MaxGauss5")

Output_train_list$PredPicoM1 <- (Output_train_list$Maximo1 - Output_train_list$media1) / Output_train_list$sd1
Output_train_list$PredPicoM2 <- (Output_train_list$Maximo2 - Output_train_list$media2) / Output_train_list$sd2
Output_train_list$PredPicoM3 <- (Output_train_list$Maximo3 - Output_train_list$media3) / Output_train_list$sd3
Output_train_list$PredPicoM4 <- (Output_train_list$Maximo4 - Output_train_list$media4) / Output_train_list$sd4
Output_train_list$PredPicoM5 <- (Output_train_list$Maximo5 - Output_train_list$media5) / Output_train_list$sd5

Output_train_list$PredPicoG1 <- (Output_train_list$MaxGauss1 - Output_train_list$media1) / Output_train_list$sd1
Output_train_list$PredPicoG2 <- (Output_train_list$MaxGauss2 - Output_train_list$media2) / Output_train_list$sd2
Output_train_list$PredPicoG3 <- (Output_train_list$MaxGauss3 - Output_train_list$media3) / Output_train_list$sd3
Output_train_list$PredPicoG4 <- (Output_train_list$MaxGauss4 - Output_train_list$media4) / Output_train_list$sd4
Output_train_list$PredPicoG5 <- (Output_train_list$MaxGauss5 - Output_train_list$media5) / Output_train_list$sd5

Output_train_predM_df <- Output_train_list %>% select(c('filename', 'set',
                                                        'PredPicoM1', 'pico1', 
                                                        'PredPicoM2', 'pico2',
                                                        'PredPicoM3', 'pico3',
                                                        'PredPicoM4', 'pico4',
                                                        'PredPicoM5', 'pico5'))
Output_train_predG_df <- Output_train_list %>% select(c('filename', 'set',
                                                        'PredPicoG1', 'peakG1', 
                                                        'PredPicoG2', 'peakG2',
                                                        'PredPicoG3', 'peakG3',
                                                        'PredPicoG4', 'peakG4',
                                                        'PredPicoG5', 'peakG5'))

# ------------- building the algorithm on train set ----------------------
x <- Output_train_predM_df$PredPicoM1
y <- as.character(train_set$zona.1)
test_set$zona.1 <- as.character(test_set$zona.1)
y_hat <- ifelse(x < ispeak, "FALSE", "TRUE")
y_hat <- as.factor(y_hat)


mean(y == y_hat) 
confusionMatrix(data = y_hat, reference = as.factor(test_set$zona.1))
# here we have the evidence of what was clear during exploration, data has a huge prevalence of NO PEAK outcomes (here counted as positive outcome)
# So clearly we can't use accuracy, need to switch to F-1 and balanced accuracy

# maximizing F-score
ispeak <- seq(from = 1.75, to = 6, by = 0.1)
F_1 <- map_dbl(ispeak, function(x){ 
  y_hat <- ifelse(Output_train_predM_df$PredPicoM1 < x, "FALSE", "TRUE")
  y_hat <- as.factor(y_hat)
  F_meas(data = y_hat, reference = as.factor(train_set$zona.1))
  })

ggplot() + aes(ispeak, F_1) + geom_point() + geom_line()
max(F_1)
best_ispeak <- ispeak[which.max(F_1)]
best_ispeak
# DONE!!


# ------------ construction "test set" as above with "train set" -----------------
datatestlist <- lapply(as.character(test_set$fullfile), function(x)read.table(x))
datatestfr <- bind_cols(datatestlist)

Output_test_list <- test_set$Archivo
Output_test_list <- as.data.frame(Output_test_list) %>% mutate(set = test_set$set)
Output_test_list <- rename(Output_test_list, filename = Output_test_list)

diff1 <- datatestfr

zona1 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[1] - validos_filt$min[1] /2) & V1...1 < (validos_filt$twothetaM[1] + validos_filt$min[1] /2))
zona2 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[2] - validos_filt$min[2] /2) & V1...1 < (validos_filt$twothetaM[2] + validos_filt$min[2] /2))
zona3 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[3] - validos_filt$min[3] /2) & V1...1 < (validos_filt$twothetaM[3] + validos_filt$min[3] /2))
zona4 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[4] - validos_filt$min[4] /2) & V1...1 < (validos_filt$twothetaM[4] + validos_filt$min[4] /2))
zona5 <- diff1 %>% filter(V1...1 > (validos_filt$twothetaM[5] - validos_filt$min[5] /2) & V1...1 < (validos_filt$twothetaM[5] + validos_filt$min[5] /2))
intsd1 <- zona1 %>% filter(V1...1 < (validos_filt$twothetaM[1] - validos_filt$min[1] /4) | V1...1 > (validos_filt$twothetaM[1] + validos_filt$min[1] /4)) 
intsd1 <- intsd1[, c(FALSE, TRUE)] 
intsd2 <- zona2 %>% filter(V1...1 < (validos_filt$twothetaM[2] - validos_filt$min[2] /4) | V1...1 > (validos_filt$twothetaM[2] + validos_filt$min[2] /4))
intsd2 <- intsd2[, c(FALSE, TRUE)]
intsd3 <- zona3 %>% filter(V1...1 < (validos_filt$twothetaM[3] - validos_filt$min[3] /4) | V1...1 > (validos_filt$twothetaM[3] + validos_filt$min[3] /4))
intsd3 <- intsd3[, c(FALSE, TRUE)]
intsd4 <- zona4 %>% filter(V1...1 < (validos_filt$twothetaM[4] - validos_filt$min[4] /4) | V1...1 > (validos_filt$twothetaM[4] + validos_filt$min[4] /4))
intsd4 <- intsd4[, c(FALSE, TRUE)]
intsd5 <- zona5 %>% filter(V1...1 < (validos_filt$twothetaM[5] - validos_filt$min[5] /4) | V1...1 > (validos_filt$twothetaM[5] + validos_filt$min[5] /4))
intsd5 <- intsd5[, c(FALSE, TRUE)]


Output_test_list <- Output_test_list %>% mutate(sd1 = transpose(data.frame(lapply(intsd1, function (x) sd(x)))), 
                                                  media1 = transpose(data.frame(lapply(intsd1, function (x) mean(x)))), 
                                                  sd2 = transpose(data.frame(lapply(intsd2, function (x) sd(x)))), 
                                                  media2 = transpose(data.frame(lapply(intsd2, function (x) mean(x)))), 
                                                  sd3 = transpose(data.frame(lapply(intsd3, function (x) sd(x)))), 
                                                  media3 = transpose(data.frame(lapply(intsd3, function (x) mean(x)))), 
                                                  sd4 = transpose(data.frame(lapply(intsd4, function (x) sd(x)))), 
                                                  media4 = transpose(data.frame(lapply(intsd4, function (x) mean(x)))), 
                                                  media5 = transpose(data.frame(lapply(intsd5, function (x) mean(x)))), 
                                                  sd5 = transpose(data.frame(lapply(intsd5, function (x) sd(x)))))
intsd1 <- zona1[, c(FALSE, TRUE)]
intsd2 <- zona2[, c(FALSE, TRUE)]
intsd3 <- zona3[, c(FALSE, TRUE)]
intsd4 <- zona4[, c(FALSE, TRUE)]
intsd5 <- zona5[, c(FALSE, TRUE)]

Output_test_list <- Output_test_list %>% mutate(Maximo1 = sapply(intsd1, function (x) max(x)), 
                                                  pico1 = diag( sapply(intsd1, function(x) { ifelse(max(x) > media1 + ispeak * sd1, TRUE, FALSE)} )), 
                                                  Maximo2 = sapply(intsd2, function(x) max(x)),
                                                  pico2 = diag( sapply(intsd2, function(x) { ifelse(max(x) > media2 + ispeak * sd2, TRUE, FALSE)} )),
                                                  Maximo3 = sapply(intsd3, function(x) max(x)),
                                                  pico3 = diag( sapply(intsd3, function(x) { ifelse(max(x) > media3 + ispeak * sd3, TRUE, FALSE)} )),
                                                  Maximo4 = sapply(intsd4, function(x) max(x)),
                                                  pico4 = diag( sapply(intsd4, function(x) { ifelse(max(x) > media4 + ispeak * sd4, TRUE, FALSE)} )),
                                                  Maximo5 = sapply(intsd5, function(x) max(x)),
                                                  pico5 = diag( sapply(intsd5, function(x) { ifelse(max(x) > media5 + ispeak * sd5, TRUE, FALSE)} )))
Output_test_list[, c('fitGmu1', 'fitGsigma1', 'fitGscale1', 'MaxGauss1', 'peakG1',
                      'fitGmu2', 'fitGsigma2', 'fitGscale2', 'MaxGauss2', 'peakG2',
                      'fitGmu3', 'fitGsigma3', 'fitGscale3', 'MaxGauss3', 'peakG3',
                      'fitGmu4', 'fitGsigma4', 'fitGscale4', 'MaxGauss4', 'peakG4',
                      'fitGmu5', 'fitGsigma5', 'fitGscale5', 'MaxGauss5', 'peakG5')
] <- NA

for(i in 1:round(ncol(diff1)/2)) {
  
  fitGauss <- as.list(fitG(zona1[[2*i-1]], zona1[[2*i]], validos_filt$twothetaM[1], 1, 1, Output_test_list$media1[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_test_list$fitGmu1[i] <- as.numeric(peakpar[1])
  Output_test_list$fitGsigma1[i] <- as.numeric(peakpar[2])
  Output_test_list$fitGscale1[i] <- as.numeric(peakpar[3])
  
  #  plot(zona1[[2*i-1]], zona1[[2*i]], main = Output_train_list$Output_train_list[i])
  #  lines(zona1[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona1[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media1[[1]][i], col = "red", pch = 20)
  
  
  fitGauss <- as.list(fitG(zona2[[2*i-1]], zona2[[2*i]], validos_filt$twothetaM[2], 1, 1, Output_test_list$media2[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_test_list$fitGmu2[i] <- as.numeric(peakpar[1])
  Output_test_list$fitGsigma2[i] <- as.numeric(peakpar[2])
  Output_test_list$fitGscale2[i] <- as.numeric(peakpar[3])
  
  #  plot(zona2[[2*i-1]], zona2[[2*i]], main = Output_train_list$Output_train_list[i])
  #  lines(zona2[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona2[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media2[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona3[[2*i-1]], zona3[[2*i]], validos_filt$twothetaM[3], 1, 1, Output_test_list$media3[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_test_list$fitGmu3[i] <- as.numeric(peakpar[1])
  Output_test_list$fitGsigma3[i] <- as.numeric(peakpar[2])
  Output_test_list$fitGscale3[i] <- as.numeric(peakpar[3])
  
  #  plot(zona3[[2*i-1]], zona3[[2*i]], main = Output_train_list$Output_train_list[i])
  #  lines(zona3[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona3[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media3[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona4[[2*i-1]], zona4[[2*i]], validos_filt$twothetaM[4], 1, 1, Output_test_list$media4[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_test_list$fitGmu4[i] <- as.numeric(peakpar[1])
  Output_test_list$fitGsigma4[i] <- as.numeric(peakpar[2])
  Output_test_list$fitGscale4[i] <- as.numeric(peakpar[3])
  
  #  plot(zona4[[2*i-1]], zona4[[2*i]], main = Output_train_list$Output_train_list[i])
  #  lines(zona4[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona4[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media4[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona5[[2*i-1]], zona5[[2*i]], validos_filt$twothetaM[5], 1, 1, Output_test_list$media5[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Output_test_list$fitGmu5[i] <- as.numeric(peakpar[1])
  Output_test_list$fitGsigma5[i] <- as.numeric(peakpar[2])
  Output_test_list$fitGscale5[i] <- as.numeric(peakpar[3])
  
  #  plot(zona5[[2*i-1]], zona5[[2*i]], main = Output_train_list$Output_train_list[i])
  #  lines(zona5[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona5[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Output_train_list$media5[[1]][i], col = "red", pch = 20)
  
}

Output_test_list$MaxGauss1 <- Output_test_list$fitGscale1 * dnorm(Output_test_list$fitGmu1, Output_test_list$fitGmu1, Output_test_list$fitGsigma1) + Output_test_list$media1
Output_test_list$peakG1 <- ifelse(Output_test_list$MaxGauss1 > Output_test_list$media1 + ispeak * Output_test_list$sd1, TRUE, FALSE)
Output_test_list$MaxGauss2 <- Output_test_list$fitGscale2 * dnorm(Output_test_list$fitGmu2, Output_test_list$fitGmu2, Output_test_list$fitGsigma2) + Output_test_list$media2
Output_test_list$peakG2 <- ifelse(Output_test_list$MaxGauss2 > Output_test_list$media2 + ispeak * Output_test_list$sd2, TRUE, FALSE)
Output_test_list$MaxGauss3 <- Output_test_list$fitGscale3 * dnorm(Output_test_list$fitGmu3, Output_test_list$fitGmu3, Output_test_list$fitGsigma3) + Output_test_list$media3
Output_test_list$peakG3 <- ifelse(Output_test_list$MaxGauss3 > Output_test_list$media3 + ispeak * Output_test_list$sd3, TRUE, FALSE)
Output_test_list$MaxGauss4 <- Output_test_list$fitGscale4 * dnorm(Output_test_list$fitGmu4, Output_test_list$fitGmu4, Output_test_list$fitGsigma4) + Output_test_list$media4
Output_test_list$peakG4 <- ifelse(Output_test_list$MaxGauss4 > Output_test_list$media4 + ispeak * Output_test_list$sd4, TRUE, FALSE)
Output_test_list$MaxGauss5 <- Output_test_list$fitGscale5 * dnorm(Output_test_list$fitGmu5, Output_test_list$fitGmu5, Output_test_list$fitGsigma5) + Output_test_list$media5
Output_test_list$peakG5 <- ifelse(Output_test_list$MaxGauss5 > Output_test_list$media5 + ispeak * Output_test_list$sd5, TRUE, FALSE)

Output_test_list <- add_column(Output_test_list, PredPicoM1 = NA, .after = "Maximo1")
Output_test_list <- add_column(Output_test_list, PredPicoM2 = NA, .after = "Maximo2")
Output_test_list <- add_column(Output_test_list, PredPicoM3 = NA, .after = "Maximo3")
Output_test_list <- add_column(Output_test_list, PredPicoM4 = NA, .after = "Maximo4")
Output_test_list <- add_column(Output_test_list, PredPicoM5 = NA, .after = "Maximo5")

Output_test_list <- add_column(Output_test_list, PredPicoG1 = NA, .after = "MaxGauss1")
Output_test_list <- add_column(Output_test_list, PredPicoG2 = NA, .after = "MaxGauss2")
Output_test_list <- add_column(Output_test_list, PredPicoG3 = NA, .after = "MaxGauss3")
Output_test_list <- add_column(Output_test_list, PredPicoG4 = NA, .after = "MaxGauss4")
Output_test_list <- add_column(Output_test_list, PredPicoG5 = NA, .after = "MaxGauss5")

Output_test_list$PredPicoM1 <- (Output_test_list$Maximo1 - Output_test_list$media1) / Output_test_list$sd1
Output_test_list$PredPicoM2 <- (Output_test_list$Maximo2 - Output_test_list$media2) / Output_test_list$sd2
Output_test_list$PredPicoM3 <- (Output_test_list$Maximo3 - Output_test_list$media3) / Output_test_list$sd3
Output_test_list$PredPicoM4 <- (Output_test_list$Maximo4 - Output_test_list$media4) / Output_test_list$sd4
Output_test_list$PredPicoM5 <- (Output_test_list$Maximo5 - Output_test_list$media5) / Output_test_list$sd5

Output_test_list$PredPicoG1 <- (Output_test_list$MaxGauss1 - Output_test_list$media1) / Output_test_list$sd1
Output_test_list$PredPicoG2 <- (Output_test_list$MaxGauss2 - Output_test_list$media2) / Output_test_list$sd2
Output_test_list$PredPicoG3 <- (Output_test_list$MaxGauss3 - Output_test_list$media3) / Output_test_list$sd3
Output_test_list$PredPicoG4 <- (Output_test_list$MaxGauss4 - Output_test_list$media4) / Output_test_list$sd4
Output_test_list$PredPicoG5 <- (Output_test_list$MaxGauss5 - Output_test_list$media5) / Output_test_list$sd5


Output_test_predM_df <- Output_test_list %>% select(c('filename', 'set',
                                                        'PredPicoM1', 'pico1', 
                                                        'PredPicoM2', 'pico2',
                                                        'PredPicoM3', 'pico3',
                                                        'PredPicoM4', 'pico4',
                                                        'PredPicoM5', 'pico5'))
Output_test_predG_df <- Output_test_list %>% select(c('filename', 'set',
                                                        'PredPicoG1', 'peakG1', 
                                                        'PredPicoG2', 'peakG2',
                                                        'PredPicoG3', 'peakG3',
                                                        'PredPicoG4', 'peakG4',
                                                        'PredPicoG5', 'peakG5'))

# done constructing the test_set as we've done with train_set


# -------------------- Evaluation on test set ----------------------------
#lets evaluate sensitivity and specificity at confusion matrix from this best_ispeak found at test_set
y_hat <- ifelse(Output_test_predM_df$PredPicoM1 < best_ispeak, "FALSE", "TRUE")
y_hat <- as.factor(y_hat)
confusionMatrix(data = y_hat, reference = as.factor(test_set$zona.1))

# Here we're still seeing an unbalance between Sensitivity and specificity, being the later one much too low
# Will try to weight differently how much important is to detect TRUE a peak than detect FALSE a no peak
# incorporate "beta" into the calculation of F_1

# ------------------ Introducing weight into F_1 score ------------------
beta_prop <- 0.25
rm(ispeak)
rm(F_1)
ispeak <- seq(from = 1.75, to = 6, by = 0.1)
F_1 <- map_dbl(ispeak, function(x){
  y_hat <- ifelse(Output_train_predM_df$PredPicoM1 < x, "FALSE", "TRUE")
  y_hat <- as.factor(y_hat)
  F_meas(data = y_hat, reference = as.factor(train_set$zona.1), beta = beta_prop)
})

ggplot() + aes(ispeak, F_1) + geom_point() + geom_line()
max(F_1)
best_ispeak <- ispeak[which.max(F_1)]
best_ispeak

y_hat <- ifelse(Output_test_predM_df$PredPicoM1 < best_ispeak, "FALSE", "TRUE")
y_hat <- as.factor(y_hat)
confusionMatrix(data = y_hat, reference = as.factor(test_set$zona.1))

# ----------- Comparing predicted and true results in test set -----------------
test_set <- add_column(test_set, PredPicoM1 = y_hat, .after = "zona.1")

index <- test_set %>% with(which(zona.1 == "TRUE" & PredPicoM1 == "TRUE")) # list of indexes of where TRUE peaks where found TRUE
lapply(index, function(x) plot(zona1[[2*x-1]], zona1[[2*x]], main = test_set$Archivo[x]))

# ----------- evaluating F_1 for every ispeak and every beta -------------------
beta_prop_vec <- seq(from = 0.2, to = 1, by = 0.1)
rm(F_1)

# F_1_2 function to calculate F_1 score as a function of ispeak and beta
F_1_2 <- function(x,y) {
  y_hat1 <- ifelse(Output_train_predM_df$PredPicoM1 < x, "FALSE", "TRUE")
  y_hat1 <- as.factor(y_hat1)
  F_meas(data = y_hat1, reference = as.factor(train_set$zona.1), beta = y)
}

F_1 <- sapply(beta_prop_vec, function(y) mapply(F_1_2,ispeak,y))
maximos <- apply(F_1, 2, which.max)    # indexes of where the F_1 maximum for each beta are
F_1 <- as.data.frame(F_1)
names(F_1) <- c(beta_prop_vec)
F_1 <- F_1 %>% add_column(valor = ispeak, .before = 1)
ispeak_max <- F_1$valor[maximos]     # ispeak of the F_1 maximum for each beta

results <- as.data.frame(beta_prop_vec)
results <- results %>% mutate(threshold1 = ispeak_max)
results <- results %>% mutate(maximum1 = apply(F_1[-1], 2, max))

F_1 <- melt(F_1 ,  id.vars = 'valor', variable.name = 'beta')
ggplot(F_1, aes(valor, value)) + geom_line(aes(colour = beta)) + geom_point(aes(colour = beta))

# -------------- Comparing method with just guessing ---------------------------
# ROC 
probs <- seq(0, 1, length.out = 11)
guessing <- map_df(probs, function(p) {
  y_hat_prob <- sample(c("FALSE", "TRUE"), length(test_index), 
                       replace = TRUE, prob = c(p, 1-p)) %>%
    factor(levels = c("FALSE", "TRUE"))
  list(method = "Guess", 
       recall = sensitivity(y_hat_prob, as.factor(test_set$zona.1)),
       precision = precision(y_hat_prob, as.factor(test_set$zona.1)))
})
ggplot() + aes(guessing$recall, guessing$precision) + geom_point() + geom_line() + geom_label(aes(label = probs))

max_method <- map_df(ispeak, function(peak) {
  y_hat_maxROC <- ifelse(Output_train_predM_df$PredPicoM1 < peak, "FALSE", "TRUE") %>%
    factor(levels = c("FALSE", "TRUE"))
  list(method = "Maximum", 
       recall = sensitivity(y_hat_maxROC, as.factor(test_set$zona.1)), 
       precision = precision(y_hat_maxROC, as.factor(test_set$zona.1)))
})
ggplot() + aes(max_method$recall, max_method$precision) + geom_point() + geom_line() + geom_label(aes(label = ispeak))

guessing <- mutate(guessing, label = probs)
max_method <- mutate(max_method, label = ispeak)
methods_comp_plot <- rbind(guessing, max_method)
ggplot(methods_comp_plot, aes(x=recall, y=precision, group=method, col=method, fill=method)) + geom_point() + geom_line() + 
  geom_text(aes(recall, precision, label = methods_comp_plot$label), nudge_x = 0.01, nudge_y = 0.005)

# ------------- Applying ML code to data ---------------------------------------
# first lets know where data actually are
fileEvalpath <- ifelse(osystem == "ubuntu", 
                       "/media/federico/WD300_2/Fisica/Otros/Cesar/Evaluation/",
                       "f:/Fisica/Otros/Cesar/Evaluation/")

fileEvalfiles <- read.table(file = paste(fileEvalpath, "folders.txt", sep = ""))
fileEvalfiles <- rename(fileEvalfiles, fullfolder = V1)
fileEvalfiles <- fileEvalfiles %>% mutate(namefolder = NA)
fileEvalfiles$namefolder <- substr(fileEvalfiles$fullfolder, 
                                   start = ifelse(osystem == "ubuntu", 55, 34),
                                   stop = str_length(fileEvalfiles$fullfolder))
fileEvalfiles$fullfolder <- paste(fileEvalfiles$fullfolder, "/", sep="")

filenames_full <- fileEvalfiles$fullfolder %>% 
  lapply(function(x) list.files(x, pattern = "*.txt", full.names = TRUE))  #list with all filenames and paths

lst_Eval_resultsM <- vector(mode = 'list', length = length(filenames_full)) #list with all results (empty for now)

# Now lets calculate the mean and sd on both edges of zones 1 through 5 for every dataset
#--example with only 1 folder


folder_number <- seq(1:length(filenames_full))  # tengo que secuenciarlo con un seq(1:length(filenames_full))

for (ii in folder_number) {
  print(sprintf("ciclo %i", ii))
  dataEvallist <- lapply(filenames_full[[ii]], function(x) read.table(x))
  seq_long_folder <- seq(1:length(dataEvallist))

  zona1 <- lapply(seq_long_folder, 
                  function(x){
                    dataEvallist[[x]] %>% filter(V1 > (validos_filt$twothetaM[1] - validos_filt$min[1]/2) & 
                                                 V1 < (validos_filt$twothetaM[1] + validos_filt$min[1]/2))
                  })
  subzona1 <- lapply(seq_long_folder, 
                     function(x){
                       zona1[[x]] %>% filter(V1 < (validos_filt$twothetaM[1] - validos_filt$min[1]/4) | 
                                             V1 > (validos_filt$twothetaM[1] + validos_filt$min[1]/4)) %>% 
                       select(V2)
                })
  max_M1 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) zona1[[x]]$V2 %>% max()))))
  sd1 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona1[[x]]$V2 %>% sd()))))
  media1 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona1[[x]]$V2 %>% mean()))))

  print(sprintf("Done zone 1 folder %i", ii))
    
  zona2 <- lapply(seq_long_folder, 
                  function(x){
                    dataEvallist[[x]] %>% filter(V1 > (validos_filt$twothetaM[2] - validos_filt$min[2]/2) &
                                                 V1 < validos_filt$twothetaM[2] + validos_filt$min[2]/2)
                    
                  })
  
  subzona2 <- lapply(seq_long_folder,
                    function(x){
                      zona2[[x]] %>% filter(V1 < (validos_filt$twothetaM[2] - validos_filt$min[2]/4) |
                                            V1 > (validos_filt$twothetaM[2] + validos_filt$min[2]/4)) %>%
                      select(V2)
                    })
  max_M2 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) zona2[[x]]$V2 %>% max()))))
  sd2 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona2[[x]]$V2 %>% sd()))))
  media2 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona2[[x]]$V2 %>% mean()))))
  
  print(sprintf("Done zone 2 folder %i", ii))
  
  zona3 <- lapply(seq_long_folder,
                  function(x){
                    dataEvallist[[x]] %>% filter(V1 > (validos_filt$twothetaM[3] - validos_filt$min[3]/2) &
                                                 V1 < (validos_filt$twothetaM[3] + validos_filt$min[3]/2))
                  })
  subzona3 <- lapply(seq_long_folder, 
                     function(x){
                       zona3[[x]] %>% filter(V1 < (validos_filt$twothetaM[3] - validos_filt$min[3]/4) |
                                             V1 > (validos_filt$twothetaM[3] + validos_filt$min[3]/4)) %>%
                       select(V2)
                     })
  max_M3 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) zona3[[x]]$V2 %>% max()))))
  sd3 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona3[[x]]$V2 %>% sd()))))
  media3 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona3[[x]]$V2 %>% mean()))))
  
  print(sprintf("Done zone 3 folder %i", ii))
  
  zona4 <- lapply(seq_long_folder, 
                  function(x){
                    dataEvallist[[x]] %>% filter(V1 > (validos_filt$twothetaM[4] - validos_filt$min[4]/2) & 
                                                   V1 < (validos_filt$twothetaM[4] + validos_filt$min[4]/2))
                  })
  subzona4 <- lapply(seq_long_folder,
                     function(x){
                       zona4[[x]] %>% filter(V1 < (validos_filt$twothetaM[4] - validos_filt$min[4]/4) |
                                             V1 > (validos_filt$twothetaM[4] + validos_filt$min[4]/4)) %>%
                       select(V2)
                     })
  max_M4 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) zona4[[x]]$V2 %>% max()))))
  sd4 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona4[[x]]$V2 %>% sd()))))
  media4 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona4[[x]]$V2 %>% mean()))))
  
  print(sprintf("Done zone 4 folder %i", ii))
  
  zona5 <- lapply(seq_long_folder, 
                  function(x){
                    dataEvallist[[x]] %>% filter(V1 > (validos_filt$twothetaM[5] - validos_filt$min[5]/2) &
                                                 V1 < (validos_filt$twothetaM[5] + validos_filt$min[5]/2))
                  })
  subzona5 <- lapply(seq_long_folder,
                     function(x){
                       zona5[[x]] %>% filter(V1 < (validos_filt$twothetaM[5] - validos_filt$min[5]/4) | 
                                             V1 > (validos_filt$twothetaM[5] + validos_filt$min[5]/4)) %>%
                       select(V2)
                     })
  max_M5 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) zona5[[x]]$V2 %>% max()))))
  sd5 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona5[[x]]$V2 %>% sd()))))
  media5 <- as.vector(t(as.data.frame(lapply(seq_long_folder, function(x) subzona5[[x]]$V2 %>% mean()))))
  
  print(sprintf("Done zone 5 folder %i", ii))
  
  # now lets build the output dataframe with the results
  df_Eval_resultsM <- data.frame("full_name" = filenames_full[[ii]], 
                                 "zona1_max" = max_M1, "zona1_mean" = media1, "zona1_sd" = sd1, prominence1 = NA, peak1 = NA,
                                 "zona2_max" = max_M2, "zona2_mean" = media2, "zona2_sd" = sd2, prominence2 = NA, peak2 = NA,  
                                 "zona3_max" = max_M3, "zona3_mean" = media3, "zona3_sd" = sd3, prominence3 = NA, peak3 = NA,
                                 "zona4_max" = max_M4, "zona4_mean" = media4, "zona4_sd" = sd4, prominence4 = NA, peak4 = NA,
                                 "zona5_max" = max_M5, "zona5_mean" = media5, "zona5_sd" = sd5, prominence5 = NA, peak5 = NA) 
  
  df_Eval_resultsM$full_name <- substring(df_Eval_resultsM$full_name, first = ifelse(osystem == "ubuntu", 55, 34))
  
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(prominence1 = (zona1_max - zona1_mean) / zona1_sd)
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(prominence2 = (zona2_max - zona2_mean) / zona2_sd)
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(prominence3 = (zona3_max - zona3_mean) / zona3_sd)
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(prominence4 = (zona4_max - zona4_mean) / zona4_sd)
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(prominence5 = (zona5_max - zona5_mean) / zona5_sd)
  
  peak_threshold <-  results$threshold1[which.max(results$maximum1)]
  
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(peak1 = ifelse(prominence1 > peak_threshold, "TRUE", "FALSE"))
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(peak2 = ifelse(prominence2 > peak_threshold, "TRUE", "FALSE"))
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(peak3 = ifelse(prominence3 > peak_threshold, "TRUE", "FALSE"))
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(peak4 = ifelse(prominence4 > peak_threshold, "TRUE", "FALSE"))
  df_Eval_resultsM <- df_Eval_resultsM %>% mutate(peak5 = ifelse(prominence5 > peak_threshold, "TRUE", "FALSE"))
  
  
  lst_Eval_resultsM[[ii]] <- df_Eval_resultsM

  print(sprintf("Done building df result folder %i", ii))
  
  write.csv(df_Eval_resultsM, paste(fileEvalfiles$fullfolder[ii], fileEvalfiles$namefolder[ii],".csv", sep=""), row.names=TRUE)

  print(sprintf("Done writing output file folder %i", ii))
}




