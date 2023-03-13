library(tidyverse)
library(ggplot2)
library(data.table)
library(pracma)
library(dplyr)
library(caret)

osystem <- "windows"

#----------- Determinacion reflexion como predictor -------------
filepath <- ifelse(osystem == "ubuntu", "/media/federico/WD300_2/Fisica/Otros/Cesar/SimulacionXRD/" ,"f:/Fisica/Otros/Cesar/SimulacionXRD/")
marten <- read.csv(file = paste(filepath, "Martensita.csv", sep = ""), header = TRUE) %>% select(2:4)
aust <- read.csv(file = paste(filepath, "Austenita.csv", sep = ""), header = TRUE) %>% select(2:4)
gammaphase <- read.csv(file = paste(filepath, "Gamma.csv", sep = ""), header = TRUE) %>% select(2:4)

validos <- marten

diff_aust_validos <- sapply(validos$twothetaM, function(x) aust$twothetaA[ which.min(abs(x - aust$twothetaA))] - x)
validos <- validos %>% mutate(diff_aust = diff_aust_validos)
rm(diff_aust_validos)
diff_gamma_validos <- sapply(validos$twothetaM, function(x) gammaphase$twothetaG[ which.min((abs(x - gammaphase$twothetaG)))] - x)
validos <- validos %>% mutate(diff_gamma = diff_gamma_validos)
rm(diff_gamma_validos)

validos <- validos %>% rowwise() %>% mutate(min = min(abs(diff_aust), abs(diff_gamma))) %>% 
  filter(abs(diff_aust) > 0.1 & abs(diff_gamma) > 0.1) %>% 
  mutate(distancia = sqrt(diff_aust**2 + diff_gamma**2))

#validos <- validos %>% filter(abs(diff_aust) > 0.1 & abs(diff_gamma) > 0.1)
#validos <- validos %>% mutate(distancia = sqrt(diff_aust**2 + diff_gamma**2))

validos %>% ggplot() + aes(abs(diff_aust), abs(diff_gamma), label = round(twothetaM, 2)) + geom_point() + geom_text(hjust=0, vjust=-0.5)
validos %>% ggplot() + aes(twothetaM, distancia) + geom_point()

distancia_min <- 0.5
validos_filt <- validos %>% filter(distancia > distancia_min & abs(diff_gamma) > 0.25 & abs(diff_aust) > 0.25)

#----------- constatacion manual de posicion de picos y redefinicion --------
validos_filt[3,1] <- 8.70


#---------- Declaracion de parametros ---------------------------
ispeak <- 4 # cuantas veces la SD tiene que tener el maximo de un pico para ser considerado como tal

#--------- Importacion de archivos y preparacion de dataframes
pathfiles <- ifelse(osystem == "ubuntu", "/media/federico/WD300_2/Fisica/Otros/Cesar/Exploracion/", "f:/Fisica/Otros/Cesar/Exploracion/")
filelist <- paste(pathfiles, list.files(pathfiles), sep="")
datalist <- lapply(filelist, function(x)read.table(x))
datafr <- bind_cols(datalist)
# Temp1 <- datafr[1, c(FALSE, TRUE)] %>% transpose() %>% mutate(nombre = list.files(pathfiles)) esto es en caso de tener una fila encabezado con parametro
# datafr <- datafr[-1,] idem anterior

Temp1 <- list.files(pathfiles)

#diff1 = data.frame(lapply(datafr, function (x) as.double(as.character(x))))
diff1 <- datafr

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


Temp1 <- as.data.frame(Temp1) %>% mutate(sd1 = transpose(data.frame(lapply(intsd1, function (x) sd(x)))), 
                                         media1 = transpose(data.frame(lapply(intsd1, function (x) mean(x)))), 
                                         sd2 = transpose(data.frame(lapply(intsd2, function (x) sd(x)))), 
                                         media2 = transpose(data.frame(lapply(intsd2, function (x) mean(x)))), 
                                         sd3 = transpose(data.frame(lapply(intsd3, function (x) sd(x)))), 
                                         media3 = transpose(data.frame(lapply(intsd3, function (x) mean(x)))), 
                                         sd4 = transpose(data.frame(lapply(intsd4, function (x) sd(x)))), 
                                         media4 = transpose(data.frame(lapply(intsd4, function (x) mean(x)))), 
                                         media5 = transpose(data.frame(lapply(intsd5, function (x) mean(x)))), 
                                         sd5 = transpose(data.frame(lapply(intsd5, function (x) sd(x))))
                                         )

#--------------- determinacion de picos ------------------------
intsd1 <- zona1[, c(FALSE, TRUE)]
intsd2 <- zona2[, c(FALSE, TRUE)]
intsd3 <- zona3[, c(FALSE, TRUE)]
intsd4 <- zona4[, c(FALSE, TRUE)]
intsd5 <- zona5[, c(FALSE, TRUE)]

Temp1 <- Temp1 %>% mutate(Maximo1 = sapply(intsd1, function (x) max(x)), 
                          pico1 = diag( sapply(intsd1, function(x) { ifelse(max(x) > media1 + ispeak * sd1, TRUE, FALSE)} )), 
                          Maximo2 = sapply(intsd2, function(x) max(x)),
                          pico2 = diag( sapply(intsd2, function(x) { ifelse(max(x) > media2 + ispeak * sd2, TRUE, FALSE)} )),
                          Maximo3 = sapply(intsd3, function(x) max(x)),
                          pico3 = diag( sapply(intsd3, function(x) { ifelse(max(x) > media3 + ispeak * sd3, TRUE, FALSE)} )),
                          Maximo4 = sapply(intsd4, function(x) max(x)),
                          pico4 = diag( sapply(intsd4, function(x) { ifelse(max(x) > media4 + ispeak * sd4, TRUE, FALSE)} )),
                          Maximo5 = sapply(intsd5, function(x) max(x)),
                          pico5 = diag( sapply(intsd5, function(x) { ifelse(max(x) > media5 + ispeak * sd5, TRUE, FALSE)} )),
                          )


fitG =
  function(x,y,mu,sig,scale,background){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2]) + background
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }

Temp1[, c('fitGmu1', 'fitGsigma1', 'fitGscale1', 'MaxGauss1', 'peakG1',
          'fitGmu2', 'fitGsigma2', 'fitGscale2', 'MaxGauss2', 'peakG2',
          'fitGmu3', 'fitGsigma3', 'fitGscale3', 'MaxGauss3', 'peakG3',
          'fitGmu4', 'fitGsigma4', 'fitGscale4', 'MaxGauss4', 'peakG4',
          'fitGmu5', 'fitGsigma5', 'fitGscale5', 'MaxGauss5', 'peakG5')
      ]<- NA


for(i in 1:round(ncol(diff1)/2)) {
  
  fitGauss <- as.list(fitG(zona1[[2*i-1]], zona1[[2*i]], validos_filt$twothetaM[1], 1, 1, Temp1$media1[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Temp1$fitGmu1[i] <- as.numeric(peakpar[1])
  Temp1$fitGsigma1[i] <- as.numeric(peakpar[2])
  Temp1$fitGscale1[i] <- as.numeric(peakpar[3])

  plot(zona1[[2*i-1]], zona1[[2*i]], main = Temp1$Temp1[i])
  lines(zona1[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona1[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Temp1$media1[[1]][i], col = "red", pch = 20)
  
  
  fitGauss <- as.list(fitG(zona2[[2*i-1]], zona2[[2*i]], validos_filt$twothetaM[2], 1, 1, Temp1$media2[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Temp1$fitGmu2[i] <- as.numeric(peakpar[1])
  Temp1$fitGsigma2[i] <- as.numeric(peakpar[2])
  Temp1$fitGscale2[i] <- as.numeric(peakpar[3])
  
  plot(zona2[[2*i-1]], zona2[[2*i]], main = Temp1$Temp1[i])
  lines(zona2[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona2[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Temp1$media2[[1]][i], col = "red", pch = 20)

  fitGauss <- as.list(fitG(zona3[[2*i-1]], zona3[[2*i]], validos_filt$twothetaM[3], 1, 1, Temp1$media3[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Temp1$fitGmu3[i] <- as.numeric(peakpar[1])
  Temp1$fitGsigma3[i] <- as.numeric(peakpar[2])
  Temp1$fitGscale3[i] <- as.numeric(peakpar[3])
  
  plot(zona3[[2*i-1]], zona3[[2*i]], main = Temp1$Temp1[i])
  lines(zona3[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona3[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Temp1$media3[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona4[[2*i-1]], zona4[[2*i]], validos_filt$twothetaM[4], 1, 1, Temp1$media4[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Temp1$fitGmu4[i] <- as.numeric(peakpar[1])
  Temp1$fitGsigma4[i] <- as.numeric(peakpar[2])
  Temp1$fitGscale4[i] <- as.numeric(peakpar[3])
  
  plot(zona4[[2*i-1]], zona4[[2*i]], main = Temp1$Temp1[i])
  lines(zona4[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona4[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Temp1$media4[[1]][i], col = "red", pch = 20)
  
  fitGauss <- as.list(fitG(zona5[[2*i-1]], zona5[[2*i]], validos_filt$twothetaM[5], 1, 1, Temp1$media5[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Temp1$fitGmu5[i] <- as.numeric(peakpar[1])
  Temp1$fitGsigma5[i] <- as.numeric(peakpar[2])
  Temp1$fitGscale5[i] <- as.numeric(peakpar[3])
  
  plot(zona5[[2*i-1]], zona5[[2*i]], main = Temp1$Temp1[i])
  lines(zona5[[2*i-1]], as.numeric(peakpar[3]) * dnorm(zona5[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Temp1$media5[[1]][i], col = "red", pch = 20)
  
  
  
  }


Temp1$MaxGauss1 <- Temp1$fitGscale1 * dnorm(Temp1$fitGmu1, Temp1$fitGmu1, Temp1$fitGsigma1) + Temp1$media1
Temp1$peakG1 <- ifelse(Temp1$MaxGauss1 > Temp1$media1 + ispeak * Temp1$sd1, TRUE, FALSE)
Temp1$MaxGauss2 <- Temp1$fitGscale2 * dnorm(Temp1$fitGmu2, Temp1$fitGmu2, Temp1$fitGsigma2) + Temp1$media2
Temp1$peakG2 <- ifelse(Temp1$MaxGauss2 > Temp1$media2 + ispeak * Temp1$sd2, TRUE, FALSE)
Temp1$MaxGauss3 <- Temp1$fitGscale3 * dnorm(Temp1$fitGmu3, Temp1$fitGmu3, Temp1$fitGsigma3) + Temp1$media3
Temp1$peakG3 <- ifelse(Temp1$MaxGauss3 > Temp1$media3 + ispeak * Temp1$sd3, TRUE, FALSE)
Temp1$MaxGauss4 <- Temp1$fitGscale4 * dnorm(Temp1$fitGmu4, Temp1$fitGmu4, Temp1$fitGsigma4) + Temp1$media4
Temp1$peakG4 <- ifelse(Temp1$MaxGauss4 > Temp1$media4 + ispeak * Temp1$sd4, TRUE, FALSE)
Temp1$MaxGauss5 <- Temp1$fitGscale5 * dnorm(Temp1$fitGmu5, Temp1$fitGmu5, Temp1$fitGsigma5) + Temp1$media5
Temp1$peakG5 <- ifelse(Temp1$MaxGauss5 > Temp1$media5 + ispeak * Temp1$sd5, TRUE, FALSE)

# -------------- creating train and test set --------------

y <- Temp1$Maximo1

test_index <- createDataPartition(y, times = 1, p = 0.5, list = TRUE)



# ------------ filtro FIR --------------
h <- fir1(10, 0.1, "low")
h1 <- fir1(5,0.1,"low")
f1 <- filter(h1, zona1[[36]])
f2 <- filter(h, zona1[[36]])
plot(zona1[[35]], zona1[[36]])
lines(zona1[[35]], f1, col = "red", lwd = 2)
lines(zona1[[35]], f2, col = "blue", lwd = 2)