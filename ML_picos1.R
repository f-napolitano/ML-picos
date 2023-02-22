library(tidyverse)
library(ggplot2)
library(data.table)
library(pracma)

#---------- Declaracion de parametros ---------------------------
ispeak <- 10 # cuantas veces la SD tiene que tener el maximo de un pico para ser considerado como tal

#--------- Importacion de archivos y preparacion de dataframes
pathfiles <- "f:/Fisica/ML/Prueba_picos/"
filelist <- paste(pathfiles, list.files(pathfiles), sep="")
datalist <- lapply(filelist, function(x)read.table(x))
datafr <- bind_cols(datalist)
Temp1 <- datafr[1, c(FALSE, TRUE)] %>% transpose() %>% mutate(nombre = list.files(pathfiles))

datafr <- datafr[-1,]

diff1 = data.frame(lapply(datafr, function (x) as.double(as.character(x))))

diff1 <- diff1 %>% filter(V1...1 > 29.75 & V1...1 < 31.25)

# ------------- calculo de sd y mean ---------------------------
intsd <- diff1 %>% filter(V1...1 < 30.3 | V1...1 > 30.8) 
intsd <- intsd[, c(FALSE, TRUE)] 

Temp1 <- Temp1 %>% mutate(sd = transpose(data.frame(lapply(intsd, function (x) sd(x)))))

Temp1 <- Temp1 %>% mutate(media = transpose(data.frame(lapply(intsd, function (x) mean(x)))))

#--------------- determinacion de picos ------------------------
intsd <- diff1[, c(FALSE, TRUE)]

Temp1 <- Temp1 %>% mutate(Maximo = sapply(intsd, function (x) max(x)))
Temp1 <- Temp1 %>% mutate(pico = diag( sapply(intsd, function(x) { ifelse(max(x) > Temp1$media + ispeak * Temp1$sd, TRUE, FALSE)} )) )

fitG =
  function(x,y,mu,sig,scale,background){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2]) + background
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }


Temp1$fitGmu <- ""
Temp1$fitGsigma <- ""
Temp1$fitGscale <- ""

for(i in 1:round(ncol(diff1)/2)) {
  
  fitGauss <- as.list(fitG(diff1[[2*i-1]], diff1[[2*i]], 30.55, 1, 100, Temp1$media[[1]][i]))
  peakpar <- as.list(fitGauss$par)
  Temp1$fitGmu[i] <- as.numeric(peakpar[1])
  Temp1$fitGsigma[i] <- as.numeric(peakpar[2])
  Temp1$fitGscale[i] <- as.numeric(peakpar[3])
  
  plot(diff1[[2*i-1]], diff1[[2*i]])
  lines(diff1[[2*i-1]], as.numeric(peakpar[3]) * dnorm(diff1[[2*i-1]], as.numeric(peakpar[1]), as.numeric(peakpar[2])) + Temp1$media[[1]][i])
  
}
Temp1$fitGmu <- as.numeric(Temp1$fitGmu)
Temp1$fitGscale <- as.numeric(Temp1$fitGscale)
Temp1$fitGsigma <- as.numeric(Temp1$fitGsigma)

Temp1 <- Temp1 %>% mutate(MaxGauss = fitGscale * dnorm(fitGmu, fitGmu, fitGsigma) + media)
Temp1 <- Temp1 %>% mutate(picoG = ifelse(MaxGauss > media + ispeak * sd, TRUE, FALSE))

grwgwg