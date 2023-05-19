setwd("C:/Users/Faye/桌面/assign")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Advanced Life insurance mathematics: #
#        Assignment - Jing Yang        #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

###### 0. Settings ######
start_year <- 1876
KULbg <- "#116E8A"

library(tidyverse)
library(demography)
library(forecast)
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(astsa)
library(MultiMoMo)
library(abind)

females<- read.table("females.csv", sep=",",head=TRUE)
males<- read.table("males.csv", sep=",",head=TRUE)
uni<- read.table("uni.csv", sep=",",head=TRUE)

names(females)# "Year" "Age"  "mx"   "qx"   "ax"   "lx"   "dx"   "Lx"   "Tx"   "ex"  

str(females)

min(females$Year)# 1876
max(females$Year)# 2021


min(females$Age)# 0
max(females$Age)# 110

###### 1. check the data in 2019 per gender ######

fem_2019<- filter(females, Year == 2019)
male_2019<- filter(males, Year == 2019)

#number of death dx----------------------
par(mfcol=c(2,2))

g_fem_dx <- ggplot(fem_2019, aes(Age, dx)) + 
  geom_step(col = KULbg) + 
  theme_bw() +
  ggtitle("Switzerland - females, 2019") 
g_fem_dx

g_male_dx <- ggplot(male_2019, aes(Age, dx)) + 
  geom_step(col = "red") + 
  theme_bw() +
  ggtitle("Switzerland - males, 2019")
g_male_dx

#combined
g_dx <- g_male_dx +
  geom_step(data = fem_2019, aes(Age, dx), col = KULbg) + 
  theme_bw() +
  ggtitle("Switzerland, 2019") +
  labs(y = "dx")
g_dx


#life expectancy on 2019 for all ages ex----------------------

g_fem_ex <- ggplot(fem_2019, aes(Age, ex)) + 
  geom_line(col = KULbg, size = 1) + 
  theme_bw() +
  ggtitle("Switzerland - females, 2019") +
  labs(y = "ex")
g_fem_ex

g_male_ex <- ggplot(male_2019, aes(Age, ex)) + 
  geom_line(col = "red", size = 1) + 
  theme_bw() +
  ggtitle("Switzerland - males, 2019") +
  labs(y = "ex")
g_male_ex

#combined
g_ex <- g_male_ex +
  geom_line(data = fem_2019, aes(Age, ex), col = KULbg, size = 1) + 
  theme_bw() +
  ggtitle("Switzerland, 2019") +
  labs(y = "ex")
g_ex

#log(qx)---------------------

g_fem_qx <- ggplot(fem_2019, aes(Age, log(qx))) + 
  geom_point(col = KULbg) + 
  geom_line(col = KULbg) + 
  theme_bw() +
  ggtitle("Switzerland - females, 2019") + 
  labs(y = bquote(ln(q[x])))
g_fem_qx

g_male_qx <- ggplot(male_2019, aes(Age, log(qx))) + 
  geom_point(col = "red") + 
  geom_line(col = "red") + 
  theme_bw() +
  ggtitle("Switzerland - males, 2019") + 
  labs(y = bquote(ln(q[x])))
g_male_qx

#combined
g_qx <- g_male_qx + 
  geom_point(data = fem_2019, aes(Age, log(qx)), col = KULbg) + 
  geom_line(data = fem_2019, aes(Age, log(qx)), col = KULbg)
theme_bw() +
  ggtitle("Switzerland, 2019") + 
  labs(y = bquote(ln(q[x]))) 
g_qx



#using the earliest year as start year to initial the process
start_year <- 1971
end_year <- 2019
years <- start_year:end_year
ages <- 0:90

females_est <- filter(females, Year <= end_year & Year >= start_year)
females_est_sub  <- subset(females_est, Age <= 90)


names(females_est_sub)
str(females_est_sub)

#### 
dtx <- matrix(females_est_sub$dx, nrow = length(years), byrow = TRUE)
for(i in 1:length(start_year:end_year)){
  for(j in 1:length(0:90))
    if(dtx[i,j] == 0){
      dtx[i,j] <- dtx[i,j-1]
    }
}
dim(dtx)

etx <- dtx / matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
etx[is.infinite(etx)] <- 0
etx[is.na(etx)] <- 0
for(i in 1:length(start_year:end_year)){
  for(j in 1:length(0:90))
  if(etx[i,j] == 0){
    etx[i,j] <- etx[i,j-1]
  }
}
dim(etx)

mtx <- matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
for(i in 1:length(start_year:end_year)){
  for(j in 1:length(0:90))
    if(mtx[i,j] == 0){
      mtx[i,j] <- mtx[i,j-1]
    }
}
dim(mtx)

df <- expand.grid(Year = years, Age = ages)
df$rates <- as.vector(mtx)
df <- filter(df, df$rates > 0)
df$rates <- log(df$rates)
names(df) <- c("Year", "Age", "logRate")

p_mx_1876 <- ggplot(df, aes(x = Age, y = logRate, group = Year)) + 
  geom_line(aes(colour = Year), linewidth = 1, linetype = 1) +
  scale_colour_gradientn(colours = rainbow(10)) +
  scale_x_continuous(breaks = seq(ages[1], tail(ages, 1) + 1, 10)) +
  theme_bw() + ylab(expression("log" ~ m[x])) + xlab("Age (x)") +
  ggtitle("Switzerland, 1876-2019")
p_mx_1876


##after the determination of the start year
p_mx_1971 <- ggplot(df, aes(x = Age, y = logRate, group = Year)) + 
  geom_line(aes(colour = Year), size = 1, linetype = 1) +
  scale_colour_gradientn(colours = rainbow(10)) +
  scale_x_continuous(breaks = seq(ages[1], tail(ages, 1) + 1, 10)) +
  theme_bw() + ylab(expression("log" ~ m[x])) + xlab("Age (x)") +
  ggtitle("Switzerland, 1971-2019") 
p_mx_1971


#### 2. Poisson likelihood and NR scheme ####


source('fitModels.R')

LCfit701 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
names(LCfit701)
sum(LCfit701$beta2)
sum(LCfit701$kappa2)


data_period <- tibble(year = years, fit = LCfit701$kappa2)
data_age <- tibble(age = ages, fit_alpha = LCfit701$beta1, fit_beta = LCfit701$beta2)

g_1 <- ggplot(data_age) + geom_point(aes(age, fit_alpha)) + 
  geom_line(aes(age, fit_alpha), col = "black") +
  theme_bw() +
  ggtitle("Switzerland - females, 1876 - 2019, Lee Carter, Poisson") + 
  labs(y = bquote(hat(beta)[x]^"(1)")) 

g_2 <- ggplot(data_age) + geom_point(aes(age, fit_beta)) + 
  geom_line(aes(age, fit_beta), col = "black") +
  theme_bw() + ggtitle("") +
  labs(y = bquote(hat(beta)[x]^"(2)")) 

g_3 <- ggplot(data_period) + geom_point(aes(year, fit)) + 
  geom_line(aes(year, fit), col = "black") +
  theme_bw() + ggtitle("") + 
  labs(y = bquote(hat(kappa)[t]^"(2)")) 

grid.arrange(g_1, g_2, g_3, ncol = 2)


## ------------------------------------------------------------------------------------

LCfit701_1876 <- LCfit701

#2.1 looking for the starting year by using chosen model: random walking and AR(1)

## 2.1.1 project with the model AR(1) for different starting year: check if it is stable  
list3 <- list()

for (k in 1876:2000) {
  start_year <- k
  years <- start_year:end_year
  ages <- 0:90
  females_est <- filter(females, Year <= end_year & Year >= k)
  females_est_sub  <- subset(females_est, Age <= 90)
  dtx <- matrix(females_est_sub$dx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(dtx[i,j] == 0){
        dtx[i,j] <- dtx[i,j-1]
      }
  }
  
  etx <- dtx / matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
  etx[is.infinite(etx)] <- 0
  etx[is.na(etx)] <- 0
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(etx[i,j] == 0){
        etx[i,j] <- etx[i,j-1]
      }
  }
  
  mtx <- matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(mtx[i,j] == 0){
        mtx[i,j] <- mtx[i,j-1]
      }
  }
  
  LCfit701 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
  
  time_series <- Arima(LCfit701$kappa2, 
                          order = c(1, 0, 0), 
                          include.drift = TRUE)
  
  list3 <- append(list3, time_series$coef)
}
ar1_coef_fem <- data.frame(matrix(list3, ncol = 3, byrow = TRUE))
names(ar1_coef_fem) <- c("ar1", "intercept", "drift")
rownames(ar1_coef_fem) <- c(1876:2000)

ar1_coef_fem[,1] <- as.numeric(ar1_coef_fem[,1])
ar1_coef_fem[,2] <- as.numeric(ar1_coef_fem[,2])
ar1_coef_fem[,3] <- as.numeric(ar1_coef_fem[,3])

max(abs(as.numeric(ar1_coef_fem$ar1)))
#[1] 0.9548458, all smaller than 1, have stationary, pointless to compare

#same as males, not coded yet

## 2.1.2 project with the model random walking with drift for different starting year: check if it is stable
#for females
list4 <- list()

for (k in 1876:2000) {
  start_year <- k
  years <- start_year:end_year
  ages <- 0:90
  females_est <- filter(females, Year <= end_year & Year >= k)
  females_est_sub  <- subset(females_est, Age <= 90)
  dtx <- matrix(females_est_sub$dx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(dtx[i,j] == 0){
        dtx[i,j] <- dtx[i,j-1]
      }
  }
  
  etx <- dtx / matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
  etx[is.infinite(etx)] <- 0
  etx[is.na(etx)] <- 0
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(etx[i,j] == 0){
        etx[i,j] <- etx[i,j-1]
      }
  }
  
  mtx <- matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(mtx[i,j] == 0){
        mtx[i,j] <- mtx[i,j-1]
      }
  }
  
  LCfit701 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
  
  time_series <- Arima(LCfit701$kappa2, 
                          order = c(0, 1, 0), 
                          include.drift = TRUE)
  
  list4 <- append(list4, time_series$aic)
}

rm_aic_fem <- as.numeric(list4)
rm_aic_fem <- diff(rm_aic_fem, lag = 1)
which(rm_aic_fem %in% head(sort(rm_aic_fem, decreasing = FALSE), n = 10))
#[1] 42 43 45 47 54 59 65 68 74 86
which(rm_aic_fem %in% head(sort(rm_aic_fem, decreasing = FALSE), n = 5))
#[1] 42 43 45 68 86 -- 1918, 1919, 1944, 1921, 1952


#for males
list5 <- list()

for (k in 1876:2000) {
  start_year <- k
  years <- start_year:end_year
  ages <- 0:90
  males_est <- filter(males, Year <= end_year & Year >= k)
  males_est_sub  <- subset(males_est, Age <= 90)
  dtx <- matrix(males_est_sub$dx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(dtx[i,j] == 0){
        dtx[i,j] <- dtx[i,j-1]
      }
  }
  
  etx <- dtx / matrix(males_est_sub$mx, nrow = length(years), byrow = TRUE)
  etx[is.infinite(etx)] <- 0
  etx[is.na(etx)] <- 0
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(etx[i,j] == 0){
        etx[i,j] <- etx[i,j-1]
      }
  }
  
  mtx <- matrix(males_est_sub$mx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(mtx[i,j] == 0){
        mtx[i,j] <- mtx[i,j-1]
      }
  }
  
  LCfit701 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
  
  time_series <- Arima(LCfit701$kappa2, 
                          order = c(0, 1, 0), 
                          include.drift = TRUE)
  
  list5 <- append(list5, time_series$aic)
}

rm_aic_male <- as.numeric(list5)
rm_aic_male <- diff(rm_aic_male, lag = 1)
which(rm_aic_male %in% head(sort(rm_aic_male, decreasing = FALSE), n = 10))
#[1] 42 43 45 53 59 65 68 75 86 88
which(rm_aic_male %in% head(sort(rm_aic_male, decreasing = FALSE), n = 5))
#[1] 42 43 45 68 86 -- 1918, 1919, 1944, 1921, 1952

#hence 1944 and 1952 are chosen


#2.2 looking for the ideal ARIMA　model for chosen starting years: 1944, 1952, 1971
#1971 is the year found during the random tests
#will compare p_value, AIC

#fit the model for different starting years
start_year <- 1971#change here
years <- start_year:end_year
ages <- 0:90

females_est <- filter(females, Year <= end_year & Year >= start_year)
females_est_sub  <- subset(females_est, Age <= 90)

dtx <- matrix(females_est_sub$dx, nrow = length(years), byrow = TRUE)
for(i in 1:length(start_year:end_year)){
  for(j in 1:length(0:90))
    if(dtx[i,j] == 0){
      dtx[i,j] <- dtx[i,j-1]
    }
}

etx <- dtx / matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
etx[is.infinite(etx)] <- 0
etx[is.na(etx)] <- 0
for(i in 1:length(start_year:end_year)){
  for(j in 1:length(0:90))
    if(etx[i,j] == 0){
      etx[i,j] <- etx[i,j-1]
    }
}

mtx <- matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
for(i in 1:length(start_year:end_year)){
  for(j in 1:length(0:90))
    if(mtx[i,j] == 0){
      mtx[i,j] <- mtx[i,j-1]
    }
}

LCfit701_1944_fem <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
LCfit701_1952_fem <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
LCfit701_1971_fem <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))


## 2.2.1 first enlighten by auto.arima over all years

list6 <- list()

for (k in 1876:2000) {
  start_year <- k
  years <- start_year:end_year
  ages <- 0:90
  females_est <- filter(females, Year <= end_year & Year >= k)
  females_est_sub  <- subset(females_est, Age <= 90)
  dtx <- matrix(females_est_sub$dx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(dtx[i,j] == 0){
        dtx[i,j] <- dtx[i,j-1]
      }
  }
  
  etx <- dtx / matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
  etx[is.infinite(etx)] <- 0
  etx[is.na(etx)] <- 0
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(etx[i,j] == 0){
        etx[i,j] <- etx[i,j-1]
      }
  }
  
  mtx <- matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
  for(i in 1:length(k:end_year)){
    for(j in 1:length(0:90))
      if(mtx[i,j] == 0){
        mtx[i,j] <- mtx[i,j-1]
      }
  }
  
  LCfit701 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))
  
  time_series <- auto.arima(LCfit701$kappa2)
  
  list6 <- append(list6, arimaorder(time_series))
}
auto_model_fem <- data.frame(t(matrix(list6, nrow = 3, byrow = TRUE)))
names(auto_model_fem) <- c("p", "d", "q")
rownames(auto_model_fem) <- c(1876:2000)

auto_model_fem[,1] <- as.numeric(auto_model_fem[,1])
auto_model_fem[,2] <- as.numeric(auto_model_fem[,2])
auto_model_fem[,3] <- as.numeric(auto_model_fem[,3])

auto_model_fem
#mostly suggests ARIMA models over all starting years: p,d,q <= 2


## 2.2.2 d = 0

sarima(LCfit701_1944_fem$kappa2,0,0,1)#terrible
sarima(LCfit701_1944_fem$kappa2,0,0,2)#terrible
sarima(LCfit701_1944_fem$kappa2,1,0,0)#terrible
sarima(LCfit701_1944_fem$kappa2,1,0,1)#terrible
sarima(LCfit701_1944_fem$kappa2,1,0,2)#....
Arima(LCfit701_1944_fem$kappa2, order = c(1,0,2), include.drift = TRUE)$aic #[1] 363.6355
sarima(LCfit701_1944_fem$kappa2,2,0,0)#terrible
sarima(LCfit701_1944_fem$kappa2,2,0,1)#terrible
sarima(LCfit701_1944_fem$kappa2,2,0,2)#terrible

sarima(LCfit701_1952_fem$kappa2,0,0,1)#terrible
sarima(LCfit701_1952_fem$kappa2,0,0,2)#terrible
sarima(LCfit701_1952_fem$kappa2,1,0,0)#terrible
sarima(LCfit701_1952_fem$kappa2,1,0,1)#....
sarima(LCfit701_1952_fem$kappa2,1,0,2)#....
sarima(LCfit701_1952_fem$kappa2,2,0,0)#terrible
sarima(LCfit701_1952_fem$kappa2,2,0,1)#terrible
sarima(LCfit701_1952_fem$kappa2,2,0,2)#terrible

sarima(LCfit701_1971_fem$kappa2,0,0,1)#terrible
sarima(LCfit701_1971_fem$kappa2,0,0,2)#terrible
sarima(LCfit701_1971_fem$kappa2,1,0,0)#terrible
sarima(LCfit701_1971_fem$kappa2,1,0,1)#....
sarima(LCfit701_1971_fem$kappa2,1,0,2)#....
sarima(LCfit701_1971_fem$kappa2,2,0,0)#terrible
sarima(LCfit701_1971_fem$kappa2,2,0,1)#....
sarima(LCfit701_1971_fem$kappa2,2,0,2)#....

## 2.2.3 d = 1
sarima(LCfit701_1944_fem$kappa2,0,1,0)#terrible
sarima(LCfit701_1944_fem$kappa2,0,1,1)#nice
sarima(LCfit701_1944_fem$kappa2,0,1,2)#nice
sarima(LCfit701_1944_fem$kappa2,1,1,0)#nice
sarima(LCfit701_1944_fem$kappa2,1,1,1)#nice
sarima(LCfit701_1944_fem$kappa2,1,1,2)#nice
sarima(LCfit701_1944_fem$kappa2,2,1,0)#nice
sarima(LCfit701_1944_fem$kappa2,2,1,1)#nice
Arima(LCfit701_1944_fem$kappa2, order = c(2,1,1), include.drift = TRUE)$aic #[1] 361.5096
sarima(LCfit701_1944_fem$kappa2,2,1,2)#nice

sarima(LCfit701_1952_fem$kappa2,0,1,0)#not good
sarima(LCfit701_1952_fem$kappa2,0,1,1)#not good
sarima(LCfit701_1952_fem$kappa2,0,1,2)#not good
sarima(LCfit701_1952_fem$kappa2,1,1,0)#not good
sarima(LCfit701_1952_fem$kappa2,1,1,1)#not good
sarima(LCfit701_1952_fem$kappa2,1,1,2)#better
sarima(LCfit701_1952_fem$kappa2,2,1,0)#better
sarima(LCfit701_1952_fem$kappa2,2,1,1)#better
Arima(LCfit701_1952_fem$kappa2, order = c(2,1,1), include.drift = TRUE)$aic #[1] 298.4103
sarima(LCfit701_1952_fem$kappa2,2,1,2)#better than d=0, but worse than 1944 at d=1 in overall

sarima(LCfit701_1971_fem$kappa2,0,1,0)#
sarima(LCfit701_1971_fem$kappa2,0,1,1)#
sarima(LCfit701_1971_fem$kappa2,0,1,2)#
sarima(LCfit701_1971_fem$kappa2,1,1,0)#
sarima(LCfit701_1971_fem$kappa2,1,1,1)#
sarima(LCfit701_1971_fem$kappa2,1,1,2)#
sarima(LCfit701_1971_fem$kappa2,2,1,0)#
Arima(LCfit701_1971_fem$kappa2, order = c(2,1,0), include.drift = TRUE)$aic #[1] 212.0876
sarima(LCfit701_1971_fem$kappa2,2,1,1)#best one among these
Arima(LCfit701_1971_fem$kappa2, order = c(2,1,1), include.drift = TRUE)$aic #[1] 210.4082
sarima(LCfit701_1971_fem$kappa2,2,1,2)#


## 2.2.3 d = 2
sarima(LCfit701_1944_fem$kappa2,0,2,0)#terrible
sarima(LCfit701_1944_fem$kappa2,0,2,1)#terrible
sarima(LCfit701_1944_fem$kappa2,0,2,2)#not bad
sarima(LCfit701_1944_fem$kappa2,1,2,0)#not good
sarima(LCfit701_1944_fem$kappa2,1,2,1)#not good
sarima(LCfit701_1944_fem$kappa2,1,2,2)#not bad
sarima(LCfit701_1944_fem$kappa2,2,2,0)#better
sarima(LCfit701_1944_fem$kappa2,2,2,1)#nice
Arima(LCfit701_1944_fem$kappa2, order = c(2,2,1), include.drift = TRUE)$aic #[1] 362.1141
sarima(LCfit701_1944_fem$kappa2,2,2,2)#not bad

sarima(LCfit701_1952_fem$kappa2,0,2,0)#terrible
sarima(LCfit701_1952_fem$kappa2,0,2,1)#terrible
sarima(LCfit701_1952_fem$kappa2,0,2,2)#worse
sarima(LCfit701_1952_fem$kappa2,1,2,0)#terrible
sarima(LCfit701_1952_fem$kappa2,1,2,1)#terrible
sarima(LCfit701_1952_fem$kappa2,1,2,2)#terrible
sarima(LCfit701_1952_fem$kappa2,2,2,0)#worse
sarima(LCfit701_1952_fem$kappa2,2,2,1)#worse
sarima(LCfit701_1952_fem$kappa2,2,2,2)#best among these
Arima(LCfit701_1952_fem$kappa2, order = c(2,2,2), include.drift = TRUE)$aic #[1] 300.1776

sarima(LCfit701_1971_fem$kappa2,0,2,0)#terrible
sarima(LCfit701_1971_fem$kappa2,0,2,1)#terrible
sarima(LCfit701_1971_fem$kappa2,0,2,2)#terrible
sarima(LCfit701_1971_fem$kappa2,1,2,0)#terrible
sarima(LCfit701_1971_fem$kappa2,1,2,1)#terrible
sarima(LCfit701_1971_fem$kappa2,1,2,2)#terrible
sarima(LCfit701_1971_fem$kappa2,2,2,0)#terrible
sarima(LCfit701_1971_fem$kappa2,2,2,1)#terrible
sarima(LCfit701_1971_fem$kappa2,2,2,2)#only one that is a bit better
Arima(LCfit701_1971_fem$kappa2, order = c(2,2,2), include.drift = TRUE)$aic #[1] 210.6566

#generally speaking, have better performance in d = 1, according to aic, 
#choose starting year =1971 and ARIMA(2,1,1) model


time_series <- Arima(LCfit701_1971_fem$kappa2, order = c(2, 1, 1), include.drift = TRUE)


#predicting
yv_spec <- years
alpha.x	  <- LCfit701_1971_fem$beta1
beta.x 	  <- LCfit701_1971_fem$beta2		  
kappa.t   <- LCfit701_1971_fem$kappa2	 	

LCfit701_1971_fem$mtx <- t(exp(rep(alpha.x, each = 1, times = length(yv_spec)) + beta.x%o%kappa.t))


source('simModels.R')

sim_LC = sim2001(
  xx = LCfit701_1971_fem$x,
  yy = LCfit701_1971_fem$y,
  beta1v = LCfit701_1971_fem$beta1,
  beta2v = LCfit701_1971_fem$beta2,
  kappa2v = LCfit701_1971_fem$kappa2,
  nsim = 10000,
  tmax = 172, # 52 + 120
  nyears = length(years)
)

names(sim_LC)

#### 2.3. life expectancy plots ####

m_rate <- sim_LC$qaa
M.tx <- LCfit701_1971_fem$mtx
repeated_A <- replicate(10001, M.tx, simplify = FALSE)
M.tx <- array(unlist(repeated_A), dim = c(dim(M.tx), length(repeated_A)))

M.tx <- aperm(M.tx, c(2, 3, 1))
testing <- aperm(m_rate, c(1, 3, 2))
C <- abind(M.tx, testing)
m_rate <- aperm(C, c(1, 3, 2))
dim(m_rate)

le_yv   <- 1971:2071
le_ages <- c(0,65)
le_type <- c("per", "coh")
dimnames(m_rate) <- list(0:90, 1971:2191, 1:10001)
le_fem <- life_exp(le_yv, le_type, le_ages, m_rate)


dimnames(dtx) <- list(1971:2019, 0:90)
dimnames(etx) <- list(1971:2019, 0:90)

kannisto_nages <- 30
kannisto_nobs  <- 11
m_obs_cl <- close_obs_death_rates(dtx, etx, kannisto_nages, kannisto_nobs)
le_e0_19 <- plot_life_exp(le_yv, 0, le_fem, "Female", le_type, "99%", m_obs = m_obs_cl)
le_e65_19 <- plot_life_exp(le_yv, 65, le_fem, "Female", le_type, "99%", m_obs = m_obs_cl)
le_e0_19
le_e65_19




#### 2.3. Residual plots ####
grid <- expand.grid(period = years, age = ages)
grid$res <- as.vector(LCfit701_1971_fem$epsilon)
names(grid) <- c("Year", "Age", "Residual")


## ------------------------------------------------------------------------------------
head(grid)

## ------------------------------------------------------------------------------------

p <- ggplot(grid, aes(x = Year, y = Age)) + geom_tile(aes(fill = Residual)) +
  scale_fill_gradientn(colours =  topo.colors(7)) +
  theme_bw() + theme(legend.position = "bottom",
                     legend.title = element_text(size = 15))
p



#### 2.4. Inspecting goodness of fit: observed mx versus fitted mx ####

age <- 25
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_25 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.00125, label = "Age 25", size = 10) +
  ggtitle("Switzerland - males, 1971 - 2019, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[25,][t])) + labs(x = bquote(Year (t)))

age <- 45
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_45 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.004, label = "Age 45", size = 10) +
  ggtitle("Switzerland - males, 1971 - 2019, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[45,][t])) + labs(x = bquote(Year (t)))

age <- 65
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_65 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.03, label = "Age 65", size = 10) +
  ggtitle("Switzerland - males, 1971 - 2019, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[65,][t])) + labs(x = bquote(Year (t)))

age <- 85
rates <- dtx/etx
df <- tibble(Year = years, obs = rates[, age - min(ages) + 1], fit = exp(LCfit701$beta1[age - min(ages) + 1] + LCfit701$beta2[age - min(ages) + 1] * LCfit701$kappa2))

g_85 <- ggplot(df) + geom_point(aes(Year, obs), col = "black") +
  geom_line(aes(Year, fit), col = "black", linetype = "dashed") +
  theme_bw() + geom_text(x = 2010, y = 0.18, label = "Age 85", size = 10) +
  ggtitle("Switzerland - males, 1971 - 2019, Lee Carter, Poisson") + 
  labs(y = bquote(hat(m)[85,][t])) + labs(x = bquote(Year (t)))

grid.arrange(g_25, g_45, g_65, g_85, ncol = 2)










#### 2.5. Forecasting with time series ####


## ------------------------------------------------------------------------------------
plot(forecast(time_series, level = c(80, 85, 95)))


#### 2.6. Forecasting with the model: longevity charts ####

source('simModels.R')

sim_LC = sim2001(
  xx = LCfit701$x,
  yy = LCfit701$y,
  beta1v = LCfit701$beta1,
  beta2v = LCfit701$beta2,
  kappa2v = LCfit701$kappa2,
  nsim = 10000,
  tmax = 50,
  nyears = length(years)
)

names(sim_LC)


## ------------------------------------------------------------------------------------



## ------------------------------------------------------------------------------------
# time series kappa_t
plot(
  LCfit701$y,
  LCfit701$kappa2,
  pch = 20,
  xlim = c(1971, 2068),
  ylim = range(c(
    range(LCfit701$kappa2), range(sim_LC$dda[, , ])
  )),
  main = bquote(paste("Projection of ", kappa[t])),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(kappa[t])
)

# fan chart
fan(sim_LC$y, sim_LC$dda[, , ], color = "red")


## ------------------------------------------------------------------------------------
age = 65
minage = min(ages)

plot(
  LCfit701$y,
  exp(LCfit701$beta1[age - minage + 1] + LCfit701$beta2[age - minage + 1] * LCfit701$kappa2),
  lwd = 2,
  col = "red",
  type = "l",
  ylim = c(0, 0.04),
  main = bquote(paste(
    "Switzerland: Lee-Carter, projection of ", mu[65](t)
  )),
  xlim = c(1971, 2068),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(mu[65](t)))
)
fan(sim_LC$y,
    exp(LCfit701$beta1[age - minage + 1] + LCfit701$beta2[age - minage + 1] * sim_LC$dda[, , ]),
    col = "red")
points(LCfit701$y, rates[, age - minage + 1], col = "black", pch = 20)


## ------------------------------------------------------------------------------------
# plot of q_{65,75,85}
ages_sel  = c(65, 75, 85)
nages = length(ages_sel)
color = c("red", "green", "blue")

# y axis is logarithmic
plot(
  c(1971, 2068),
  c(0.005, 0.2),
  type = "n",
  log = "y",
  lwd = 2,
  col = "black",
  main = bquote(paste(
    "Switzerland: projection of ", q[65](t), ", ", q[75](t), ", ", q[85](t)
  )),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(q[65](t), ", ", q[75](t), ", ", q[85](t)))
)

for (j in 1:nages) {
  lines(LCfit701$y, 1 - exp(-exp(
    LCfit701$beta1[ages_sel[j] - minage + 1] + LCfit701$beta2[ages_sel[j] - minage + 1] *
      LCfit701$kappa2
  )), col = color[j], lwd = 2)
  matlines(sim_LC$y, sim_LC$qaa[ages_sel[j] - minage + 1, , 1:1000], col =
             "grey")
}
p = seq(0.05, 0.95, 0.05)

# Construct the corresponding empirical quantiles
for(j in 1:nages) {
  q_int = matrix(data = 0,
                 nrow = 50,
                 ncol = length(p))
  for (i in 1:50)
    q_int[i, ] = quantile(sim_LC$qaa[ages_sel[j] - minage + 1, i, ], p = p)
  # highlight 0.05 and 0.95 quantile
  lines(2019:2067,
        q_int[1:50, 1],
        col = color[j],
        lty = 1,
        lwd = 2)
  lines(2019:2067,
        q_int[1:50, length(p)],
        col = color[j],
        lty = 1,
        lwd = 2)
}

## ------------------------------------------------------------------------------------
# Alternative: plot of q_{65,75,85} using fan function in simModels.R
plot(
  c(1971, 2058),
  c(0.005, 0.2),
  log = "y",
  main = bquote(paste(
    "Switzerland: projection of ", q[65](t), ", ", q[75](t), ", ", q[85](t)
  )),
  sub = "ARIMA(0,1,0), 10000 sim",
  xlab = "Year (t)",
  ylab = bquote(paste(q[65](t), ", ", q[75](t), ", ", q[85](t))),
  type = "n"
)


for (j in 1:nages) {
  points(LCfit701$y, 1 - exp(-exp(
    LCfit701$beta1[ages_sel[j] - minage + 1] + LCfit701$beta2[ages_sel[j] - minage + 1] *
      LCfit701$kappa2
  )), col = color[j], pch = 20)
  fan(sim_LC$y, sim_LC$qaa[ages_sel[j] - minage + 1, ,], color = color[j])
}

