setwd("C:/Users/r0823957/Downloads")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Advanced Life insurance mathematics: #
#        Assignment - Jing Yang        #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

###### 0. Settings ######
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
library(spaMM)
library(foreach)

source('fitModels.R')
source('simModels.R')
source('life_exp.R')

females<- read.table("females.csv", sep=",",head=TRUE)
colnames(females)[1] <- c("Year")
males<- read.table("males.csv", sep=",",head=TRUE)

####--2019, RW, original model####

start_year <- 1971
end_year <- 2019
years <- start_year:2019
ages <- 0:90

females_est <- filter(females, Year <= 2019 & Year >= start_year)
females_est_sub  <- subset(females_est, Age <= 90)

dtx <- matrix(females_est_sub$dx, nrow = length(years), byrow = TRUE)
for(i in 1:length(start_year:2019)){
  for(j in 1:length(0:90))
    if(dtx[i,j] == 0){
      dtx[i,j] <- dtx[i,j-1]
    }
}
dim(dtx)

etx <- dtx / matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
etx[is.infinite(etx)] <- 0
etx[is.na(etx)] <- 0
for(i in 1:length(start_year:2019)){
  for(j in 1:length(0:90))
    if(etx[i,j] == 0){
      etx[i,j] <- etx[i,j-1]
    }
}
dim(etx)

mtx <- matrix(females_est_sub$mx, nrow = length(years), byrow = TRUE)
for(i in 1:length(start_year:2019)){
  for(j in 1:length(0:90))
    if(mtx[i,j] == 0){
      mtx[i,j] <- mtx[i,j-1]
    }
}
dim(mtx)

LCfit701_1971_fem <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))

time_series <- Arima(LCfit701_1971_fem$kappa2, order = c(2, 1, 1), include.drift = TRUE)

yv_spec <- years
alpha.x	  <- LCfit701_1971_fem$beta1
beta.x 	  <- LCfit701_1971_fem$beta2		  
kappa.t   <- LCfit701_1971_fem$kappa2	 	

LCfit701_1971_fem$mtx <- t(exp(rep(alpha.x, each = 1, times = length(yv_spec)) + beta.x%o%kappa.t))


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
le_fem_19 <- life_exp(le_yv, le_type, le_ages, m_rate)


dimnames(dtx) <- list(1971:2019, 0:90)
dimnames(etx) <- list(1971:2019, 0:90)

kannisto_nages <- 30
kannisto_nobs  <- 11
m_obs_cl <- close_obs_death_rates(dtx, etx, kannisto_nages, kannisto_nobs)
le_e0_19 <- plot_life_exp(le_yv, 0, le_fem_19, "Female", le_type, "99%", m_obs = m_obs_cl)


####--2020, RW####

start_year <- 1971
end_year <- 2020
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

LCfit701_1971_fem_2020 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))

yv_spec <- years
alpha.x	  <- LCfit701_1971_fem_2020$beta1
beta.x 	  <- LCfit701_1971_fem_2020$beta2		  
kappa.t   <- LCfit701_1971_fem_2020$kappa2	 	

LCfit701_1971_fem_2020$mtx <- t(exp(rep(alpha.x, each = 1, times = length(yv_spec)) + beta.x%o%kappa.t))

sim_LC = sim2001(
  xx = LCfit701_1971_fem_2020$x,
  yy = LCfit701_1971_fem_2020$y,
  beta1v = LCfit701_1971_fem_2020$beta1,
  beta2v = LCfit701_1971_fem_2020$beta2,
  kappa2v = LCfit701_1971_fem_2020$kappa2,
  nsim = 10000,
  tmax = 172, # 52 + 120
  nyears = length(years)
)


m_rate <- sim_LC$qaa
M.tx <- LCfit701_1971_fem_2020$mtx
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
dimnames(m_rate) <- list(0:90, 1971:2192, 1:10001)
le_fem_20 <- life_exp(le_yv, le_type, le_ages, m_rate)


dimnames(dtx) <- list(1971:2020, 0:90)
dimnames(etx) <- list(1971:2020, 0:90)

kannisto_nages <- 30
kannisto_nobs  <- 11
m_obs_cl <- close_obs_death_rates(dtx, etx, kannisto_nages, kannisto_nobs)
le_e0_20 <- plot_life_exp(le_yv, 0, le_fem_20, "Female", le_type, "99%", m_obs = m_obs_cl)
le_e65_20 <- plot_life_exp(le_yv, 65, le_fem_20, "Female", le_type, "99%", m_obs = m_obs_cl)
le_e0_20
le_e65_20

####--2021, RW####

start_year <- 1971
end_year <- 2021
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

LCfit701_1971_fem_2021 <- fit701(ages, years, etx, dtx, matrix(1, length(years), length(ages)))

yv_spec <- years
alpha.x	  <- LCfit701_1971_fem_2021$beta1
beta.x 	  <- LCfit701_1971_fem_2021$beta2		  
kappa.t   <- LCfit701_1971_fem_2021$kappa2	 	

LCfit701_1971_fem_2021$mtx <- t(exp(rep(alpha.x, each = 1, times = length(yv_spec)) + beta.x%o%kappa.t))

sim_LC = sim2001(
  xx = LCfit701_1971_fem_2021$x,
  yy = LCfit701_1971_fem_2021$y,
  beta1v = LCfit701_1971_fem_2021$beta1,
  beta2v = LCfit701_1971_fem_2021$beta2,
  kappa2v = LCfit701_1971_fem_2021$kappa2,
  nsim = 10000,
  tmax = 172, # 52 + 120
  nyears = length(years)
)


m_rate <- sim_LC$qaa
M.tx <- LCfit701_1971_fem_2021$mtx
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
dimnames(m_rate) <- list(0:90, 1971:2193, 1:10001)
le_fem_21 <- life_exp(le_yv, le_type, le_ages, m_rate)


dimnames(dtx) <- list(1971:2021, 0:90)
dimnames(etx) <- list(1971:2021, 0:90)

kannisto_nages <- 30
kannisto_nobs  <- 11
m_obs_cl <- close_obs_death_rates(dtx, etx, kannisto_nages, kannisto_nobs)
le_e0_21 <- plot_life_exp(le_yv, 0, le_fem_21, "Female", le_type, "99%", m_obs = m_obs_cl)
le_e65_21 <- plot_life_exp(le_yv, 65, le_fem_21, "Female", le_type, "99%", m_obs = m_obs_cl)
le_e0_21
le_e65_21

##
saveRDS(le_fem_19, file = "le_19.rds")
AA <- readRDS("le_19.rds")

