######This script is for VAR and VECM models for risk narratives measures
### Change the parameter values to generate different figures and tables
## Author: Yuting Chen
library(tsDyn)
library(vars)
# JOHANSEN COINTEGRATION 
library(urca)
library(forecast)
library(tidyverse)
rm(list = ls(all.names = TRUE))
library("readxl")
my_data <- read_excel("Data_VECM_RiskNarra_VIX.xlsx",guess_max = 10000)

# create a list to contain the manual inputs and print this out in the end
list_manualcheck <- list()

# convert to monthly when necessary
dailyormonthly = 1 # 1 for daily data
list_manualcheck[["dailyormonthly"]] <- dailyormonthly


# if monthly, aggregate to get monthly data, at the moment, monthly data is used only for narra-vix system
if (dailyormonthly !=1){
  library(dplyr)
  Narramonthly <- my_data %>%mutate(year = year(date), month = month(date)) %>%
    group_by(year, month) %>% summarise(total = sum(RiskNarratives))
  Narramonthly <- rename(Narramonthly, RiskNarratives = total)
  
  VIXmonthly   <- my_data %>%mutate(year = year(date), month = month(date)) %>%
    group_by(year, month)%>%summarise(total = mean(VIX,na.rm=TRUE))
  VIXmonthly   <-rename(VIXmonthly, VIX = total)
  
  SPmonthly   <- my_data %>%mutate(year = year(date), month = month(date)) %>%
    group_by(year, month)%>%summarise(total = mean(SP500,na.rm=TRUE))
  SPmonthly   <-rename(SPmonthly, SP500 = total)
  
  my_data <- data.frame(Narramonthly$month, Narramonthly$RiskNarratives, VIXmonthly$VIX, SPmonthly$SP500)
  colnames(my_data) <- c("date", "RiskNarratives", "VIX","SP500")
  
}


# drop missing rows
# my_data <- my_data[,"VIX","RiskNarratives","SP500"] # I use this line to check whether adding the exogens cause more missing rows
my_data <- my_data %>% drop_na()
# check the standar errors of the narrative measure
narra_std <- apply(my_data, 2, sd)
# standardize data to be zero mean and 1 std

my_main_data <- scale(my_data[, c("RiskNarratives","NarraVirality","SP500","VIX","Volume")])
# create the series used for analysis
risknarra <- my_main_data[,"RiskNarratives"]
narraviral <- my_main_data[,"NarraVirality"]
vix <- my_main_data[,"VIX"]
sp <- my_main_data[,"SP500"]
volume <- my_main_data[,"Volume"]
# need to convert them to timeseries object to be able to plot as a preview, check whether they are persistent
# autoplot(cbind(ts(risknarra),ts(vix)))


### create the system, can have different options and orderings
modelused <- 2 # 1 for (risknarra,vix), 2 for (risknarra,volume,sp), 3 for (narraviral,volume,sp)
list_manualcheck[["modelused"]] <- modelused
if (modelused == 1){
  dset <- cbind(risknarra,vix,sp)
}else if (modelused == 2){
  exogens <- c("VIXl1","VIXl2","VIXl3","VIXl4","VIXl5","Monday","Tuesday","Wednesday","Thursday","January")
  exogens <- my_data[,exogens]
  dset <- cbind(risknarra,volume,sp) # I also get the results of (risknarra,sp), (risknarra,volume)
}else if (modelused == 3){
  exogens <- c("VIXl1","VIXl2","VIXl3","VIXl4","VIXl5","Monday","Tuesday","Wednesday","Thursday","January")
  exogens <- my_data[,exogens]
  dset <- cbind(narraviral,volume,sp) # I also get the results of (risknarra,sp), (risknarra,volume)
}

# dset <- cbind(vix,risknarra)
# 


### select the optimal number of lags, it can be used for both VAR and VECM
if (dailyormonthly == 1){
  lagmax <- 5
}else{
  lagmax <- 6
}
list_manualcheck[["lagmax"]] <- lagmax

lagselect <- VARselect(dset,lag.max = lagmax, type = "const", exogen = exogens)
selected_lags <- lagselect$selection
# define a function to get the mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# get the mode of the number of lags.
selected_lag <- getmode(selected_lags) - 1 # need to have "minus one" for VECM, because it's difference
list_manualcheck[["selected_lag"]] <- selected_lag # sometimes I will select the number arbitrarily

print(selected_lags)
print(selected_lag)
# lagselect$criteria # uncomment this if I m interested in the result details


# stationarity test first, if all of them pass the tests, then no need for cointegration test, no need for VECM
library(tseries)
## adf test requires the declaration of the number of lags
# adf.test(risknarra)
# adf.test(sp)
# adf.test(vix)
# pp.test(risknarra)
# pp.test(sp)
# pp.test(volume)
# pp.test(vix)

stationaryornot = 1 # one if all of the variables are stationary both theoritically or based on tests
list_manualcheck[["stationaryornot"]] <- stationaryornot 

# if not stationary, do VECM
if (stationaryornot != 1){
  cointest_trace <- ca.jo(dset, type = "trace", ecdet = "const", K = selected_lag)
  summary(cointest_trace)
  cointest_eigen <- ca.jo(dset, type = "eigen", ecdet = "const", K = selected_lag)
  summary(cointest_eigen)
  
  # check the reslts and type to change the number of cointegrations manually
  num_coint = 1
  list_manualcheck[["num_coint"]] <- num_coint 
  
  # build the VECM model
  model1 <- VECM(dset,selected_lag,r=num_coint,estim = ("2OLS"))
  summary(model1)
  
  # change vecm to var
  model1var <- vec2var(cointest_trace,r = num_coint)
  
  # IRF
  nahead <- 10
  list_manualcheck[["nahead"]] <- nahead
  IRFvix <- irf(model1var, impulse = "vix", response = "risknarra", n.ahead =nahead, boot = TRUE)
  plot(IRFvix,ylab = "Risk Narrative Intensity", main = "VIX's shock to Risk Narrative Intensity")
  IRFnarra <- irf(model1var, impulse = "risknarra", response = "vix", n.ahead =nahead, boot = TRUE)
  plot(IRFnarra,ylab = "VIX", main = "Risk Narratives' shock to VIX")
  # IRFsp <- irf(model1var, impulse = "risknarra", response = "sp", n.ahead =25, boot = TRUE)
  # plot(IRFsp,ylab = "Risk Narrative Intensity", main = "SP500' shock to Risk Narrative Intensity")
  
  # variance decomposition
  FEVD <- fevd(model1var,n.ahead = 10)
  plot(FEVD)
  
  # standard VAR for granger causality test
  model2 <- VAR(dset,p= selected_lag+1, type = "const", season = NULL, exog = NULL)
  
}else{
  # else, standard VAR, 
  # first, a simple ols regression to check the relationship
  ols1 <- lm(vix~risknarra)
  summary(ols1)
  ols2 <- lm(vix~ lag(risknarra))
  summary(ols2)
  
  # standard VAR 
  model2 <- VAR(dset,p= selected_lag+1, type = "const", season = NULL, exogen = exogens)
  
  # IRFs
  if(dailyormonthly == 1){
    nahead <- 15
  }else{
    nahead <- 15
  }
  
  list_manualcheck[["nahead"]] <- nahead
  
  if (modelused == 1){
    variable1 <- "vix"
    variable2 <- "risknarra"
    ylab1 <- "Risk Narrative Intensity"
    ylab2 <- "VIX"
    main1 <- "VIX's shock to Risk Narrative Intensity"
    main2 <- "Risk Narratives' shock to VIX"
  }else if(modelused ==2){
    variable1 <- "sp"
    variable2 <- "risknarra"
    ylab1 <- "Risk Narrative Intensity"
    ylab2 <- "S&P500 Returns"
    main1 <- "Returns' shock to Risk Narrative Intensity"
    main2 <- "Risk Narrative Intensity's shock to Returns"
  }else{
    variable1 <- "sp"
    variable2 <- "narraviral"
    ylab1 <- "Risk Narrative Virality"
    ylab2 <- "S&P500 Returns"
    main1 <- "Returns' shock to Risk Narrative Virality"
    main2 <- "Risk Narrative Virality's shock to Returns"
  }
  
  if (dailyormonthly == 1){
    figname =gsub(" ", "", paste("/Output/IRF",variable1,"to",variable2,".png"))
  }else{
    figname =gsub(" ", "", paste("/Output/IRF",variable1,"to",variable2,"_monthly.png"))
  }
  IRFvix <- irf(model2, impulse = variable1, response = variable2, n.ahead =nahead, boot = TRUE)
  png(figname, width = 800, height = 350)
  plot(IRFvix,ylab = ylab1, main = main1, sub = NA,mar = c(2, 5, 0, 5),oma=c(2,0,3,0),cex.main=3)
  dev.off()
  
  if (dailyormonthly == 1){
    figname =gsub(" ", "", paste("/Output/IRF",variable2,"to",variable1,".png"))
  }else{
    figname =gsub(" ", "", paste("/Output/IRF",variable2,"to",variable1,"_monthly.png"))
  }
  IRFnarra <- irf(model2, impulse = variable2, response = variable1, n.ahead =nahead, boot = TRUE)
  png(figname, width = 800, height = 350)
  plot(IRFnarra,ylab = ylab2, main = main2,sub = NA,mar = c(2, 5, 0, 5),oma=c(2,0,3,0),cex.main=3)
  dev.off()
  
  ### Combine the figures using a python program
  
}

# Granger causality
# coeftest(model2, vcov = NeweyWest(model2, lag = selected_lag))
Grangervix <- causality(model2, cause = variable1, vcov.=NeweyWest(model2, lag = selected_lag)) # vcov.=vcovHC(model2)
Grangervix
Grangernarra <- causality(model2, cause = variable2, vcov.=NeweyWest(model2, lag = selected_lag))
Grangernarra

# the coeffients
coeftest(model2, vcov. = sandwich ) #, vcov. = sandwich 


# check the inputs
list_manualcheck

