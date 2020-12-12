
################################################################
##                    Multiple City-Rscript                    #
################################################################

# COVID 19 multiple city pipeline code
# Aim:
# 1. Convert code for V 1.3
# 2. Change code to allow for multiple cities.

# first load all the libraries needed.

options(warn = -1)
options(message = -1)
# install.packages("drat")
# drat:::add("epiforecasts")
# install.packages("rstan")
# install.packages("EpiNow2")

suppressMessages(library("lubridate"))
suppressMessages(library("tidyverse"))
suppressMessages(library("EpiNow2"))
suppressMessages(library("rstan"))
suppressMessages(library(EpiEstim))
suppressMessages(library(ggplot2))
suppressMessages(library("gridExtra"))
suppressMessages(library(incidence))
suppressMessages(library(magrittr))
suppressMessages(library(readr)) # for read_csv
suppressMessages(library(knitr)) # for kable
suppressMessages(library(readxl))


# to make the code acceptable for multiple
# cites, we need to keep the name of the df constant
# and then keep a placeholder for that specific city.
# then run the code where x will be replaced by the name of the city.
# column names need to be used for the code, but need to be kept the same always.
# am going to use the colnames that are available in the mumbai df.
# to check the code I have created a toy dataset where the two cities are mumbai and pune.
# I have randomly assigned the two cities into another column - city.
# this code can be used to calculate rt and dt for each city separately
# need to provide value for x and then run the code without any changes.
# we can provide a path to the folder with the city name to make sure that the
# results are saved in that folder.


# if we setwd for each time that we run the Rcode, then we will 
# not need to decide where the file will be saved.
# file will be saved directly to that working dir.

# setwd() --- we can set the working dir and then the 
# results will automatically be saved here.





# decide city name.

x <- cityname # insert city name here without " "

x <- 'Mumbai'

# load the dataframe.
# df <- read_csv("path/here.csv") #


df = 
read_excel('https://docs.google.com/spreadsheets/d/1HeTZKEXtSYFDNKmVEcRmF573k2ZraDb6DzgCOSXI0f0/edit#gid=0.xlsx',
sheet = 'city_stats')

# filter to keep data from only city of interest.
# here in the toy dataset have created a col - city with two values - mumbai , pune.
# we can change this small code according to the data-set colname that we are actually using

df2 <- df %>% filter(city == x) # now to only keep city of interest.

# check the df2 

glimpse(df2)

delta_case <- df2$Confirmed + df2$Active - df2$Recovered - df2$Deceased
date <- as_date(df2$Date)

confirm <- delta_case

df3 <- tibble(date, confirm) # the colnames need to be
# exactly this for the epinow function to work.

glimpse(df3)


# now the 2 columns needed for the Rt calculations are decided.
# to make the Rt calculation for that cityname
# need to get the generation_time and then the incubation_period.


# get the generation and incubation time from the new EpiNow2 package.

generation_time <- 
get_generation_time(disease = "SARS-CoV-2", source = "ganyani")


incubation_period <- 
get_incubation_period(disease = "SARS-CoV-2", source = "lauer")




## model parameters as default##
## note that parameters about generation_time,
# incubation_period, and reporting_delay are set as default in the package.

reporting_delay <- EpiNow2::bootstrapped_dist_fit(rlnorm(100, log(6), 1))

## Set max allowed delay to 30 days to truncate computation

reporting_delay$max <- 30

# values for generation time and incubation period have been defined now.
# the code below is for v 1.3.0 package.
# set credible interval as 0.95 


 rt <- 
  EpiNow2::epinow(reported_cases = df3, 
  generation_time = generation_time,
  delays = delay_opts(incubation_period, reporting_delay),
  rt = rt_opts(prior = list(mean = 2, sd = 0.2)), 
  stan = stan_opts(cores = 4, samples = 100),
  verbose = TRUE,
  CrIs = 0.95)


#
# get the summary estimates with the credible intervals.

rt <- summary(rt, type = "parameters", params = "R")  


# this is the summary estimate for the city specified in the beginning.
# here I have pasted the results back to my folder, but we can paste them back to google sheets
# directly.


write_csv(rt, "E:/rt.csv")


# doubling time function does not depend upon any package.
# so that can be used it is. 
# paste the doubling time function here first...


compute_doubling_time <- function(total_cases, case_dates, time.gap, alpha=0.05){
  suppressMessages(require(dplyr))
  
  data_tab = data.frame(date = case_dates, tot_cases = total_cases )
  
  delta_case = data_tab[-1,2] - data_tab[-dim(data_tab)[1],2] 
  
  dat = data_tab
  
  t.gap = time.gap
  
  #	dbl_timr <- function(dat,  t.gap = time.gap) {
  
  #if (is.null(end_date)) {
  #    end_date <- max(dat$date)
  #  }
  
  #t.start <-  dat %>% filter(date == as.Date(as.Date(end_date, origin="1970-01-01") - t.gap)) %>% pull(tot_cases)
  #n = length(data$date)
  
  #t.start = as.Date(data$date[-seq(n-time + 1, n)], origin="1970-01-01")
  #  if (length(t.start) == 0) {
  #    NA
  # } else if (t.start == 0) {
  #    NA
  #  } else {
  #    t.end   <- data %>% filter(date == as.Date(end_date, origin="1970-01-01")) %>% pull(tot_cases)
  #t.end <- as.Date(data$date[-seq(1, time)], origin="1970-01-01")
  # }
  
  
  end.time   <- dat$date + t.gap
  
  end.time   <- end.time[which(end.time %in% dat$date)]
  
  t.end   <- dat$tot_cases[which(dat$date %in%  end.time)]
  
  start.time <- dat$date[seq(1, length(t.end))]
  
  t.start <- dat$tot_cases[seq(1, length(t.end))]
  
  if(length(t.start) != length(t.end)){
    message("check the date")
    break
  }
  
  
  r <- ((t.end - t.start) / t.start) 
  dt <- time.gap * (log(2) / log(1 + (r)))
  
  r_d <- (dat$tot_cases[-1] - dat$tot_cases[-length(dat$tot_cases)] )/dat$tot_cases[-length(dat$tot_cases)]
  dt_d <-  (log(2) / log(1 + (r_d)))
  
  sd_r <-c()
  sd_dt <-c()
  for(t in 1:(length(t.start)-1)){
    
    sd_r <- c(sd_r, sd(r_d[which(dat$date %in% seq(start.time[t], end.time[t], 1))]))
    sd_dt <- c(sd_dt, sd(dt_d[which(dat$date %in% seq(start.time[t], end.time[t], 1))]))
    
  }
  sd_r <- c(sd_r, sd_r[(length(t.start)-1)])
  sd_dt <- c(sd_dt, sd_dt[(length(t.start)-1)])
  
  
  return(data.frame(date=as.Date(end.time, origin="1970-01-01"),r=r, r_ci_low = r + qnorm(alpha/2)*sd_r, r_ci_up = r + qnorm(1-alpha/2)*sd_r,
                    dt=dt, dt_ci_low = dt + qnorm(alpha/2)*sd_dt, dt_ci_up = dt + qnorm(1-alpha/2)*sd_dt))
  
}


total_cases <- df3$confirm
cases_dates <- df3$date



db <- compute_doubling_time(total_cases, cases_dates, time.gap = 7, alpha = 0.95)


write_csv(db, "E:/db.csv")

# now rt and db are the results for that city which is provided in x.
# rt and db are saved.
# we can create a path at the beginning of the code which can the decide where we
# want the results for both to be saved.

