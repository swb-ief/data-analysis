##COVID19 mumbai for pipline ##


options(warn=-1)
options(message=-1)
#install.packages("drat")
#drat:::add("epiforecasts")
#install.packages("rstan")
#install.packages("EpiNow2")


suppressMessages(library("EpiNow2"))
suppressMessages(library("rstan"))
suppressMessages(library(EpiEstim))
suppressMessages(library(ggplot2))
suppressMessages(library("gridExtra"))
suppressMessages(library(incidence))
suppressMessages(library(magrittr))
suppressMessages(library(readr))# for read_csv
suppressMessages(library(knitr)) # for kable
suppressMessages(library(dplyr))
#myfile <- "https://raw.githubusercontent.com/saurabhmj/etl-pipeline/metrics-branch/data_pipeline/output/metrics.csv"
#suppressMessages(mumbai<-read_csv(myfile))

                                                                                                                                        
mumbai_recent <- read.csv("ME.csv", header=T)

mumbai_working_tab <- data.frame(date= as.Date(mumbai_recent$as_of,  "%m/%d/%y", origin="01/01/70"), mumbai_recent[,c(3:6)])
mumbai_working_tab <- mumbai_working_tab[order(mumbai_working_tab[,1]), ]

delta_case_d <- mumbai_working_tab[,2] + mumbai_working_tab[,5]  - mumbai_working_tab[,3] -mumbai_working_tab[,4]

plot(delta_case_d, ty="o") ##check how the series looks like##

mumbai_recent_tab <- data.frame(date=as.Date(mumbai_working_tab[,1], origin="1970-01-01"), confirm=delta_case_d)

##model parameters as default##
##note that parameters about generation_time, incubation_period, reporting_delay are 
reporting_delay <- EpiNow2::bootstrapped_dist_fit(rlnorm(100, log(6), 1))
## Set max allowed delay to 30 days to truncate computation
reporting_delay$max <- 30
generation_time <- list(mean = EpiNow2::covid_generation_times[1, ]$mean,
                        mean_sd = EpiNow2::covid_generation_times[1, ]$mean_sd,
                        sd = EpiNow2::covid_generation_times[1, ]$sd,
                        sd_sd = EpiNow2::covid_generation_times[1, ]$sd_sd,
                        max = 30)
incubation_period <- list(mean = EpiNow2::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow2::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)


estimates_mumbai_recent <- EpiNow2::epinow(reported_cases = mumbai_recent_tab, generation_time = generation_time,
                                           delays = list(incubation_period, reporting_delay), horizon = 7, samples = 1000, 
                                           warmup = 200, cores = 4, chains = 4, verbose = TRUE, adapt_delta = 0.95)

Rt_EpiNow2 <- estimates_mumbai_recent$estimates$summarised

Rt_mean_sd <- Rt_EpiNow2[which(Rt_EpiNow2[,"variable"]=="R" & Rt_EpiNow2[,"type"]=="estimate") ,c(10,11)]


Rt_tab = cbind(Rt_EpiNow2[which(Rt_EpiNow2[,"variable"]=="R" &  Rt_EpiNow2[,"type"]=="estimate") ,c(1, 9, 7,8)], mean=Rt_mean_sd[,1],
     CI_lower = Rt_mean_sd[,1]-1.96* Rt_mean_sd[,2],  CI_upper =Rt_mean_sd[,1]+ 1.96* Rt_mean_sd[,2]) 

Rt_tab

write.csv(Rt_tab, "Rt_mumbai_0823_1115.csv")

#################################################
#################################################
####ignore below as of 11.10.2020################
#################################################
#################################################

##set date 
#kable(head(mumbai))

case_series_mumbai<-as.numeric(unlist(mumbai[,7])) ## take out delta case##
tot_cases_mumbai<-as.numeric(unlist(mumbai[,4])) ## take out delta case##
case_dates_mumbai <- unlist(mumbai[,1])


#length(case_series_mumbai)

mumbai_tab <- data.frame(date= as.Date(case_dates_mumbai,  origin = "1970-01-01"), confirm=case_series_mumbai)


mumbai_tab2 <- mumbai_tab[-1,]
mumbai_tab3 <- data.frame(date= as.Date(case_dates_mumbai,  origin = "1970-01-01"), tot_cases=tot_cases_mumbai)

##this part is from {incidence}##
mumbai_tab2$dates.x <- (case_dates_mumbai[-1] -  case_dates_mumbai[-length(case_dates_mumbai)])/2
lm1 <- stats::lm(log(confirm) ~ dates.x, data = mumbai_tab2)

r <- stats::coef(lm1)["dates.x"]
r.conf <- stats::confint(lm1, "dates.x", 0.95)
new.data <- data.frame(dates.x = sort(unique(lm1$model$dates.x)))
pred     <- exp(stats::predict(lm1, newdata = new.data, interval = "confidence", level = 0.95))
pred <- cbind.data.frame(new.data, pred)
info_list <- list(
  tab = round(c(r = r,
  r.conf = r.conf,
  doubling = log(2) / r,
  doubling.conf = log(2) / r.conf),4),
  pred = pred
)
#info_list


##this part is from 
dbl_timr <- function(data, end_date = NULL, time = 7) {
  
  if (is.null(end_date)) {
    end_date <- max(data$date)
  }
  
  start <-  data %>% filter(date == as.Date(as.Date(end_date, origin="1970-01-01") - time)) %>% pull(tot_cases)
  
  if (length(start) == 0) {
    NA
  } else if (start == 0) {
    NA
  } else {
    end   <- data %>% filter(date == as.Date(end_date, origin="1970-01-01")) %>% pull(tot_cases)
    
    r <- ((end - start) / start) * 100
    
    dt <- time * (log(2) / log(1 + (r / 100)))
    return(c(r=r, dt=dt))
  }
}

dbl_times <- NA

  
  tmp_v     <- matrix(NA, ncol=2, nrow=length(case_dates_mumbai))
    for(j in seq_along(case_dates_mumbai)) {
      task <- dbl_timr(data = mumbai_tab3, end_date = case_dates_mumbai[j], time = 7)
   if(is.na(task)==T) {
     tmp_v[j, ] <-c(NA, NA)
   } else {
     tmp_v[j, ]  <- task
    }
              
  }
  
  colnames(tmp_v) <- c("r", "doubling time")
  dt_mumbai <-data.frame(date=as.Date(case_dates_mumbai, origin="1970-01-01"), tmp_v)
  dt_mumbai <-dt_mumbai[is.na(dt_mumbai[,2])==F, ]

  tab_dt_mumbai <- c(r = mean(dt_mumbai[,2]/100), r_CI = c(mean(dt_mumbai[,2]/100) + qnorm(0.025)*sd(dt_mumbai[,2]/100), mean(dt_mumbai[,2]/100) + qnorm(1-0.025)*sd(dt_mumbai[,2]/100)),
                     doubling_time = mean(dt_mumbai[,3]), dt_CI = c(mean(dt_mumbai[,3]) + qnorm(0.025)*sd(dt_mumbai[,3]), mean(dt_mumbai[,3]) + qnorm(1-0.025)*sd(dt_mumbai[,3])))



##old Rt: EpiEstim##
t_start <- seq(6, 87 - 6)
t_end   <- t_start + 6

Rt_covid_mumbai_0 <- EpiEstim::estimate_R(incid = case_series_mumbai, method = "parametric_si",
                                                config = make_config(list(mean_si = 3.96, std_si = 4.75, si_parametric_distr = "G",
                                                                          t_start = t_start, t_end = t_end, seed = 123)))

Rt_covid_mumbai_1 <- EpiEstim::estimate_R(incid = case_series_mumbai, method = "parametric_si",
                                          config = make_config(list(mean_si = 2, std_si = 4.75, si_parametric_distr = "G",
                                                                    t_start = t_start, t_end = t_end, seed = 123)))

Rt_covid_mumbai_2 <- EpiEstim::estimate_R(incid = case_series_mumbai, method = "parametric_si",
                                          config = make_config(list(mean_si = 8, std_si = 4.75, si_parametric_distr = "G",
                                                                    t_start = t_start, t_end = t_end, seed = 123)))


Rt_covid_mumbai_3 <- EpiEstim::estimate_R(incid = case_series_mumbai, method = "parametric_si",
                                          config = make_config(list(mean_si = 4, std_si = 10, si_parametric_distr = "G",
                                                                    t_start = t_start, t_end = t_end, seed = 123)))

Rt_covid_mumbai_4 <- EpiEstim::estimate_R(incid = case_series_mumbai, method = "parametric_si",
                                          config = make_config(list(mean_si = 4, std_si = 2.5, si_parametric_distr = "G",
                                                                    t_start = t_start, t_end = t_end, seed = 123)))

write.csv(cbind(mumbai_tab[unlist(Rt_covid_mumbai_0$R[2]),1], round(cbind(Rt_covid_mumbai_0$R[,c (8, 5, 11)], Rt_covid_mumbai_1$R[,c (8, 5, 11)],
      Rt_covid_mumbai_2$R[,c (8, 5, 11)], Rt_covid_mumbai_3$R[,c (8, 5, 11)], Rt_covid_mumbai_4$R[,c (8, 5, 11)]), 3)), "Epiestim_1010.csv")


#plot(Rt_covid_mumbai) #see the result##

##R_sim_CI <- sample_posterior_R(Rt_covid19_mumbai, n = 10000, window=77:81) ##need to fit model moving window##


generation_time2 <- list(mean = 3.96,
                        mean_sd = EpiNow2::covid_generation_times[1, ]$mean_sd,
                        sd = 4.75,
                        sd_sd = EpiNow2::covid_generation_times[1, ]$sd_sd,
                        max = 30)







estimates_mumbai <- EpiNow2::epinow(reported_cases = mumbai_tab, generation_time = generation_time,
                             delays = list(incubation_period, reporting_delay), horizon = 7, samples = 1000, 
                             warmup = 200, cores = 4, chains = 4, verbose = TRUE, adapt_delta = 0.95)
estimates_mumbai2 <- EpiNow2::epinow(reported_cases = mumbai_tab, generation_time = generation_time2,
                                    delays = list(incubation_period, reporting_delay), horizon = 7, samples = 1000, 
                                    warmup = 200, cores = 4, chains = 4, verbose = TRUE, adapt_delta = 0.95)

estimates_mumbai_recent <- EpiNow2::epinow(reported_cases = mumbai_recent_tab, generation_time = generation_time,
                                     delays = list(incubation_period, reporting_delay), horizon = 7, samples = 1000, 
                                     warmup = 200, cores = 4, chains = 4, verbose = TRUE, adapt_delta = 0.95)

Rt_EpiNow2 <- estimates_mumbai_recent$estimates$summarised[which(estimates_mumbai$estimates$summarised[,"variable"]=="R" & estimates_mumbai$estimates$summarised[,"type"]=="estimate"),]
Rt_EpiNow2 <- Rt_EpiNow2[which(unlist(Rt_EpiNow2[,1]) %in% unlist(Rt_Epiestim[,1])) ,c(1, 9, 7,8)] 


##to see the result##
##estimates_mumbai$summary
###estimates_mumbai$plot

#compare result##
Rt_Epiestim <- cbind(mumbai_tab[unlist(Rt_covid_mumbai$R[ 2]),1],Rt_covid_mumbai$R[,c (8, 5, 11)])

Rt_EpiNow2 <- estimates_mumbai$estimates$summarised[which(estimates_mumbai$estimates$summarised[,"variable"]=="R" & estimates_mumbai$estimates$summarised[,"type"]=="estimate"),]
Rt_EpiNow2 <- Rt_EpiNow2[which(unlist(Rt_EpiNow2[,1]) %in% unlist(Rt_Epiestim[,1])) ,c(1, 9, 7,8)] 


Rt_EpiNow2.2 <- estimates_mumbai2$estimates$summarised[which(estimates_mumbai2$estimates$summarised[,"variable"]=="R" & estimates_mumbai2$estimates$summarised[,"type"]=="estimate"),]
Rt_EpiNow2.2 <- Rt_EpiNow2.2[which(unlist(Rt_EpiNow2.2[,1]) %in% unlist(Rt_Epiestim[,1])) ,c(1, 9, 7,8)] 


tab_Rt <- cbind(Rt_Epiestim, Rt_EpiNow2[,-1], Rt_EpiNow2.2[,-1])
colnames(tab_Rt) <- c("date", "R_med_EpiEstim", "R_low_EpiEstim",  "R_up_EpiEstim", 
                   "R_med_EpiNow2", "R_low_EpiNow2",  "R_up_EpiNow2",
                   "R_med_EpiNow2.2", "R_low_EpiNow2.2",  "R_up_EpiNow2.2")
tab_Rt

write.csv(tab_Rt, "EpiNow2_Rt.csv")
tab_dt <- rbind(incidence = info_list$tab, Epinow = c(unlist(estimates_mumbai$summary[4,]$numeric_estimate)[1:3], unlist(estimates_mumbai$summary[5,]$numeric_estimate)),
               covid1i_india =   tab_dt_mumbai )
colnames(tab_dt) <- c("r", "r_low",  "r_up",  "doubling time", "dt_low",  "dt_up")

tab_dt
