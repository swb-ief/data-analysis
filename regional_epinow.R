# using regional epinow.
# from df2, which is the loaded dataset


df2 <- df

delta_case <- df2$Confirmed + df2$Active - df2$Recovered - df2$Deceased

date <- as_date(df2$Date)

confirm <- delta_case

region <- df2$city

df_r = tibble(date, confirm, region)

df_r %>% count(region)


rt_regional <- 
  EpiNow2::regional_epinow(reported_cases = df_r, 
                  generation_time = generation_time,
                  delays = delay_opts(incubation_period, reporting_delay),
                  rt = rt_opts(prior = list(mean = 2, sd = 0.2)), 
                  stan = stan_opts(cores = 4, samples = 100),
                  verbose = TRUE,
                  output = c("regions"),
                  CrIs = 0.95)

# mumbai results

rt_regional_mumbai <- rt_regional$regional$Pune$estimates$summarised %>% tbl_df()

# pune results

rt_regional_pune <- rt_regional$regional$Pune$estimates$summarised %>% tbl_df()


# these dataframes can now be saved back to the folder that is made as the working dir.
# they can also be saved directly to google sheets.

