#analysis for predicting lead at the individual location level

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(readr)

#####loading and cleaning####

ws <- read_csv("~/Watershed_aggregagted_attributes.csv")%>%
  select("site", "pb_mean", "area_m2", "tss_med", "age_mode", "age_mean", "ws_cpy", "str_cpy", "str_dens", "pop_dens", "aadt", "res")

sites <- read_csv("~/stormwater_conc_lead_v_1.2.csv")%>%
  select("site", "month", "year", "pb_ppb") %>%
  filter(pb_ppb != "NA")%>%
  filter(site != "mprb-9")%>%
  filter(site != "crwd-40")

#combine

sites_ws <- left_join(sites, ws, by = "site")


####model comparison####

#scale variables and filter

sites_2 <- sites_ws%>%
  mutate(area_m2_scaled = scale(area_m2)[,1],
         tss_med_scaled = scale(tss_med)[,1],
         age_mean_scaled = scale(age_mean)[,1],
         ws_cpy_scaled = scale(ws_cpy)[,1],
         str_cpy_scaled = scale(str_cpy)[,1],
         str_dens_scaled = scale (str_dens)[,1],
         pop_dens_scaled = scale (pop_dens)[,1],
         aadt_scaled = scale (aadt)[,1],
         res_scaled = scale (res)[,1])%>%
  filter(aadt != "NA",
         ws_cpy != "NA",
         str_cpy != "NA",
         tss_med != "NA")
  
#runoff Pb ~ avg parcel age 

m_1 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_1 <- add_criterion(m_1, "loo")

#compare to a model that includes size of watershed

m_2 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_2 <- add_criterion(m_2, "loo")

loo(m_1, m_2)

#runoff Pb ~ avg parcel age + median TSS

m_3 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_3 <- add_criterion(m_3, "loo")

loo_compare(m_1, m_2, m_3)

#runoff Pb ~ avg parcel age + median TSS + pop density 

m_4 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + pop_dens_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_4 <- add_criterion(m_4, "loo")

loo_compare(m_1, m_2, m_3, m_4)

#runoff Pb ~ avg parcel age + median TSS + street density

m_9 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + str_dens_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_9 <- add_criterion(m_9, "loo")

loo_compare(m_1, m_2, m_3, m_4, m_9)

#runoff Pb ~ avg parcel age + median TSS + traffic volume 

m_5 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + aadt_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_5 <- add_criterion(m_5, "loo")

loo_compare(m_1, m_2, m_3, m_4, m_5)

#runoff Pb ~ avg parcel age + median TSS + tree canopy

m_6 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + ws_cpy_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_6 <- add_criterion(m_6, "loo")

loo_compare(m_1, m_2, m_3, m_4, m_5, m_6)

m_7 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + str_cpy_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_7 <- add_criterion(m_7, "loo")

loo_compare(m_1, m_2, m_3, m_4, m_5, m_6, m_7)

#runoff Pb ~ avg parcel age + median TSS + proportion land residential

m_8 <- brm(data = sites_2, family = lognormal(),
           pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + res_scaled + (1|site),
           iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

m_8 <- add_criterion(m_8, "loo")

loo_compare(m_1, m_2, m_3, m_4, m_5, m_6, m_7, m_8, m_9)


####exploration of best fit model structure and plot####

sites_1 <- sites_ws%>%
  mutate(area_m2_scaled = scale(area_m2)[,1],
         tss_med_scaled = scale(tss_med)[,1],
         age_mean_scaled = scale(age_mean)[,1],
         ws_cpy_scaled = scale(ws_cpy)[,1],
         str_cpy_scaled = scale(str_cpy)[,1],
         str_dens_scaled = scale (str_dens)[,1],
         pop_dens_scaled = scale (pop_dens)[,1],
         aadt_scaled = scale (aadt)[,1],
         res_scaled = scale (res)[,1])

m_ind_bf_log <- brm(data = sites_1, family = lognormal(),
                    pb_ppb ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + pop_dens_scaled + (1|site),
                    iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())

summary(m_ind_bf_log)


#condition effects of parcel age

p <- conditional_effects(m_ind_bf_log) %>% plot(print = FALSE)

p$age_mean_scaled$data %>%
  unscale_dataframe(original = sites_1) %>%
  ggplot(aes(x = age_mean)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = sites_1, aes(y = pb_mean))+
  xlab("Mean parcel age within watershed (year)")+
  ylab("Mean runoff Pb within watershed (ppb)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#conditional effects of tss

p$tss_med_scaled$data %>%
  unscale_dataframe(original = sites_1) %>%
  ggplot(aes(x = tss_med)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = sites_1, aes(y = pb_mean))+
  xlab("Median TSS (mg/L)")+
  ylab("Mean runoff Pb within watershed (ppb)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  
#conditional effects of area
p$area_m2_scaled$data %>%
  unscale_dataframe(original = sites_1) %>%
  ggplot(aes(x = area_m2)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = sites_1, aes(y = pb_mean)) +
  xlab("Watershed area (m^2)") +
  ylab("Mean runoff Pb within watershed (ppb)") +
  theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())

#pop dens
p$pop_dens_scaled$data %>%
  unscale_dataframe(original = sites_1) %>%
  ggplot(aes(x = pop_dens)) +
  geom_line(aes(y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data = ws, aes(y = pb_mean))+
  xlab("Population density (people/km^2")+
  ylab("Mean runoff Pb within watershed (ppb)")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

summary(m_ws_bf_log)