#analysis_comparison for other metals compared to pb distribution

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(readr)

#####loading and cleaning####

ws <- read_csv("~/Watershed_aggregagted_attributes.csv")%>%
  select("site", "pb_mean", "area_m2", "tss_med", "age_mode", "age_mean", "ws_cpy", "str_cpy", "str_dens", "pop_dens", "aadt", "res")

##need to get the dissagreggated data for Zn, Cu, and Ni

sites <- read_csv("~/stormwater_metals_database.csv")%>%
  select("site", "month", "year", "pb", "cu", "ni", "zn") %>%
  filter(site != "mprb-9")%>%
  filter(site != "crwd-40")

sites_ws <- left_join(sites, ws, by = "site")

#scale variables

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


####cu####
m_cu_bfpb <- brm(data = sites_1, family = lognormal(),
                    cu ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + pop_dens_scaled + (1|site),
                    iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())
summary(m_cu_bfpb)

print(summary(m_cu_bfpb), digits = 4) 

####ni####
m_ni_bfpb <- brm(data = sites_1, family = lognormal(),
                 ni ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + pop_dens_scaled + (1|site),
                 iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())
summary(m_ni_bfpb)

print(summary(m_ni_bfpb), digits = 4) 
####zn####
m_zn_bfpb <- brm(data = sites_1, family = lognormal(),
                 zn ~ age_mean_scaled + area_m2_scaled + tss_med_scaled + pop_dens_scaled + (1|site),
                 iter = 10000, warmup = 3000, chains = 3, cores = future::availableCores())
summary(m_zn_bfpb)

print(summary(m_zn_bfpb), digits = 4) 

