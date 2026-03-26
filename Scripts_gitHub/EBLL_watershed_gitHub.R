#analysis for predicting EBLL at the watershed level

library(tidyr)
library(dplyr)
library(ggplot2)
library(brms)
library(readr)


####loading and cleaning####

ws <- read_csv("~/Watershed_aggregagted_attributes.csv")%>%
  select("site", "pb_mean", "area_m2", "res", "EBLL_weighted_smallerThan2_methodology", "hh_incm_lessThan2_2013to2017", "pop_dens", "aadt", "age_mean", "str_dens") %>%
  filter(site != "mprb-9")%>%
  filter(site != "crwd-40")%>%
  filter(site != "swwd-8")

####model comparison####

#scale variables and filter

ws_2 <- ws %>%
  mutate(area_m2_scaled = scale(area_m2)[,1],
         pb_mean_scaled = scale(pb_mean)[,1],
         age_mean_scaled = scale(age_mean)[,1],
         hh_incm_lessThan2_2013to2017_scaled = scale(hh_incm_lessThan2_2013to2017)[,1],
         str_dens_scaled = scale (str_dens)[,1],
         pop_dens_scaled = scale (pop_dens)[,1],
         aadt_scaled = scale (aadt)[,1],
         res_scaled = scale (res)[,1])%>%
  filter(aadt != "NA",
         pb_mean != "NA")

#%EBLL ~ runoff Pb

m_1 <- brm(data = ws_2, family = "gamma",
              EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled,
              iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_1 <- add_criterion(m_1, "loo")

#%EBLL ~ parcel age

m_2 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ age_mean_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_2 <- add_criterion(m_2, "loo")

#%EBLL ~ avg household income + mean runoff Pb

m_4 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled + hh_incm_lessThan2_2013to2017_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_4 <- add_criterion(m_4, "loo")
loo_compare(m_1, m_2, m_4)

#%EBLL ~ avg household income * mean runoff Pb

m_5 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled * hh_incm_lessThan2_2013to2017_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_5 <- add_criterion(m_5, "loo")
loo_compare(m_1, m_2, m_4, m_5)

#%EBLL ~ pop density + mean runoff Pb 

m_6 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled + pop_dens_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_6 <- add_criterion(m_6, "loo")
loo_compare(m_1, m_2, m_4, m_5, m_6)

#%EBLL ~ str density + mean runoff Pb

m_7 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled + str_dens_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_7 <- add_criterion(m_7, "loo")
loo_compare(m_1, m_2, m_4, m_5, m_6, m_7 )

#%EBLL ~ traffic volume + mean runoff Pb

m_8 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled + aadt_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_8 <- add_criterion(m_8, "loo")
loo_compare(m_1, m_2, m_4, m_5, m_6, m_7, m_8 )

#%EBLL ~ proportion of residential area + mean runoff Pb

m_9 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled + res_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_9 <- add_criterion(m_9, "loo")
loo_compare(m_1, m_2, m_4, m_5, m_6, m_7, m_8, m_9 )

#%EBLL ~ avg household income * mean runoff Pb + pop density

m_10 <- brm(data = ws_2, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled * hh_incm_lessThan2_2013to2017_scaled + pop_dens_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_10 <- add_criterion(m_10, "loo")
loo_compare(m_1, m_2, m_4, m_5, m_6, m_7, m_8, m_9, m_10)


####exploration of best fit model and plot####

ws_3 <- ws %>%
  mutate(area_m2_scaled = scale(area_m2)[,1],
         pb_mean_scaled = scale(pb_mean)[,1],
         age_mean_scaled = scale(age_mean)[,1],
         hh_incm_lessThan2_2013to2017_scaled = scale(hh_incm_lessThan2_2013to2017)[,1],
         str_dens_scaled = scale (str_dens)[,1],
         pop_dens_scaled = scale (pop_dens)[,1],
         aadt_scaled = scale (aadt)[,1],
         res_scaled = scale (res)[,1])

m_bestfit<- brm(data = ws_3, family = "gamma",
               EBLL_weighted_smallerThan2_methodology ~  pb_mean_scaled * hh_incm_lessThan2_2013to2017_scaled,
               iter = 10000, warmup = 3000, chains = 4, cores = future::availableCores(), control = list(adapt_delta = 0.999, max_treedepth = 17))

summary(m_bestfit)

#condition effects of hh income

p <- conditional_effects(m_bestfit, 
                         int_conditions = list(hh_incm_lessThan2_2013to2017_scaled = c(-1, 1))) %>% plot(print = FALSE)


ws_3$hh_incm_lessThan2_2013to2017

ws_3_plot <- ws_3 %>%
  mutate(income_group = ifelse(hh_incm_lessThan2_2013to2017_scaled >= 0, 
                               "Above mean", 
                               "Below mean"))

p$`pb_mean_scaled:hh_incm_lessThan2_2013to2017_scaled`$data %>%
  unscale_dataframe(original = ws_3_plot) %>%
  ggplot(aes(x = pb_mean)) +
  geom_line(aes(y = estimate__,
                color = factor(hh_incm_lessThan2_2013to2017_scaled),
                fill = factor(hh_incm_lessThan2_2013to2017_scaled))) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__,
                  fill = factor(hh_incm_lessThan2_2013to2017_scaled)), alpha = 0.2) +
  geom_point(data = ws_3_plot, aes(y = EBLL_weighted_smallerThan2_methodology,  shape = income_group))+
  xlab("Average Lead in Watershed Runoff (ppb)")+
  ylab("Elevated Blood Lead Level (percent elevated)")+
  labs(color = "Model predictions household income\n(standard deviations from mean)",
       fill = "Model predictions household income\n(standard deviations from mean)",
       shape = "Raw data for household income") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_cartesian(ylim = c(0, 10))

####subset soil analysis####

#additional model exploration of adding soil pb as a predictor of EBLL versus runoff Pb or parcel age

#load pb values

pb_soil <- read_csv("~/stormshedsXsoilPb_FIXED.csv")%>%
  select("site", "meanPb_20s", "n_size_20s") %>%
  filter(site != "mprb-9")%>%
  filter(site != "crwd-40")%>%
  filter(site != "swwd-8")%>%
  filter(n_size_20s > 5)

#add to ws#
pb_ws <- left_join(pb_soil, ws, "site") %>%
  filter(pb_mean != "NA")

#only 11 stormsheds where we have an estimate of pb in the soil and pb in the runoff
ws_4 <- pb_ws %>%
  mutate(area_m2_scaled = scale(area_m2)[,1],
         pb_mean_scaled = scale(pb_mean)[,1],
         age_mean_scaled = scale(age_mean)[,1],
         hh_incm_lessThan2_2013to2017_scaled = scale(hh_incm_lessThan2_2013to2017)[,1],
         str_dens_scaled = scale (str_dens)[,1],
         pop_dens_scaled = scale (pop_dens)[,1],
         aadt_scaled = scale (aadt)[,1],
         res_scaled = scale (res)[,1],
         meanPb_20s_scaled = scale (meanPb_20s)[,1])

#%EBLL ~ runoff Pb

m_1 <- brm(data = ws_4, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ pb_mean_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_1 <- add_criterion(m_1, "loo", moment_match = TRUE)

#%EBLL ~ parcel age

m_2 <- brm(data = ws_4, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ age_mean_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_2 <- add_criterion(m_2, "loo", moment_match = TRUE)

#%EBLL ~ soil pb

m_3 <- brm(data = ws_4, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ meanPb_20s_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_3 <- add_criterion(m_3, "loo", moment_match = TRUE)

#%EBLL ~ avg household income

m_4 <- brm(data = ws_4, family = "gamma",
           EBLL_weighted_smallerThan2_methodology ~ hh_incm_lessThan2_2013to2017_scaled,
           iter = 1000, warmup = 300, chains = 3, cores = future::availableCores())

m_4 <- add_criterion(m_4, "loo")


loo_compare(m_1, m_2, m_3, m_4)

summary(m_4)
summary(m_1)


