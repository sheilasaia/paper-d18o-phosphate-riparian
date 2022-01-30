# ---- script header ----
# script name: analysis_script.R
# purpose of script: statistical analysis and plotting code used in paper analysis
# author: sheila saia
# date created: 2021-03-11
# email: ssaia@ncsu.edu


# ---- notes ----
# notes:

# helpful lme tutorial: https://bodowinter.com/tutorial/bw_LME_tutorial2.pdf
# plotting lme results: https://lmudge13.github.io/sample_code/mixed_effects.html
# help on effects: https://stats.stackexchange.com/questions/233007/interpreting-effects-plots-in-r/233018


# ---- to do ----
# to do list


# ---- 1. load libraries ----
# library(devtools)
# devtools::install_github("strengejacke/sjPlot")

library(tidyverse)
library(forcats)
library(here)
library(car) # for qqPlot
library(lme4)
library(sjPlot)
library(effects)
library(renv)


# ---- 2. load data ----
# load d18O data
d18o_data <- read_csv(here::here("data", "d18o_data.csv"), col_names = TRUE)
# these data are from the pooled (Pi or Po???) fractions for each extraction

# load Hedley P fraction data
hedley_data <- read_csv(here::here("data", "hedley_data.csv"), col_names = TRUE) %>%
  mutate(extract = fct_relevel(extract, c("H2O", "NaHCO3", "NaOH", "HCl")))
# inorganic p (Pi) has triplicates (three reps)
# organic p (Po) has only one rep


# ---- 3. d18o-phos modeling and figure 3 ----
# check normality in response
hist(d18o_data$d18o_value) # hard to tell --> look at qqplot
qqPlot(d18o_data$d18o_value) # data is along qqline and within bounds --> normally distributed

# model without depth
# include id as random effect
# don't have enough samples for random slope models (i.e., 1 + transect_dist_m | bonn_id)
# use random intercept models only (i.e., 1 | bonn_id)
d18o_lme0 <- lmer(d18o_value ~ transect_dist_m + (1 | bonn_id), data = d18o_data, REML = FALSE)
summary(d18o_lme0)
sjPlot::tab_model(d18o_lme0)

# model with only depth (no distance)
# include id as random effect
d18o_lme1 <- lmer(d18o_value ~ depth_cm + (1 | bonn_id), data = d18o_data, REML = FALSE)
summary(d18o_lme1)
sjPlot::tab_model(d18o_lme1)

# model with depth and distance
# include id as random effect
d18o_lme2 <- lmer(d18o_value ~ transect_dist_m + depth_cm + (1 | bonn_id), data = d18o_data, REML = FALSE)
summary(d18o_lme2)
sjPlot::tab_model(d18o_lme2)

# model with depth and distance interaction
# include id as random effect
d18o_lme3 <- lmer(d18o_value ~ transect_dist_m * depth_cm + (1 | bonn_id), data = d18o_data, REML = FALSE)
summary(d18o_lme3)
sjPlot::tab_model(d18o_lme3)

# AIC to compare models
AIC(d18o_lme1, d18o_lme0, d18o_lme2, d18o_lme3)
# model 2 and 3 are best but not different so go with the simplest one -> lme2

# anova to look at p value of having depth and distance
# anova(d18o_lme1, d18o_lme2)
# distance impacted d18o-phosphate (X2(1) = 31.38, p = 2.12e-8), 
# lowering it by 0.15 +/- 0.02 per mil

# look at residuals of first model
# hist of residuals
hist(resid(d18o_lme2)) # looks normally distributed
qqPlot(resid(d18o_lme2)) # looks normally distributed
# check for constant variance
plot(resid(d18o_lme2) ~ fitted(d18o_lme2))
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(resid(d18o_lme2) ~ d18o_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)

# look at effects
# sjPlot::plot_model(d18o_lme2)

# get fixed effects (i.e., distance)
lme2_fixed_effects <- effects::effect(term = "transect_dist_m", mod = d18o_lme2)
summary(lme2_fixed_effects)

# save as df
lme2_fixed_effects_df <- as.data.frame(lme2_fixed_effects)

# plot
# plot data and model (figure 3)
# pdf(here::here("figures_tables", "figure_3.pdf"), width = 12, height = 10)
ggplot() +
  # observations
  geom_point(data = d18o_data,
             mapping = aes(x = transect_dist_m, y = d18o_value),
             size  = 2, alpha = 0.75) +
  geom_rect(mapping = aes(xmin = 0, xmax = 60, ymin = 16, ymax = 20), 
            alpha = 0.2, fill = "grey25") +
  # model
  geom_line(data = lme2_fixed_effects_df,
            mapping = aes(x = transect_dist_m, y = fit),
            color = "#66c2a5", lwd = 1, lty = 1) + 
  geom_ribbon(data = lme2_fixed_effects_df,
              mapping = aes(x = transect_dist_m, y = fit, 
                            ymin = lower, ymax = upper),
              fill = "#66c2a5", alpha = 0.2) +
  xlim(0, 60) +
  labs(x = "Transect Distance (m)", 
       y = expression(paste(delta^{18}*O[P], " (\u2030)"))) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        text = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
# dev.off()
# getting error if use pdf() function so export manually


# ---- 4. separate Pi and Po data ----
# transect distance key
transect_dist_key <- data.frame(transect_dist_id = seq(1, 6, 1),
                                transect_dist_m = unique(d18o_data$transect_dist_m))

# filter Pi data
inorg_hedley_data <- hedley_data %>%
  filter(fraction == "Pi") %>%
  left_join(transect_dist_key, by = "transect_dist_id") %>%
  filter(soil_moist_type == "dried")
# note: there are triplicate measuremens for pi

# filter Po data
org_hedley_data <- hedley_data %>%
  filter(fraction == "Po")  %>%
  left_join(transect_dist_key, by = "transect_dist_id")  %>%
  filter(soil_moist_type == "dried")
# note: po values calcuated as total p minus mean of pi for that extraction and depth
# total p was only measured once per sample and raw data were lost due to laboratory closure


# ---- 5. hedley data plots (figures 1 and 2) ----
# colors for plots
#my_colors <- c("white", "grey60", "grey40", "grey30")

# inorganic P fractions (figure 1)
pdf(here::here("figures_tables", "figure_1.pdf"), width = 10, height = 10, pointsize = 14)
ggplot(data = inorg_hedley_data,
       mapping = aes(x = factor(transect_dist_m), y = p_conc_mgperkg, fill = extract)) +
  geom_boxplot() +
  geom_point(alpha = 0.3, size = 2) +
  facet_wrap(~ extract, scales = "free") +
  labs(x = "Transect Distance (m)", 
       y = "Inorganic P Concentration (mg per kg)",
       fill = "Extract") +
  #scale_fill_manual(values = my_colors) +
  theme_classic()
dev.off()

# organic P fractions (figure 2)
pdf(here::here("figures_tables", "figure_2.pdf"), width = 10, height = 10, pointsize = 14)
ggplot(data = org_hedley_data,
       mapping = aes(x = factor(transect_dist_m), y = p_conc_mgperkg, fill = extract)) +
  geom_boxplot() +
  geom_point(alpha = 0.3, size = 2) +
  facet_wrap(~ extract, scales = "free") +
  labs(x = "Transect Distance (m)", 
       y = "Organic P Concentration (mg per kg)",
       fill = "Extract") +
  #scale_fill_manual(values = my_colors) +
  theme_classic()
dev.off()


# ---- 5. Pi modeling ----
# h2o extraction
# select data
h2o_pi_data <- inorg_hedley_data %>%
  filter(extract == "H2O")

# check normality in response
hist(log(h2o_pi_data$p_conc_mgperkg)) # looks normally distributed
# qqnorm(log(h2o_pi_data$p_conc_mgperkg))
# qqline(log(h2o_pi_data$p_conc_mgperkg))
qqPlot(log(h2o_pi_data$p_conc_mgperkg))
# most of data is within bounds when logged

# model without depth
h2o_pi_lm1 <- lm(log(p_conc_mgperkg) ~ transect_dist_m, data = h2o_pi_data)
summary(h2o_pi_lm1)

# model with depth
h2o_pi_lm2 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm, data = h2o_pi_data)
summary(h2o_pi_lm2)

# model with depth and interaction
h2o_pi_lm3 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = h2o_pi_data)
summary(h2o_pi_lm3)

# AIC comparison
AIC(h2o_pi_lm1, h2o_pi_lm2, h2o_pi_lm3)
# models 2 and 3 are units so take lm2

# hist of residuals
hist(h2o_pi_lm2$residuals) # looks normally distributed
# qqnorm(h2o_pi_lm2$residuals)
# qqline(h2o_pi_lm2$residuals)
qqPlot(h2o_pi_lm2$residuals) # looks normally distributed
# check for constant variance
plot(h2o_pi_lm2$residuals ~ h2o_pi_lm2$fitted.values)
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(h2o_pi_lm2$residuals ~ h2o_pi_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# nahco3 extraction
# select data
nahco3_pi_data <- inorg_hedley_data %>%
  filter(extract == "NaHCO3")

# check normality in response
hist(log(nahco3_pi_data$p_conc_mgperkg)) # more normally distributed after logging
# qqnorm(nahco3_pi_data$p_conc_mgperkg)
# qqline(nahco3_pi_data$p_conc_mgperkg)
qqPlot(log(nahco3_pi_data$p_conc_mgperkg))
# data is within bounds when logged

# model without depth
nahco3_pi_lm1 <- lm(log(p_conc_mgperkg) ~ transect_dist_m, data = nahco3_pi_data)
summary(nahco3_pi_lm1)

# model with depth
nahco3_pi_lm2 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm, data = nahco3_pi_data)
summary(nahco3_pi_lm2)

# model with depth and interaction
nahco3_pi_lm3 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = nahco3_pi_data)
summary(nahco3_pi_lm3)

# AIC comparison
AIC(nahco3_pi_lm1, nahco3_pi_lm2, nahco3_pi_lm3)
# model 3 is lowest

# hist of residuals
hist(nahco3_pi_lm3$residuals) # not sure if this looks normal
# qqnorm(nahco3_pi_lm3$residuals)
# qqline(nahco3_pi_lm3$residuals)
qqPlot(nahco3_pi_lm3$residuals) # not sure if this looks normal
# check for constant variance
plot(nahco3_pi_lm3$residuals ~ nahco3_pi_lm3$fitted.values)
# potential funnel pattern in residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(nahco3_pi_lm3$residuals ~ nahco3_pi_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# naoh extraction
# select data
naoh_pi_data <- inorg_hedley_data %>%
  filter(extract == "NaOH")

# check normality in response
hist(log(naoh_pi_data$p_conc_mgperkg)) # more normally distributed after logging
# qqnorm(log(naoh_pi_data$p_conc_mgperkg))
# qqline(log(naoh_pi_data$p_conc_mgperkg))
qqPlot(log(naoh_pi_data$p_conc_mgperkg))
# looks normally distributed when logged

# model without depth
naoh_pi_lm1 <- lm(log(p_conc_mgperkg) ~ transect_dist_m, data = naoh_pi_data)
summary(naoh_pi_lm1)

# model with depth
naoh_pi_lm2 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm, data = naoh_pi_data)
summary(naoh_pi_lm2)

# model with depth and interaction
naoh_pi_lm3 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = naoh_pi_data)
summary(naoh_pi_lm3)

# AIC comparison
AIC(naoh_pi_lm1, naoh_pi_lm2, naoh_pi_lm3)
# lm3 is lowest

# hist of residuals
hist(naoh_pi_lm3$residuals) # not sure if this looks normal
# qqnorm(naoh_pi_lm3$residuals)
# qqline(naoh_pi_lm3$residuals)
qqPlot(naoh_pi_lm3$residuals) # not sure if this looks normal
# check for constant variance
plot(naoh_pi_lm3$residuals ~ naoh_pi_lm3$fitted.values)
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(naoh_pi_lm3$residuals ~ naoh_pi_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# hcl extraction
# select data
hcl_pi_data <- inorg_hedley_data %>%
  filter(extract == "HCl")

# check normality in response
hist(log(hcl_pi_data$p_conc_mgperkg)) # not sure about being normally distributed
# qqnorm(log(hcl_pi_data$p_conc_mgperkg))
# qqline(log(hcl_pi_data$p_conc_mgperkg))
qqPlot(log(hcl_pi_data$p_conc_mgperkg))
# middle of curve looks better when logged

# model without depth
hcl_pi_lm1 <- lm(log(p_conc_mgperkg) ~ transect_dist_m, data = hcl_pi_data)
summary(hcl_pi_lm1)

# model with depth
hcl_pi_lm2 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm, data = hcl_pi_data)
summary(hcl_pi_lm2)

# model with depth and interaction
hcl_pi_lm3 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = hcl_pi_data)
summary(hcl_pi_lm3)

# AIC comparison
AIC(hcl_pi_lm1, hcl_pi_lm2, hcl_pi_lm3)
# lm3 is lowest

# hist of residuals
hist(hcl_pi_lm3$residuals) # not sure if this looks normal
# qqnorm(hcl_pi_lm3$residuals)
# qqline(hcl_pi_lm3$residuals)
qqPlot(hcl_pi_lm3$residuals) # most of data in bounds except at tails
# check for constant variance
plot(hcl_pi_lm3$residuals ~ hcl_pi_lm3$fitted.values)
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(hcl_pi_lm3$residuals ~ hcl_pi_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# ---- 6. Po modeling ----
# h2o extraction
# select data
h2o_po_data <- org_hedley_data %>%
  filter(extract == "H2O")

# check normality in response
hist(log(h2o_po_data$p_conc_mgperkg)) # more normally distributed after logging
# qqnorm(log(h2o_po_data$p_conc_mgperkg))
# qqline(log(h2o_po_data$p_conc_mgperkg))
qqPlot(log(h2o_po_data$p_conc_mgperkg))
# looks normally distributed after logging

# model without depth
h2o_po_lm1 <- lm(log(p_conc_mgperkg) ~ transect_dist_m, data = h2o_po_data)
summary(h2o_po_lm1)

# model with depth
h2o_po_lm2 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm, data = h2o_po_data)
summary(h2o_po_lm2)

# model with depth and interaction
h2o_po_lm3 <- lm(log(p_conc_mgperkg) ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = h2o_po_data)
summary(h2o_po_lm3)

# AIC comparison
AIC(h2o_po_lm1, h2o_po_lm2, h2o_po_lm3)
# lm3 has the lowest value

# hist of residuals
hist(h2o_po_lm3$residuals) # looks normally distributed
# qqnorm(h2o_po_lm3$residuals)
# qqline(h2o_po_lm3$residuals)
qqPlot(h2o_po_lm3$residuals) # looks normally distributed
# check for constant variance
plot(h2o_po_lm3$residuals ~ h2o_po_lm1$fitted.values)
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(h2o_po_lm3$residuals ~ h2o_po_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# nahco3 extraction
# select data
nahco3_po_data <- org_hedley_data %>%
  filter(extract == "NaHCO3")

# check normality in response
hist(nahco3_po_data$p_conc_mgperkg) # looks normally distributed
# qqnorm(nahco3_po_data$p_conc_mgperkg)
# qqline(nahco3_po_data$p_conc_mgperkg)
qqPlot(nahco3_po_data$p_conc_mgperkg)
# looks normally distributed without logging

# model without depth
nahco3_po_lm1 <- lm(p_conc_mgperkg ~ transect_dist_m, data = nahco3_po_data)
summary(nahco3_po_lm1)

# model with depth
nahco3_po_lm2 <- lm(p_conc_mgperkg ~ transect_dist_m + depth_cm, data = nahco3_po_data)
summary(nahco3_po_lm2)

# model with depth and interaction
nahco3_po_lm3 <- lm(p_conc_mgperkg ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = nahco3_po_data)
summary(nahco3_po_lm3)

# AIC comparison
AIC(nahco3_po_lm1, nahco3_po_lm2, nahco3_po_lm3)
# lm2 and lm3 are not sign. different so go with lm2

# hist of residuals
hist(nahco3_po_lm2$residuals) # looks normally distributed
# qqnorm(nahco3_po_lm2$residuals)
# qqline(nahco3_po_lm2$residuals)
qqPlot(nahco3_po_lm2$residuals) # looks normally distributed
# check for constant variance
plot(nahco3_po_lm2$residuals ~ nahco3_po_lm1$fitted.values)
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(nahco3_po_lm2$residuals ~ nahco3_po_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# naoh extraction
# select data
naoh_po_data <- org_hedley_data %>%
  filter(extract == "NaOH")

# check normality in response
hist(naoh_po_data$p_conc_mgperkg) # looks normally distributed
# qqnorm(naoh_po_data$p_conc_mgperkg))
# qqline(naoh_po_data$p_conc_mgperkg))
qqPlot(naoh_po_data$p_conc_mgperkg)
# looks normally distributed without logging

# model without depth
naoh_po_lm1 <- lm(p_conc_mgperkg ~ transect_dist_m, data = naoh_po_data)
summary(naoh_po_lm1)

# model with depth
naoh_po_lm2 <- lm(p_conc_mgperkg ~ transect_dist_m + depth_cm, data = naoh_po_data)
summary(naoh_po_lm2)

# model with depth and interaction
naoh_po_lm3 <- lm(p_conc_mgperkg ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = naoh_po_data)
summary(naoh_po_lm3)

# AIC comparison
AIC(naoh_po_lm1, naoh_po_lm2, naoh_po_lm3)
# lm3 is lowest

# hist of residuals
hist(naoh_po_lm3$residuals) # looks normally distributed
# qqnorm(naoh_po_lm3$residuals)
# qqline(naoh_po_lm3$residuals)
qqPlot(naoh_po_lm3$residuals) # looks normally distributed
# check for constant variance
plot(naoh_po_lm3$residuals ~ naoh_po_lm3$fitted.values)
# funnel-like
plot(naoh_po_lm3$residuals ~ naoh_po_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# hcl extraction
# select data
hcl_po_data <- org_hedley_data %>%
  filter(extract == "HCl")

# check normality in response
hist(hcl_po_data$p_conc_mgperkg) # looks normally distrubted
# qqnorm(hcl_po_data$p_conc_mgperkg)
# qqline(hcl_po_data$p_conc_mgperkg)
qqPlot(hcl_po_data$p_conc_mgperkg) # looks normally distrubted

# model without depth
hcl_po_lm1 <- lm(p_conc_mgperkg ~ transect_dist_m, data = hcl_po_data)
summary(hcl_po_lm1)

# model with depth
hcl_po_lm2 <- lm(p_conc_mgperkg ~ transect_dist_m + depth_cm, data = hcl_po_data)
summary(hcl_po_lm2)

# model with depth and interaction
hcl_po_lm3 <- lm(p_conc_mgperkg ~ transect_dist_m + depth_cm + transect_dist_m:depth_cm, data = hcl_po_data)
summary(hcl_po_lm3)

# AIC comparison
AIC(hcl_po_lm1, hcl_po_lm2, hcl_po_lm3)
# lm2 is lowest

# hist of residuals
hist(hcl_po_lm2$residuals) # looks normally distributed
# qqnorm(hcl_po_lm2$residuals)
# qqline(hcl_po_lm2$residuals)
qqPlot(hcl_po_lm2$residuals) # looks normally distributed
# check for constant variance
plot(hcl_po_lm2$residuals ~ hcl_po_lm2$fitted.values)
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(hcl_po_lm2$residuals ~ hcl_po_data$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)


# ---- 7. redox modeling and figure s3 ----
# add redox labeling and drop furthest distance
d18o_data_redox <- d18o_data %>%
  filter(transect_dist_id <= 4) %>%
  mutate(redox_status = case_when(depth_cm == 10 ~ "oxidized",
                                  depth_cm == 30 ~ "other",
                                  depth_cm == 50 ~ "other",
                                  depth_cm == 100 ~ "reduced")) %>%
  filter(redox_status != "other")

# check normality in response
hist(d18o_data_redox$d18o_value) # hard to tell --> look at qqplot
qqPlot(d18o_data_redox$d18o_value) # data is along qqline and within bounds --> normally distributed

# model without distance
# include id as random effect
# don't have enough samples for random slope models (i.e., 1 + transect_dist_m | bonn_id)
# use random intercept models only (i.e., 1 | bonn_id)
d18o_lme4 <- lmer(d18o_value ~ redox_status + (1 | bonn_id), data = d18o_data_redox, REML=FALSE)
summary(d18o_lme4)
sjPlot::tab_model(d18o_lme4)

# model with distance and redox status
# include id as random effect
d18o_lme5 <- lmer(d18o_value ~ transect_dist_m + redox_status + (1 | bonn_id), data = d18o_data_redox, REML=FALSE)
summary(d18o_lme5)
sjPlot::tab_model(d18o_lme5)

# model with distance and redox status as well as interaction
# include id as random effect
d18o_lme6 <- lmer(d18o_value ~ transect_dist_m * redox_status + (1 | bonn_id), data = d18o_data_redox, REML=FALSE)
summary(d18o_lme6)
sjPlot::tab_model(d18o_lme6)

# AIC comparison
AIC(d18o_lme4, d18o_lme5, d18o_lme6)
# model 5 and 6 are best but not different so take the simplest one -> lme5

# look at residuals of first model
# hist of residuals
hist(resid(d18o_lme5)) # looks normally distributed
qqPlot(resid(d18o_lme5)) # looks normally distributed
# check for constant variance
plot(resid(d18o_lme5) ~ fitted(d18o_lme5))
# no linear or funnel pattern in the residuals vs fitted values --> homoscedastic (i.e., constant variance)
plot(resid(d18o_lme5) ~ d18o_data_redox$transect_dist_m)
# no linear or funnel pattern in the residuals vs distance --> homoscedastic (i.e., constant variance)

# look at effects
# sjPlot::plot_model(d18o_lme5)

# get fixed effects (i.e., distance)
lme5_fixed_effects <- effects::Effect(focal.predictors = c("transect_dist_m", "redox_status"), mod = d18o_lme5, 
                                      xlevels = list(transect_dist_m = c(1, 10, 20, 30, 40, 45), redox_status = c("oxidized", "reduced")))
summary(lme5_fixed_effects)

# save as df
lme5_fixed_effects_df <- as.data.frame(lme5_fixed_effects)

# pick colors
my_redox_colors <- c("#fc8d62", "#8da0cb")

# plot data and model
# pdf(here::here("figures_tables", "figure_s3.pdf"), width = 12, height = 10)
ggplot() +
  # observations
  geom_point(data = d18o_data_redox,
             mapping = aes(x = transect_dist_m, y = d18o_value, shape = redox_status),
             size  = 3, alpha = 0.75) +
  geom_rect(mapping = aes(xmin = 0, xmax = 45, ymin = 16, ymax = 20), 
            alpha = 0.2, fill = "grey25") +
  # oxidized model
  geom_line(data = lme5_fixed_effects_df %>% filter(redox_status == "oxidized"),
            mapping = aes(x = transect_dist_m, y = fit),
            color = my_redox_colors[1], lwd = 1, lty = 2) + 
  geom_ribbon(data = lme5_fixed_effects_df %>% filter(redox_status == "oxidized"),
              mapping = aes(x = transect_dist_m, y = fit, 
                            ymin = lower, ymax = upper),
              fill = my_redox_colors[1], alpha = 0.3) +
  # reduced model
  geom_line(data = lme5_fixed_effects_df %>% filter(redox_status == "reduced"),
            mapping = aes(x = transect_dist_m, y = fit),
            color = my_redox_colors[2], lwd = 1, lty = 1) + 
  geom_ribbon(data = lme5_fixed_effects_df %>% filter(redox_status == "reduced"),
              mapping = aes(x = transect_dist_m, y = fit, 
                            ymin = lower, ymax = upper),
              fill = my_redox_colors[2], alpha = 0.3) +
  xlim(0, 45) +
  labs(x = "Transect Distance (m)", 
       y = expression(paste(delta^{18}*O[P], " (\u2030)"))) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        text = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
# dev.off()
# getting error if use pdf() function so export manually

