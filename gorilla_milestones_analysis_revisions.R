# Gorilla milestones analysis - revisions ####

# 1 Description ####

# Based on feedback in reviews, we need to run separate models for the 10
# behaviours where we have a decent sample size in each species.

# We also need to illustrate complexity better

# Dataset:
# chimp_gorilla: milestone data extracted from our study and that of Brundl for plotting
# fig1_d: means and standard deviations of all our milestones in each species
# m1d: milestone data organised for between species comparisons - 
# this includes only milestones present in both species; we will subset this data 
# so each milestone is modelled individually
# m2d: just the bwindi data for sex comparisons.

# 2 Packages ####

library(tidyverse)
library(beepr)
library(GGally)
library(brms)
library(cmdstanr)
library(car)
library(reshape2)
library(MetBrewer)
library(ggforce)
library(patchwork)

# 3 Add in complexity variable ####

milestone <- sort(unique(m2d$milestone))
milestone <- as.data.frame(milestone)

complexity <- c(
  2,
  1,
  1,
  1,
  2,
  1,
  1,
  1,
  1,
  1,
  2,
  1,
  2,
  2,
  0,
  0,
  1,
  2,
  2,
  2,
  1,
  0,
  0,
  0,
  0,
  0,
  0,
  0)
  

milestone$complexity <- complexity

test <- m1d
test <- test %>%
  left_join(milestone, by = "milestone")

m1d <- test

test <- m2d
test <- test %>%
  left_join(milestone, by = "milestone")

m2d <- test

# figure data worded more cleanly so needs its own adding in

milestone2 <- sort(unique(fig1_d$milestone))
milestone2 <- as.data.frame(milestone2)

complexity <- c(
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  1,
  1,
  1,
  1,
  1,
  1,
  1,
  1,
  2,
  1,
  2,
  2,
  2,
  0,
  0,
  1,
  2,
  2,
  2,
  2,
  1)

milestone2$complexity <- complexity
milestone2$milestone <- milestone2$milestone2
milestone2$milestone2 <- NULL

test <- fig1_d
test <- test %>%
  left_join(milestone2, by = "milestone")

fig1_d <- test

rm(milestone, milestone2,test,complexity)

# 3 Reviewers wanted sample size per sex for the milestones ####

xx <- aggregate(m1d$age, by = list(m1d$milestone, m1d$population, m1d$sex), FUN=length)
colnames(xx) <- c("milestone", "species", "sex", "sample_size")

xx <- xx %>%
  pivot_wider(
    names_from = c(species, sex),
    values_from = sample_size
  )

xx <- aggregate(m2d$age, by = list(m2d$milestone, m2d$population, m2d$sex), FUN=length)
colnames(xx) <- c("milestone", "species", "sex", "sample_size")

xx <- xx %>%
  pivot_wider(
    names_from = c(species, sex),
    values_from = sample_size
  )

# 4 Adjust the ages ####

library(dplyr)

m1d <- m1d %>%
  mutate(
    age_adjust = case_when(
      population == "Bwindi" & sex == "male" ~ age / 15,
      population == "Bwindi" & sex == "female" ~ age / 10.5,
      population == "Loango" & sex == "male" ~ age / 18,
      population == "Loango" & sex == "female" ~ age / 12.2,
      TRUE ~ NA_real_   
    )
  )

# 5 Subset data for specific milestones ####

sort(unique(m1d$milestone))

dorsal_d <- filter(m1d, milestone == "dorsal_50")
spatial_d <- filter(m1d, milestone == "spatial_independence_partner_2_50")
nest_d <- filter(m1d, milestone == "build_nest")
solo_d <- filter(m1d, milestone == "solitary_play")
chest_d <- filter(m1d, milestone == "chest_beat")
mgroom_d <- filter(m1d, milestone == "mother_grooming")
groom_d <- filter(m1d, milestone == "social_grooming")
play_d <- filter(m1d, milestone == "social_play_non_mother")
sex_d <- filter(m1d, milestone == "copulation_play_actor")
agg_d <- filter(m1d, milestone == "aggress_mild")

# 6 Dorsal model a ####

# do data checks

sort(unique(dorsal_d$sex))
length(unique(dorsal_d$id_code))
length(unique(dorsal_d$milestone))
sort(unique(dorsal_d$milestone))
table(dorsal_d$population)
dorsal_d$group_size <- as.numeric(as.character(dorsal_d$group_size))
range(dorsal_d$group_size)
hist(dorsal_d$age)
hist(dorsal_d$age_adjust)
dorsal_d$zsize = (dorsal_d$group_size - mean(dorsal_d$group_size))/sd(dorsal_d$group_size)
range(dorsal_d$zsize)

covees = as.data.frame(
  cbind(
    dorsal_d$zsize, 
    dorsal_d$population, 
    dorsal_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = m1d)
vif(mod)

vif_table <- as.data.frame(vif(mod))
vif_table$model <- "Model 1a"
vif_table$variable <- rownames(vif_table)
vif_table <- vif_table[,c(2,3,1)]

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

m1prior <- get_prior(age ~ zsize + population,
                     data = dorsal_d, family = gaussian())

m1prior = c(prior(normal(0,1), class=b),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population,
              data = dorsal_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
prior_check <- brm(
  age ~ zsize + population,
  data = dorsal_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975), 
  sample_prior = "only",
  silent=2)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 51.90464 secs
mcmc_plot(prior_check)
rm(prior_check)

start_time <- Sys.time()
m1 <- brm(
  age ~ zsize + population,
  data = dorsal_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m1)
summary(m1)
windows()
pp_check(m1) # fine
ppp1 <- pp_check(m1)
ppp1 <- ppp1+labs(title="Model 1: Weakly regularising priors")

mod1_table <- as.data.frame(round(posterior_summary(m1, probs = c(0.05, 0.95)),3))
mod1_table <- mod1_table[c(1:3),]
mod1_table$variable <- rownames(mod1_table)
rownames(mod1_table) <- NULL
mod1_valid <- as.data.frame(summary(m1)$fixed)
mod1_valid$model <- "Model 1a"
mod1_valid <- mod1_valid[,c(6:8)]

m1_coef <- mcmc_plot(m1, variable = "^b_", regex=TRUE)
m1_coef <- m1_coef+labs(title="Model 1: Weakly regularising priors")

# prior sensitivity checks

m1prior_uniform = c(prior(normal(0,100), class=b),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population,
              data = dorsal_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m1uni <- brm(
  age ~ zsize + population,
  data = dorsal_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m1uni)
summary(m1uni)
windows()
pp_check(m1uni) # fine
ppp1uni <- pp_check(m1uni)
ppp1uni <- ppp1uni+labs(title="Model 1: Uniform priors")

m1uni_coef <- mcmc_plot(m1uni, variable = "^b_", regex=TRUE)
m1uni_coef <- m1uni_coef+labs(title="Model 1: Uniform priors")

##

m1prior_strong = c(prior(normal(0,0.5), class=b),
                    prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population,
              data = dorsal_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m1str <- brm(
  age ~ zsize + population,
  data = dorsal_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m1str)
summary(m1str)
windows()
pp_check(m1str) # fine
ppp1str <- pp_check(m1str)
ppp1str <- ppp1str+labs(title="Model 1: Strong priors")


m1str_coef <- mcmc_plot(m1str, variable = "^b_", regex=TRUE)
m1str_coef <- m1str_coef+labs(title="Model 1: Strong priors")

m1_fixed_comparisons <- m1_coef+
  m1uni_coef+
  m1str_coef
windows()
m1_fixed_comparisons

m1_ppcs <- ppp1+ppp1uni+ppp1str
windows()
m1_ppcs

# model passes all checks

# 7 Dorsal model b ####

m1bprior <- get_prior(age_adjust ~ zsize + population,
                     data = dorsal_d, family = Beta())

m1bprior = c(prior(normal(0,1), class=b))

make_stancode(age_adjust ~ zsize + population,
              data = dorsal_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m1b <- brm(
  age_adjust ~ zsize + population,
  data = dorsal_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m1b)
summary(m1b)
windows()
pp_check(m1b) # fine
ppp1b <- pp_check(m1b)
ppp1b <- ppp1+labs(title="Model 11: Weakly regularising priors")

mod1b_table <- as.data.frame(round(posterior_summary(m1b, probs = c(0.05, 0.95)),3))
mod1b_table <- mod1b_table[c(1:3),]
mod1b_table$variable <- rownames(mod1b_table)
rownames(mod1b_table) <- NULL
mod1b_valid <- as.data.frame(summary(m1b)$fixed)
mod1b_valid$model <- "Model 1b"
mod1b_valid <- mod1b_valid[,c(6:8)]

m1b_coef <- mcmc_plot(m1b, variable = "^b_", regex=TRUE)
m1b_coef <- m1b_coef+labs(title="Model 11: Weakly regularising priors")

# prior sensitivity checks

m1bprior_uniform = c(prior(normal(0,100), class=b))

make_stancode(age_adjust ~ zsize + population,
              data = dorsal_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m1buni <- brm(
  age_adjust ~ zsize + population,
  data = dorsal_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m1buni)
summary(m1buni)
windows()
pp_check(m1buni) # fine
ppp1buni <- pp_check(m1buni)
ppp1buni <- ppp1buni+labs(title="Model 11: Uniform priors")

m1buni_coef <- mcmc_plot(m1buni, variable = "^b_", regex=TRUE)
m1buni_coef <- m1buni_coef+labs(title="Model 11: Uniform priors")

##

m1bprior_strong = c(prior(normal(0,0.5), class=b))

make_stancode(age_adjust ~ zsize + population,
              data = dorsal_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m1bstr <- brm(
  age_adjust ~ zsize + population,
  data = dorsal_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m1bstr)
summary(m1bstr)
windows()
pp_check(m1bstr) # fine
ppp1bstr <- pp_check(m1bstr)
ppp1bstr <- ppp1str+labs(title="Model 11: Strong priors")


m1bstr_coef <- mcmc_plot(m1bstr, variable = "^b_", regex=TRUE)
m1bstr_coef <- m1bstr_coef+labs(title="Model 11: Strong priors")

m1b_fixed_comparisons <- m1b_coef+
  m1buni_coef+
  m1bstr_coef
windows()
m1b_fixed_comparisons

m1b_ppcs <- ppp1b+ppp1buni+ppp1bstr
windows()
m1b_ppcs

# model passes all checks

# 8 Spatial model a ####

# do data checks

sort(unique(spatial_d$sex))
length(unique(spatial_d$id_code))
length(unique(spatial_d$milestone))
sort(unique(spatial_d$milestone))
table(spatial_d$population)
spatial_d$group_size <- as.numeric(as.character(spatial_d$group_size))
range(spatial_d$group_size)
hist(spatial_d$age)
hist(spatial_d$age_adjust)
spatial_d$zsize = (spatial_d$group_size - mean(spatial_d$group_size))/sd(spatial_d$group_size)
range(spatial_d$zsize)

covees = as.data.frame(
  cbind(
    spatial_d$zsize, 
    spatial_d$population, 
    spatial_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = spatial_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 2a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

m2prior <- get_prior(age ~ zsize + population,
                     data = spatial_d, family = gaussian())

m2prior = c(prior(normal(0,1), class=b),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population,
              data = spatial_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m2 <- brm(
  age ~ zsize + population,
  data = spatial_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m2)
summary(m2)
windows()
pp_check(m2) # fine
ppp2 <- pp_check(m2)
ppp2 <- ppp2+labs(title="Model 2: Weakly regularising priors")

mod2_table <- as.data.frame(round(posterior_summary(m2, probs = c(0.05, 0.95)),3))
mod2_table <- mod2_table[c(1:3),]
mod2_table$variable <- rownames(mod2_table)
rownames(mod2_table) <- NULL
mod2_valid <- as.data.frame(summary(m2)$fixed)
mod2_valid$model <- "Model 2a"
mod2_valid <- mod2_valid[,c(6:8)]

m2_coef <- mcmc_plot(m2, variable = "^b_", regex=TRUE)
m2_coef <- m2_coef+labs(title="Model 2: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = spatial_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m2uni <- brm(
  age ~ zsize + population,
  data = spatial_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m2uni)
summary(m2uni)
windows()
pp_check(m2uni) # fine
ppp2uni <- pp_check(m2uni)
ppp2uni <- ppp2uni+labs(title="Model 2: Uniform priors")

m2uni_coef <- mcmc_plot(m2uni, variable = "^b_", regex=TRUE)
m2uni_coef <- m2uni_coef+labs(title="Model 2: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = spatial_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m2str <- brm(
  age ~ zsize + population,
  data = spatial_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m2str)
summary(m2str)
windows()
pp_check(m2str) # fine
ppp2str <- pp_check(m2str)
ppp2str <- ppp2str+labs(title="Model 2: Strong priors")


m2str_coef <- mcmc_plot(m2str, variable = "^b_", regex=TRUE)
m2str_coef <- m1str_coef+labs(title="Model 2: Strong priors")

m2_fixed_comparisons <- m2_coef+
  m2uni_coef+
  m2str_coef
windows()
m2_fixed_comparisons

m2_ppcs <- ppp2+ppp2uni+ppp2str
windows()
m2_ppcs

# model passes all checks

# 9 Spatial model b ####

m2bprior <- get_prior(age_adjust ~ zsize + population,
                      data = spatial_d, family = Beta())

m2bprior = c(prior(normal(0,1), class=b))

make_stancode(age_adjust ~ zsize + population,
              data = spatial_d, 
              family = Beta(),
              prior=m2bprior)

start_time <- Sys.time()
m2b <- brm(
  age_adjust ~ zsize + population,
  data = spatial_d, 
  family = Beta(),
  prior=m2bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m2b)
summary(m2b)
windows()
pp_check(m2b) # fine
ppp2b <- pp_check(m2b)
ppp2b <- ppp2+labs(title="Model 12: Weakly regularising priors")

mod2b_table <- as.data.frame(round(posterior_summary(m2b, probs = c(0.05, 0.95)),3))
mod2b_table <- mod2b_table[c(1:3),]
mod2b_table$variable <- rownames(mod2b_table)
rownames(mod2b_table) <- NULL
mod2b_valid <- as.data.frame(summary(m2b)$fixed)
mod2b_valid$model <- "Model 2b"
mod2b_valid <- mod2b_valid[,c(6:8)]

m2b_coef <- mcmc_plot(m2b, variable = "^b_", regex=TRUE)
m2b_coef <- m2b_coef+labs(title="Model 12: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = spatial_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m2buni <- brm(
  age_adjust ~ zsize + population,
  data = spatial_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m2buni)
summary(m2buni)
windows()
pp_check(m2buni) # fine
ppp2buni <- pp_check(m2buni)
ppp2buni <- ppp2buni+labs(title="Model 12: Uniform priors")

m2buni_coef <- mcmc_plot(m2buni, variable = "^b_", regex=TRUE)
m2buni_coef <- m2buni_coef+labs(title="Model 12: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = spatial_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m2bstr <- brm(
  age_adjust ~ zsize + population,
  data = spatial_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m2bstr)
summary(m2bstr)
windows()
pp_check(m2bstr) # fine
ppp2bstr <- pp_check(m2bstr)
ppp2bstr <- ppp2str+labs(title="Model 12: Strong priors")


m2bstr_coef <- mcmc_plot(m2bstr, variable = "^b_", regex=TRUE)
m2bstr_coef <- m2bstr_coef+labs(title="Model 12: Strong priors")

m2b_fixed_comparisons <- m2b_coef+
  m2buni_coef+
  m2bstr_coef
windows()
m2b_fixed_comparisons

m2b_ppcs <- ppp2b+ppp2buni+ppp2bstr
windows()
m2b_ppcs

# model passes all checks

# 10 Nest model a ####

# do data checks

sort(unique(nest_d$sex))
length(unique(nest_d$id_code))
length(unique(nest_d$milestone))
sort(unique(nest_d$milestone))
table(nest_d$population)
nest_d$group_size <- as.numeric(as.character(nest_d$group_size))
range(nest_d$group_size)
hist(nest_d$age)
hist(nest_d$age_adjust)
nest_d$zsize = (nest_d$group_size - mean(nest_d$group_size))/sd(nest_d$group_size)
range(nest_d$zsize)

covees = as.data.frame(
  cbind(
    nest_d$zsize, 
    nest_d$population, 
    nest_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = nest_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 3a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

m3prior <- get_prior(age ~ zsize + population,
                     data = nest_d, family = gaussian())

m3prior = c(prior(normal(0,1), class=b),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population,
              data = nest_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m3 <- brm(
  age ~ zsize + population,
  data = nest_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m3)
summary(m3)
windows()
pp_check(m3) # fine
ppp3 <- pp_check(m3)
ppp3 <- ppp3+labs(title="Model 3: Weakly regularising priors")

mod3_table <- as.data.frame(round(posterior_summary(m3, probs = c(0.05, 0.95)),3))
mod3_table <- mod3_table[c(1:3),]
mod3_table$variable <- rownames(mod3_table)
rownames(mod3_table) <- NULL
mod3_valid <- as.data.frame(summary(m3)$fixed)
mod3_valid$model <- "Model 3a"
mod3_valid <- mod3_valid[,c(6:8)]

m3_coef <- mcmc_plot(m3, variable = "^b_", regex=TRUE)
m3_coef <- m3_coef+labs(title="Model 3: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = nest_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m3uni <- brm(
  age ~ zsize + population,
  data = nest_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m3uni)
summary(m3uni)
windows()
pp_check(m3uni) # fine
ppp3uni <- pp_check(m3uni)
ppp3uni <- ppp2uni+labs(title="Model 3: Uniform priors")

m3uni_coef <- mcmc_plot(m3uni, variable = "^b_", regex=TRUE)
m3uni_coef <- m3uni_coef+labs(title="Model 3: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = nest_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m3str <- brm(
  age ~ zsize + population,
  data = nest_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m3str)
summary(m3str)
windows()
pp_check(m3str) # fine
ppp3str <- pp_check(m3str)
ppp3str <- ppp3str+labs(title="Model 3: Strong priors")


m3str_coef <- mcmc_plot(m3str, variable = "^b_", regex=TRUE)
m3str_coef <- m3str_coef+labs(title="Model 3: Strong priors")

m3_fixed_comparisons <- m3_coef+
  m3uni_coef+
  m3str_coef
windows()
m3_fixed_comparisons

m3_ppcs <- ppp3+ppp3uni+ppp3str
windows()
m3_ppcs

# model passes all checks

# 11 Nest model b ####

m3bprior <- get_prior(age_adjust ~ zsize + population,
                      data = nest_d, family = Beta())

m3bprior = c(prior(normal(0,1), class=b))

make_stancode(age_adjust ~ zsize + population,
              data = nest_d, 
              family = Beta(),
              prior=m3bprior)

start_time <- Sys.time()
m3b <- brm(
  age_adjust ~ zsize + population,
  data = nest_d, 
  family = Beta(),
  prior=m3bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m3b)
summary(m3b)
windows()
pp_check(m3b) # fine
ppp3b <- pp_check(m3b)
ppp3b <- ppp3+labs(title="Model 13: Weakly regularising priors")

mod3b_table <- as.data.frame(round(posterior_summary(m3b, probs = c(0.05, 0.95)),3))
mod3b_table <- mod3b_table[c(1:3),]
mod3b_table$variable <- rownames(mod3b_table)
rownames(mod3b_table) <- NULL
mod3b_valid <- as.data.frame(summary(m3b)$fixed)
mod3b_valid$model <- "Model 3b"
mod3b_valid <- mod3b_valid[,c(6:8)]

m3b_coef <- mcmc_plot(m3b, variable = "^b_", regex=TRUE)
m3b_coef <- m3b_coef+labs(title="Model 13: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = nest_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m3buni <- brm(
  age_adjust ~ zsize + population,
  data = nest_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m3buni)
summary(m3buni)
windows()
pp_check(m3buni) # fine
ppp3buni <- pp_check(m3buni)
ppp3buni <- ppp3buni+labs(title="Model 13: Uniform priors")

m3buni_coef <- mcmc_plot(m3buni, variable = "^b_", regex=TRUE)
m3buni_coef <- m3buni_coef+labs(title="Model 13: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = nest_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m3bstr <- brm(
  age_adjust ~ zsize + population,
  data = nest_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m3bstr)
summary(m3bstr)
windows()
pp_check(m3bstr) # fine
ppp3bstr <- pp_check(m3bstr)
ppp3bstr <- ppp3str+labs(title="Model 13: Strong priors")


m3bstr_coef <- mcmc_plot(m3bstr, variable = "^b_", regex=TRUE)
m3bstr_coef <- m3bstr_coef+labs(title="Model 13: Strong priors")

m3b_fixed_comparisons <- m3b_coef+
  m3buni_coef+
  m3bstr_coef
windows()
m3b_fixed_comparisons

m3b_ppcs <- ppp3b+ppp3buni+ppp3bstr
windows()
m3b_ppcs

# model passes all checks

# 12 Solo play model a ####

# do data checks

sort(unique(solo_d$sex))
length(unique(solo_d$id_code))
length(unique(solo_d$milestone))
sort(unique(solo_d$milestone))
table(solo_d$population)
solo_d$group_size <- as.numeric(as.character(solo_d$group_size))
range(solo_d$group_size)
hist(solo_d$age)
hist(solo_d$age_adjust)
solo_d$zsize = (solo_d$group_size - mean(solo_d$group_size))/sd(solo_d$group_size)
range(solo_d$zsize)

covees = as.data.frame(
  cbind(
    solo_d$zsize, 
    solo_d$population, 
    solo_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = solo_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 4a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

m4prior <- get_prior(age ~ zsize + population,
                     data = solo_d, family = gaussian())

m4prior = c(prior(normal(0,1), class=b),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population,
              data = solo_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m4 <- brm(
  age ~ zsize + population,
  data = solo_d, 
  family = gaussian(),
  prior=m4prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m4)
summary(m4)
windows()
pp_check(m4) # fine
ppp4 <- pp_check(m4)
ppp4 <- ppp4+labs(title="Model 4: Weakly regularising priors")

mod4_table <- as.data.frame(round(posterior_summary(m4, probs = c(0.05, 0.95)),3))
mod4_table <- mod4_table[c(1:3),]
mod4_table$variable <- rownames(mod4_table)
rownames(mod4_table) <- NULL
mod4_valid <- as.data.frame(summary(m4)$fixed)
mod4_valid$model <- "Model 4a"
mod4_valid <- mod4_valid[,c(6:8)]

m4_coef <- mcmc_plot(m4, variable = "^b_", regex=TRUE)
m4_coef <- m4_coef+labs(title="Model 4: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = solo_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m4uni <- brm(
  age ~ zsize + population,
  data = solo_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m4uni)
summary(m4uni)
windows()
pp_check(m4uni) # fine
ppp4uni <- pp_check(m4uni)
ppp4uni <- ppp4uni+labs(title="Model 4: Uniform priors")

m4uni_coef <- mcmc_plot(m4uni, variable = "^b_", regex=TRUE)
m4uni_coef <- m4uni_coef+labs(title="Model 4: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = solo_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m4str <- brm(
  age ~ zsize + population,
  data = solo_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m4str)
summary(m4str)
windows()
pp_check(m4str) # fine
ppp4str <- pp_check(m4str)
ppp4str <- ppp3str+labs(title="Model 4: Strong priors")


m4str_coef <- mcmc_plot(m4str, variable = "^b_", regex=TRUE)
m4str_coef <- m4str_coef+labs(title="Model 4: Strong priors")

m4_fixed_comparisons <- m4_coef+
  m4uni_coef+
  m4str_coef
windows()
m4_fixed_comparisons

m4_ppcs <- ppp4+ppp4uni+ppp4str
windows()
m4_ppcs

# model passes all checks

# 13 Solo play model b ####

m4bprior <- get_prior(age_adjust ~ zsize + population,
                      data = solo_d, family = Beta())

m4bprior = c(prior(normal(0,1), class=b))

make_stancode(age_adjust ~ zsize + population,
              data = solo_d, 
              family = Beta(),
              prior=m4bprior)

start_time <- Sys.time()
m4b <- brm(
  age_adjust ~ zsize + population,
  data = solo_d, 
  family = Beta(),
  prior=m4bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m4b)
summary(m4b)
windows()
pp_check(m4b) # fine
ppp4b <- pp_check(m4b)
ppp4b <- ppp4+labs(title="Model 14: Weakly regularising priors")

mod4b_table <- as.data.frame(round(posterior_summary(m4b, probs = c(0.05, 0.95)),3))
mod4b_table <- mod4b_table[c(1:3),]
mod4b_table$variable <- rownames(mod4b_table)
rownames(mod4b_table) <- NULL
mod4b_valid <- as.data.frame(summary(m4b)$fixed)
mod4b_valid$model <- "Model 4b"
mod4b_valid <- mod4b_valid[,c(6:8)]

m4b_coef <- mcmc_plot(m4b, variable = "^b_", regex=TRUE)
m4b_coef <- m4b_coef+labs(title="Model 14: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = solo_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m4buni <- brm(
  age_adjust ~ zsize + population,
  data = solo_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m4buni)
summary(m4buni)
windows()
pp_check(m4buni) # fine
ppp4buni <- pp_check(m4buni)
ppp4buni <- ppp4buni+labs(title="Model 14: Uniform priors")

m4buni_coef <- mcmc_plot(m4buni, variable = "^b_", regex=TRUE)
m4buni_coef <- m4buni_coef+labs(title="Model 14: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = solo_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m4bstr <- brm(
  age_adjust ~ zsize + population,
  data = solo_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m4bstr)
summary(m4bstr)
windows()
pp_check(m4bstr) # fine
ppp4bstr <- pp_check(m4bstr)
ppp4bstr <- ppp4str+labs(title="Model 14: Strong priors")


m4bstr_coef <- mcmc_plot(m4bstr, variable = "^b_", regex=TRUE)
m4bstr_coef <- m4bstr_coef+labs(title="Model 14: Strong priors")

m4b_fixed_comparisons <- m4b_coef+
  m4buni_coef+
  m4bstr_coef
windows()
m4b_fixed_comparisons

m4b_ppcs <- ppp3b+ppp3buni+ppp3bstr
windows()
m4b_ppcs

# model passes all checks

# 14 Chest model a ####

# do data checks

sort(unique(chest_d$sex))
length(unique(chest_d$id_code))
length(unique(chest_d$milestone))
sort(unique(chest_d$milestone))
table(chest_d$population)
chest_d$group_size <- as.numeric(as.character(chest_d$group_size))
range(chest_d$group_size)
hist(chest_d$age)
hist(chest_d$age_adjust)
chest_d$zsize = (chest_d$group_size - mean(chest_d$group_size))/sd(chest_d$group_size)
range(chest_d$zsize)

covees = as.data.frame(
  cbind(
    chest_d$zsize, 
    chest_d$population, 
    chest_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = chest_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 5a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

make_stancode(age ~ zsize + population,
              data = chest_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m5 <- brm(
  age ~ zsize + population,
  data = chest_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m5)
summary(m5)
windows()
pp_check(m5) # fine
ppp5 <- pp_check(m5)
ppp5 <- ppp5+labs(title="Model 5: Weakly regularising priors")

mod5_table <- as.data.frame(round(posterior_summary(m5, probs = c(0.05, 0.95)),3))
mod5_table <- mod5_table[c(1:3),]
mod5_table$variable <- rownames(mod5_table)
rownames(mod5_table) <- NULL
mod5_valid <- as.data.frame(summary(m5)$fixed)
mod5_valid$model <- "Model 5a"
mod5_valid <- mod5_valid[,c(6:8)]

m5_coef <- mcmc_plot(m5, variable = "^b_", regex=TRUE)
m5_coef <- m5_coef+labs(title="Model 5: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = chest_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m5uni <- brm(
  age ~ zsize + population,
  data = chest_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m5uni)
summary(m5uni)
windows()
pp_check(m5uni) # fine
ppp5uni <- pp_check(m5uni)
ppp5uni <- ppp5uni+labs(title="Model 5: Uniform priors")

m5uni_coef <- mcmc_plot(m5uni, variable = "^b_", regex=TRUE)
m5uni_coef <- m5uni_coef+labs(title="Model 5: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = chest_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m5str <- brm(
  age ~ zsize + population,
  data = chest_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m5str)
summary(m5str)
windows()
pp_check(m5str) # fine
ppp5str <- pp_check(m5str)
ppp5str <- ppp5str+labs(title="Model 5: Strong priors")


m5str_coef <- mcmc_plot(m5str, variable = "^b_", regex=TRUE)
m5str_coef <- m5str_coef+labs(title="Model 5: Strong priors")

m5_fixed_comparisons <- m5_coef+
  m5uni_coef+
  m5str_coef
windows()
m5_fixed_comparisons

m5_ppcs <- ppp5+ppp5uni+ppp5str
windows()
m5_ppcs

# model passes all checks

# 15 Chest model b ####

make_stancode(age_adjust ~ zsize + population,
              data = chest_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m5b <- brm(
  age_adjust ~ zsize + population,
  data = chest_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m5b)
summary(m5b)
windows()
pp_check(m5b) # fine
ppp5b <- pp_check(m5b)
ppp5b <- ppp5+labs(title="Model 15: Weakly regularising priors")

mod5b_table <- as.data.frame(round(posterior_summary(m5b, probs = c(0.05, 0.95)),3))
mod5b_table <- mod5b_table[c(1:3),]
mod5b_table$variable <- rownames(mod5b_table)
rownames(mod5b_table) <- NULL
mod5b_valid <- as.data.frame(summary(m5b)$fixed)
mod5b_valid$model <- "Model 5b"
mod5b_valid <- mod5b_valid[,c(6:8)]

m5b_coef <- mcmc_plot(m5b, variable = "^b_", regex=TRUE)
m5b_coef <- m5b_coef+labs(title="Model 15: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = chest_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m5buni <- brm(
  age_adjust ~ zsize + population,
  data = chest_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m5buni)
summary(m5buni)
windows()
pp_check(m5buni) # fine
ppp5buni <- pp_check(m5buni)
ppp5buni <- ppp5buni+labs(title="Model 15: Uniform priors")

m5buni_coef <- mcmc_plot(m5buni, variable = "^b_", regex=TRUE)
m5buni_coef <- m5buni_coef+labs(title="Model 15: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = chest_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m5bstr <- brm(
  age_adjust ~ zsize + population,
  data = chest_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m5bstr)
summary(m5bstr)
windows()
pp_check(m5bstr) # fine
ppp5bstr <- pp_check(m5bstr)
ppp5bstr <- ppp5str+labs(title="Model 15: Strong priors")


m5bstr_coef <- mcmc_plot(m5bstr, variable = "^b_", regex=TRUE)
m5bstr_coef <- m5bstr_coef+labs(title="Model 15: Strong priors")

m5b_fixed_comparisons <- m5b_coef+
  m5buni_coef+
  m5bstr_coef
windows()
m5b_fixed_comparisons

m5b_ppcs <- ppp5b+ppp5buni+ppp5bstr
windows()
m5b_ppcs

# 16 Mother groom model a ####

# do data checks

sort(unique(mgroom_d$sex))
length(unique(mgroom_d$id_code))
length(unique(mgroom_d$milestone))
sort(unique(mgroom_d$milestone))
table(mgroom_d$population)
mgroom_d$group_size <- as.numeric(as.character(mgroom_d$group_size))
range(mgroom_d$group_size)
hist(mgroom_d$age)
hist(mgroom_d$age_adjust)
mgroom_d$zsize = (mgroom_d$group_size - mean(mgroom_d$group_size))/sd(mgroom_d$group_size)
range(mgroom_d$zsize)

covees = as.data.frame(
  cbind(
    mgroom_d$zsize, 
    mgroom_d$population, 
    mgroom_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = mgroom_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 6a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

make_stancode(age ~ zsize + population,
              data = mgroom_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m6 <- brm(
  age ~ zsize + population,
  data = mgroom_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m6)
summary(m6)
windows()
pp_check(m6) # fine
ppp6 <- pp_check(m6)
ppp6 <- ppp6+labs(title="Model 6: Weakly regularising priors")

mod6_table <- as.data.frame(round(posterior_summary(m6, probs = c(0.05, 0.95)),3))
mod6_table <- mod6_table[c(1:3),]
mod6_table$variable <- rownames(mod6_table)
rownames(mod6_table) <- NULL
mod6_valid <- as.data.frame(summary(m6)$fixed)
mod6_valid$model <- "Model 6a"
mod6_valid <- mod6_valid[,c(6:8)]

m6_coef <- mcmc_plot(m6, variable = "^b_", regex=TRUE)
m6_coef <- m6_coef+labs(title="Model 6: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = mgroom_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m6uni <- brm(
  age ~ zsize + population,
  data = mgroom_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m6uni)
summary(m6uni)
windows()
pp_check(m6uni) # fine
ppp6uni <- pp_check(m6uni)
ppp6uni <- ppp6uni+labs(title="Model 6: Uniform priors")

m6uni_coef <- mcmc_plot(m6uni, variable = "^b_", regex=TRUE)
m6uni_coef <- m6uni_coef+labs(title="Model 6: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = mgroom_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m6str <- brm(
  age ~ zsize + population,
  data = mgroom_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m6str)
summary(m6str)
windows()
pp_check(m6str) # fine
ppp6str <- pp_check(m6str)
ppp6str <- ppp6str+labs(title="Model 6: Strong priors")


m6str_coef <- mcmc_plot(m6str, variable = "^b_", regex=TRUE)
m6str_coef <- m6str_coef+labs(title="Model 6: Strong priors")

m6_fixed_comparisons <- m6_coef+
  m6uni_coef+
  m6str_coef
windows()
m6_fixed_comparisons

m6_ppcs <- ppp6+ppp6uni+ppp6str
windows()
m6_ppcs

# model passes all checks

# 17 Mother groom model b ####

make_stancode(age_adjust ~ zsize + population,
              data = mgroom_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m6b <- brm(
  age_adjust ~ zsize + population,
  data = mgroom_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m6b)
summary(m6b)
windows()
pp_check(m6b) # fine
ppp6b <- pp_check(m6b)
ppp6b <- ppp6+labs(title="Model 16: Weakly regularising priors")

mod6b_table <- as.data.frame(round(posterior_summary(m6b, probs = c(0.05, 0.95)),3))
mod6b_table <- mod6b_table[c(1:3),]
mod6b_table$variable <- rownames(mod6b_table)
rownames(mod6b_table) <- NULL
mod6b_valid <- as.data.frame(summary(m6b)$fixed)
mod6b_valid$model <- "Model 6b"
mod6b_valid <- mod6b_valid[,c(6:8)]

m6b_coef <- mcmc_plot(m6b, variable = "^b_", regex=TRUE)
m6b_coef <- m6b_coef+labs(title="Model 16: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = mgroom_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m6buni <- brm(
  age_adjust ~ zsize + population,
  data = mgroom_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m6buni)
summary(m6buni)
windows()
pp_check(m6buni) # fine
ppp6buni <- pp_check(m6buni)
ppp6buni <- ppp6buni+labs(title="Model 16: Uniform priors")

m6buni_coef <- mcmc_plot(m6buni, variable = "^b_", regex=TRUE)
m6buni_coef <- m6buni_coef+labs(title="Model 16: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = mgroom_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m6bstr <- brm(
  age_adjust ~ zsize + population,
  data = mgroom_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m6bstr)
summary(m6bstr)
windows()
pp_check(m6bstr) # fine
ppp6bstr <- pp_check(m6bstr)
ppp6bstr <- ppp6str+labs(title="Model 16: Strong priors")


m6bstr_coef <- mcmc_plot(m6bstr, variable = "^b_", regex=TRUE)
m6bstr_coef <- m6bstr_coef+labs(title="Model 16: Strong priors")

m6b_fixed_comparisons <- m6b_coef+
  m6buni_coef+
  m6bstr_coef
windows()
m6b_fixed_comparisons

m6b_ppcs <- ppp6b+ppp6buni+ppp6bstr
windows()
m6b_ppcs

# 18 Groom model a ####

# do data checks

sort(unique(groom_d$sex))
length(unique(groom_d$id_code))
length(unique(groom_d$milestone))
sort(unique(groom_d$milestone))
table(groom_d$population)
groom_d$group_size <- as.numeric(as.character(groom_d$group_size))
range(groom_d$group_size)
hist(groom_d$age)
hist(groom_d$age_adjust)
groom_d$zsize = (groom_d$group_size - mean(groom_d$group_size))/sd(groom_d$group_size)
range(groom_d$zsize)

covees = as.data.frame(
  cbind(
    groom_d$zsize, 
    groom_d$population, 
    groom_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = groom_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 7a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

make_stancode(age ~ zsize + population,
              data = groom_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m7 <- brm(
  age ~ zsize + population,
  data = groom_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m7)
summary(m7)
windows()
pp_check(m7) # fine
ppp7 <- pp_check(m7)
ppp7 <- ppp7+labs(title="Model 7: Weakly regularising priors")

mod7_table <- as.data.frame(round(posterior_summary(m7, probs = c(0.05, 0.95)),3))
mod7_table <- mod7_table[c(1:3),]
mod7_table$variable <- rownames(mod7_table)
rownames(mod7_table) <- NULL
mod7_valid <- as.data.frame(summary(m7)$fixed)
mod7_valid$model <- "Model 7a"
mod7_valid <- mod7_valid[,c(6:8)]

m7_coef <- mcmc_plot(m7, variable = "^b_", regex=TRUE)
m7_coef <- m7_coef+labs(title="Model 7: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = groom_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m7uni <- brm(
  age ~ zsize + population,
  data = groom_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m7uni)
summary(m7uni)
windows()
pp_check(m7uni) # fine
ppp7uni <- pp_check(m7uni)
ppp7uni <- ppp7uni+labs(title="Model 7: Uniform priors")

m7uni_coef <- mcmc_plot(m7uni, variable = "^b_", regex=TRUE)
m7uni_coef <- m7uni_coef+labs(title="Model 7: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = groom_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m7str <- brm(
  age ~ zsize + population,
  data = groom_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m7str)
summary(m7str)
windows()
pp_check(m7str) # fine
ppp7str <- pp_check(m7str)
ppp7str <- ppp6str+labs(title="Model 7: Strong priors")


m7str_coef <- mcmc_plot(m7str, variable = "^b_", regex=TRUE)
m7str_coef <- m7str_coef+labs(title="Model 7: Strong priors")

m7_fixed_comparisons <- m7_coef+
  m7uni_coef+
  m7str_coef
windows()
m7_fixed_comparisons

m7_ppcs <- ppp7+ppp7uni+ppp7str
windows()
m7_ppcs

# model passes all checks

# 19 Groom model b ####

make_stancode(age_adjust ~ zsize + population,
              data = groom_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m7b <- brm(
  age_adjust ~ zsize + population,
  data = groom_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m7b)
summary(m7b)
windows()
pp_check(m7b) # fine
ppp7b <- pp_check(m7b)
ppp7b <- ppp7+labs(title="Model 17: Weakly regularising priors")

mod7b_table <- as.data.frame(round(posterior_summary(m7b, probs = c(0.05, 0.95)),3))
mod7b_table <- mod7b_table[c(1:3),]
mod7b_table$variable <- rownames(mod7b_table)
rownames(mod7b_table) <- NULL
mod7b_valid <- as.data.frame(summary(m7b)$fixed)
mod7b_valid$model <- "Model 7b"
mod7b_valid <- mod7b_valid[,c(6:8)]

m7b_coef <- mcmc_plot(m7b, variable = "^b_", regex=TRUE)
m7b_coef <- m7b_coef+labs(title="Model 17: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = groom_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m7buni <- brm(
  age_adjust ~ zsize + population,
  data = groom_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m7buni)
summary(m7buni)
windows()
pp_check(m7buni) # fine
ppp7buni <- pp_check(m7buni)
ppp7buni <- ppp7buni+labs(title="Model 17: Uniform priors")

m7buni_coef <- mcmc_plot(m7buni, variable = "^b_", regex=TRUE)
m7buni_coef <- m7buni_coef+labs(title="Model 17: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = groom_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m7bstr <- brm(
  age_adjust ~ zsize + population,
  data = groom_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m7bstr)
summary(m7bstr)
windows()
pp_check(m7bstr) # fine
ppp7bstr <- pp_check(m7bstr)
ppp7bstr <- ppp7str+labs(title="Model 17: Strong priors")


m7bstr_coef <- mcmc_plot(m7bstr, variable = "^b_", regex=TRUE)
m7bstr_coef <- m7bstr_coef+labs(title="Model 17: Strong priors")

m7b_fixed_comparisons <- m7b_coef+
  m7buni_coef+
  m7bstr_coef
windows()
m7b_fixed_comparisons

m7b_ppcs <- ppp7b+ppp7buni+ppp7bstr
windows()
m7b_ppcs

# 20 Play model a ####

# do data checks

sort(unique(play_d$sex))
length(unique(play_d$id_code))
length(unique(play_d$milestone))
sort(unique(play_d$milestone))
table(play_d$population)
play_d$group_size <- as.numeric(as.character(play_d$group_size))
range(play_d$group_size)
hist(play_d$age)
hist(play_d$age_adjust)
play_d$zsize = (play_d$group_size - mean(play_d$group_size))/sd(play_d$group_size)
range(play_d$zsize)

covees = as.data.frame(
  cbind(
    play_d$zsize, 
    play_d$population, 
    play_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = play_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 8a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

make_stancode(age ~ zsize + population,
              data = play_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m8 <- brm(
  age ~ zsize + population,
  data = play_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m8)
summary(m8)
windows()
pp_check(m8) # fine
ppp8 <- pp_check(m8)
ppp8 <- ppp8+labs(title="Model 8: Weakly regularising priors")

mod8_table <- as.data.frame(round(posterior_summary(m8, probs = c(0.05, 0.95)),3))
mod8_table <- mod8_table[c(1:3),]
mod8_table$variable <- rownames(mod8_table)
rownames(mod8_table) <- NULL
mod8_valid <- as.data.frame(summary(m8)$fixed)
mod8_valid$model <- "Model 8a"
mod8_valid <- mod8_valid[,c(6:8)]

m8_coef <- mcmc_plot(m8, variable = "^b_", regex=TRUE)
m8_coef <- m8_coef+labs(title="Model 8: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = play_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m8uni <- brm(
  age ~ zsize + population,
  data = play_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m8uni)
summary(m8uni)
windows()
pp_check(m8uni) # fine
ppp8uni <- pp_check(m8uni)
ppp8uni <- ppp8uni+labs(title="Model 8: Uniform priors")

m8uni_coef <- mcmc_plot(m8uni, variable = "^b_", regex=TRUE)
m8uni_coef <- m8uni_coef+labs(title="Model 8: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = play_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m8str <- brm(
  age ~ zsize + population,
  data = play_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m8str)
summary(m8str)
windows()
pp_check(m8str) # fine
ppp8str <- pp_check(m8str)
ppp8str <- ppp8str+labs(title="Model 8: Strong priors")


m8str_coef <- mcmc_plot(m8str, variable = "^b_", regex=TRUE)
m8str_coef <- m8str_coef+labs(title="Model 8: Strong priors")

m8_fixed_comparisons <- m8_coef+
  m8uni_coef+
  m8str_coef
windows()
m8_fixed_comparisons

m8_ppcs <- ppp8+ppp8uni+ppp8str
windows()
m8_ppcs

# model passes all checks

# 21 Play model b ####

make_stancode(age_adjust ~ zsize + population,
              data = play_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m8b <- brm(
  age_adjust ~ zsize + population,
  data = play_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m8b)
summary(m8b)
windows()
pp_check(m8b) # fine
ppp8b <- pp_check(m8b)
ppp8b <- ppp8+labs(title="Model 18: Weakly regularising priors")

mod8b_table <- as.data.frame(round(posterior_summary(m8b, probs = c(0.05, 0.95)),3))
mod8b_table <- mod8b_table[c(1:3),]
mod8b_table$variable <- rownames(mod8b_table)
rownames(mod8b_table) <- NULL
mod8b_valid <- as.data.frame(summary(m8b)$fixed)
mod8b_valid$model <- "Model 8b"
mod8b_valid <- mod8b_valid[,c(6:8)]

m8b_coef <- mcmc_plot(m8b, variable = "^b_", regex=TRUE)
m8b_coef <- m8b_coef+labs(title="Model 18: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = play_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m8buni <- brm(
  age_adjust ~ zsize + population,
  data = play_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m8buni)
summary(m8buni)
windows()
pp_check(m8buni) # fine
ppp8buni <- pp_check(m8buni)
ppp8buni <- ppp8buni+labs(title="Model 18: Uniform priors")

m8buni_coef <- mcmc_plot(m8buni, variable = "^b_", regex=TRUE)
m8buni_coef <- m8buni_coef+labs(title="Model 18: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = play_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m8bstr <- brm(
  age_adjust ~ zsize + population,
  data = play_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m8bstr)
summary(m8bstr)
windows()
pp_check(m8bstr) # fine
ppp8bstr <- pp_check(m8bstr)
ppp8bstr <- ppp8str+labs(title="Model 18: Strong priors")


m8bstr_coef <- mcmc_plot(m8bstr, variable = "^b_", regex=TRUE)
m8bstr_coef <- m8bstr_coef+labs(title="Model 18: Strong priors")

m8b_fixed_comparisons <- m8b_coef+
  m8buni_coef+
  m8bstr_coef
windows()
m8b_fixed_comparisons

m8b_ppcs <- ppp8b+ppp8buni+ppp8bstr
windows()
m8b_ppcs

# 22 Copulation model a ####

# do data checks

sort(unique(sex_d$sex))
length(unique(sex_d$id_code))
length(unique(sex_d$milestone))
sort(unique(sex_d$milestone))
table(sex_d$population)
sex_d$group_size <- as.numeric(as.character(sex_d$group_size))
range(sex_d$group_size)
hist(sex_d$age)
hist(sex_d$age_adjust)
sex_d$zsize = (sex_d$group_size - mean(sex_d$group_size))/sd(sex_d$group_size)
range(sex_d$zsize)

covees = as.data.frame(
  cbind(
    sex_d$zsize, 
    sex_d$population, 
    sex_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = sex_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 9a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

make_stancode(age ~ zsize + population,
              data = sex_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m9 <- brm(
  age ~ zsize + population,
  data = sex_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m9)
summary(m9)
windows()
pp_check(m9) # fine
ppp9 <- pp_check(m9)
ppp9 <- ppp9+labs(title="Model 9: Weakly regularising priors")

mod9_table <- as.data.frame(round(posterior_summary(m9, probs = c(0.05, 0.95)),3))
mod9_table <- mod9_table[c(1:3),]
mod9_table$variable <- rownames(mod9_table)
rownames(mod9_table) <- NULL
mod9_valid <- as.data.frame(summary(m9)$fixed)
mod9_valid$model <- "Model 9a"
mod9_valid <- mod9_valid[,c(6:8)]

m9_coef <- mcmc_plot(m9, variable = "^b_", regex=TRUE)
m9_coef <- m9_coef+labs(title="Model 9: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = sex_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m9uni <- brm(
  age ~ zsize + population,
  data = sex_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m9uni)
summary(m9uni)
windows()
pp_check(m9uni) # fine
ppp9uni <- pp_check(m9uni)
ppp9uni <- ppp9uni+labs(title="Model 9: Uniform priors")

m9uni_coef <- mcmc_plot(m9uni, variable = "^b_", regex=TRUE)
m9uni_coef <- m9uni_coef+labs(title="Model 9: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = sex_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m9str <- brm(
  age ~ zsize + population,
  data = sex_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m9str)
summary(m9str)
windows()
pp_check(m9str) # fine
ppp9str <- pp_check(m9str)
ppp9str <- ppp9str+labs(title="Model 9: Strong priors")


m9str_coef <- mcmc_plot(m9str, variable = "^b_", regex=TRUE)
m9str_coef <- m9str_coef+labs(title="Model 9: Strong priors")

m9_fixed_comparisons <- m9_coef+
  m9uni_coef+
  m9str_coef
windows()
m9_fixed_comparisons

m9_ppcs <- ppp9+ppp9uni+ppp9str
windows()
m9_ppcs

# model passes all checks

# 23 Copulation model b ####

make_stancode(age_adjust ~ zsize + population,
              data = sex_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m9b <- brm(
  age_adjust ~ zsize + population,
  data = sex_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m9b)
summary(m9b)
windows()
pp_check(m9b) # fine
ppp9b <- pp_check(m9b)
ppp9b <- ppp9+labs(title="Model 19: Weakly regularising priors")

mod9b_table <- as.data.frame(round(posterior_summary(m9b, probs = c(0.05, 0.95)),3))
mod9b_table <- mod9b_table[c(1:3),]
mod9b_table$variable <- rownames(mod9b_table)
rownames(mod9b_table) <- NULL
mod9b_valid <- as.data.frame(summary(m9b)$fixed)
mod9b_valid$model <- "Model 9b"
mod9b_valid <- mod9b_valid[,c(6:8)]

m9b_coef <- mcmc_plot(m9b, variable = "^b_", regex=TRUE)
m9b_coef <- m9b_coef+labs(title="Model 19: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = sex_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m9buni <- brm(
  age_adjust ~ zsize + population,
  data = sex_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m9buni)
summary(m9buni)
windows()
pp_check(m9buni) # fine
ppp9buni <- pp_check(m9buni)
ppp9buni <- ppp9buni+labs(title="Model 19: Uniform priors")

m9buni_coef <- mcmc_plot(m9buni, variable = "^b_", regex=TRUE)
m9buni_coef <- m9buni_coef+labs(title="Model 19: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = sex_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m9bstr <- brm(
  age_adjust ~ zsize + population,
  data = sex_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m9bstr)
summary(m9bstr)
windows()
pp_check(m9bstr) # fine
ppp9bstr <- pp_check(m9bstr)
ppp9bstr <- ppp9str+labs(title="Model 19: Strong priors")


m9bstr_coef <- mcmc_plot(m9bstr, variable = "^b_", regex=TRUE)
m9bstr_coef <- m9bstr_coef+labs(title="Model 19: Strong priors")

m9b_fixed_comparisons <- m9b_coef+
  m9buni_coef+
  m9bstr_coef
windows()
m9b_fixed_comparisons

m9b_ppcs <- ppp9b+ppp9buni+ppp9bstr
windows()
m9b_ppcs

# 24 Agg model a ####

# do data checks

sort(unique(agg_d$sex))
length(unique(agg_d$id_code))
length(unique(agg_d$milestone))
sort(unique(agg_d$milestone))
table(agg_d$population)
agg_d$group_size <- as.numeric(as.character(agg_d$group_size))
range(agg_d$group_size)
hist(agg_d$age)
hist(agg_d$age_adjust)
agg_d$zsize = (agg_d$group_size - mean(agg_d$group_size))/sd(agg_d$group_size)
range(agg_d$zsize)

covees = as.data.frame(
  cbind(
    agg_d$zsize, 
    agg_d$population, 
    agg_d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "sex")

mod = lm(age ~ zsize + 
           population, 
         data = agg_d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 10a"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# let's check prior

make_stancode(age ~ zsize + population,
              data = agg_d, 
              family = gaussian(),
              prior=m1prior)

start_time <- Sys.time()
m10 <- brm(
  age ~ zsize + population,
  data = agg_d, 
  family = gaussian(),
  prior=m1prior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m10)
summary(m10)
windows()
pp_check(m10) # fine
ppp10 <- pp_check(m10)
ppp10 <- ppp10+labs(title="Model 10: Weakly regularising priors")

mod10_table <- as.data.frame(round(posterior_summary(m10, probs = c(0.05, 0.95)),3))
mod10_table <- mod10_table[c(1:3),]
mod10_table$variable <- rownames(mod10_table)
rownames(mod10_table) <- NULL
mod10_valid <- as.data.frame(summary(m10)$fixed)
mod10_valid$model <- "Model 10a"
mod10_valid <- mod10_valid[,c(6:8)]

m10_coef <- mcmc_plot(m10, variable = "^b_", regex=TRUE)
m10_coef <- m10_coef+labs(title="Model 10: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age ~ zsize + population,
              data = agg_d,
              family = gaussian(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m10uni <- brm(
  age ~ zsize + population,
  data = agg_d, 
  family = gaussian(),
  prior=m1prior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m10uni)
summary(m10uni)
windows()
pp_check(m10uni) # fine
ppp10uni <- pp_check(m10uni)
ppp10uni <- ppp10uni+labs(title="Model 10: Uniform priors")

m10uni_coef <- mcmc_plot(m10uni, variable = "^b_", regex=TRUE)
m10uni_coef <- m10uni_coef+labs(title="Model 10: Uniform priors")

##

make_stancode(age ~ zsize + population,
              data = agg_d, 
              family = gaussian(),
              prior=m1prior_strong)

start_time <- Sys.time()
m10str <- brm(
  age ~ zsize + population,
  data = agg_d, 
  family = gaussian(),
  prior=m1prior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m10str)
summary(m10str)
windows()
pp_check(m10str) # fine
ppp10str <- pp_check(m10str)
ppp10str <- ppp10str+labs(title="Model 10: Strong priors")


m10str_coef <- mcmc_plot(m10str, variable = "^b_", regex=TRUE)
m10str_coef <- m10str_coef+labs(title="Model 10: Strong priors")

m10_fixed_comparisons <- m10_coef+
  m10uni_coef+
  m10str_coef
windows()
m10_fixed_comparisons

m10_ppcs <- ppp10+ppp10uni+ppp10str
windows()
m10_ppcs

# model passes all checks

# 25 Agg model b ####

make_stancode(age_adjust ~ zsize + population,
              data = agg_d, 
              family = Beta(),
              prior=m1bprior)

start_time <- Sys.time()
m10b <- brm(
  age_adjust ~ zsize + population,
  data = agg_d, 
  family = Beta(),
  prior=m1bprior,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m10b)
summary(m10b)
windows()
pp_check(m10b) # fine
ppp10b <- pp_check(m10b)
ppp10b <- ppp10+labs(title="Model 20: Weakly regularising priors")

mod10b_table <- as.data.frame(round(posterior_summary(m10b, probs = c(0.05, 0.95)),3))
mod10b_table <- mod10b_table[c(1:3),]
mod10b_table$variable <- rownames(mod10b_table)
rownames(mod10b_table) <- NULL
mod10b_valid <- as.data.frame(summary(m10b)$fixed)
mod10b_valid$model <- "Model 10b"
mod10b_valid <- mod10b_valid[,c(6:8)]

m10b_coef <- mcmc_plot(m10b, variable = "^b_", regex=TRUE)
m10b_coef <- m10b_coef+labs(title="Model 20: Weakly regularising priors")

# prior sensitivity checks

make_stancode(age_adjust ~ zsize + population,
              data = agg_d,
              family = Beta(),
              prior=m1bprior_uniform)

start_time <- Sys.time()
m10buni <- brm(
  age_adjust ~ zsize + population,
  data = agg_d, 
  family = Beta(),
  prior=m1bprior_uniform,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m10buni)
summary(m10buni)
windows()
pp_check(m10buni) # fine
ppp10buni <- pp_check(m10buni)
ppp10buni <- ppp10buni+labs(title="Model 20: Uniform priors")

m10buni_coef <- mcmc_plot(m10buni, variable = "^b_", regex=TRUE)
m10buni_coef <- m10buni_coef+labs(title="Model 20: Uniform priors")

##

make_stancode(age_adjust ~ zsize + population,
              data = agg_d, 
              family = Beta(),
              prior=m1bprior_strong)

start_time <- Sys.time()
m10bstr <- brm(
  age_adjust ~ zsize + population,
  data = agg_d, 
  family = Beta(),
  prior=m1bprior_strong,
  chains = 3, 
  cores = 3,
  backend = "cmdstanr",
  stan_model_args = list(stanc_options=list("O1")),
  iter = 5000, 
  warmup = 2500, 
  thin = 2,
  control = list(adapt_delta = 0.975),
  silent=0)
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of  1.401607 mins

windows()
plot(m10bstr)
summary(m10bstr)
windows()
pp_check(m10bstr) # fine
ppp10bstr <- pp_check(m10bstr)
ppp10bstr <- ppp10str+labs(title="Model 20: Strong priors")


m10bstr_coef <- mcmc_plot(m10bstr, variable = "^b_", regex=TRUE)
m10bstr_coef <- m10bstr_coef+labs(title="Model 20: Strong priors")

m10b_fixed_comparisons <- m10b_coef+
  m10buni_coef+
  m10bstr_coef
windows()
m10b_fixed_comparisons

m10b_ppcs <- ppp10b+ppp10buni+ppp10bstr
windows()
m10b_ppcs

# 26 Create summary of first set of models ####

age_mods <- rbind(
  mod1_table,
  mod2_table,
  mod3_table,
  mod4_table,
  mod5_table,
  mod6_table,
  mod7_table,
  mod8_table,
  mod9_table,
  mod10_table)

age_adjusted_mods <- rbind(
  mod1b_table,
  mod2b_table,
  mod3b_table,
  mod4b_table,
  mod5b_table,
  mod6b_table,
  mod7b_table,
  mod8b_table,
  mod9b_table,
  mod10b_table)

# 26 Summarise Models 1a and 1b ####

mod1_table
mod1a_table

# proportion support

post <- posterior_samples(m1)
names(post)

round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m1b)
names(post)

round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m8)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m8b)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m9)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m9b)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m4)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m4b)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m5b)
round(sum(post$`b_populationLoango` < 0) / 
        length(post$`b_populationLoango`),3)

post <- posterior_samples(m3)
round(sum(post$`b_populationLoango` > 0) / 
        length(post$`b_populationLoango`),3)

# 27 Figure 1 ####

xx <- met.brewer(name="Veronese", n=2, type="discrete")

plot1a <- ggplot(fig1_d,
                aes(x = reorder(milestone,-mean),
                    y = mean)) + 
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd,
                    color=species),
                position = position_dodge(width = 1)) +
  geom_point(color='black',
             shape=21,
             size=4,
             aes(fill=species),
             position=position_dodge(width=1))+
  scale_fill_manual(values=xx)+
  scale_color_manual(values=xx)+
  coord_flip()+
  theme_bw()+
  xlab("Milestone")+
  ylab("Mean age of emergence")+
  ggtitle("A")+
  facet_wrap(~species)+
  #scale_x_discrete(labels = function(x) 
  #stringr::str_wrap(x, width = 7))+
  theme(axis.title=element_text(face="bold",size=15),
        axis.text=element_text(size=15),
        axis.text.x=element_text(angle = 90, vjust = 0.5),
        plot.title=element_text(face="bold",size=15),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        strip.text=element_text(face="bold", size=15))
windows()  
plot1a

yy <- met.brewer(name="Demuth", n=5, type="discrete")

plot1b <- ggplot(fig1_d, aes(x = factor(category), 
                             y = mean, 
                             fill = factor(category))) +
  geom_boxplot(position=position_dodge(width=0.9),
              #width=1.1,
              trim=FALSE) +
  geom_jitter(data=fig1_d,
               aes(x=factor(category),
                   y=mean,
                   fill=factor(category)),
               #position=position_dodge(width=0.9),
               width=0.075,
               outlier.shape = NA)+
  scale_fill_manual(values = yy) +
  xlab("Milestone Category") +
  ylab("Age of emergence") +
  ggtitle("B")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
plot1b

zz <- met.brewer(name="Egypt", n=3, type="discrete")

plot1c <- ggplot(fig1_d, aes(x = factor(complexity), 
                             y = mean, 
                             fill = factor(complexity))) +
  geom_boxplot(position=position_dodge(width=0.9),
               #width=1.1,
               trim=FALSE) +
  geom_jitter(data=fig1_d,
              aes(x=factor(complexity),
                  y=mean,
                  fill=factor(complexity)),
              #position=position_dodge(width=0.9),
              width=0.075,
              outlier.shape = NA)+
  scale_fill_manual(values = zz) +
  xlab("Complexity") +
  ylab("Age of emergence") +
  ggtitle("C")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
plot1c


fig1 <- plot1a+plot1b+plot1c
windows()
fig1

# 28 Figure 2 ####

xx <- met.brewer(name="Veronese", n=2, type="discrete")
dorsal_d$species <- ifelse(dorsal_d$population=="Bwindi","Mountain","Western")

plot2a <- ggplot(dorsal_d, aes(x = factor(species), 
                             y = age_adjust, 
                             fill = factor(species))) +
  geom_boxplot(position=position_dodge(width=0.9),
               #width=1.1,
               trim=FALSE) +
  geom_jitter(data=dorsal_d,
              aes(x=factor(species),
                  y=age_adjust,
                  fill=factor(species)),
              #position=position_dodge(width=0.9),
              width=0.075,
              outlier.shape = NA)+
  scale_fill_manual(values = xx) +
  xlab("Species") +
  ylab("Emergence as proportion of age at sexual maturation") +
  ggtitle("(A) Dorsal Carry 50%")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.position = "none")
windows()
plot2a

solo_d$species <- ifelse(solo_d$population=="Bwindi","Mountain","Western")

plot2b <- ggplot(solo_d, aes(x = factor(species), 
                               y = age_adjust, 
                               fill = factor(species))) +
  geom_boxplot(position=position_dodge(width=0.9),
               #width=1.1,
               trim=FALSE) +
  geom_jitter(data=solo_d,
              aes(x=factor(species),
                  y=age_adjust,
                  fill=factor(species)),
              #position=position_dodge(width=0.9),
              width=0.075,
              outlier.shape = NA)+
  scale_fill_manual(values = xx) +
  xlab("Species") +
  ylab("Emergence as proportion of age at sexual maturation") +
  ggtitle("(B) Solitary Play")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.position = "none")
windows()
plot2b

chest_d$species <- ifelse(chest_d$population=="Bwindi","Mountain","Western")

plot2c <- ggplot(chest_d, aes(x = factor(species), 
                             y = age_adjust, 
                             fill = factor(species))) +
  geom_boxplot(position=position_dodge(width=0.9),
               #width=1.1,
               trim=FALSE) +
  geom_jitter(data=chest_d,
              aes(x=factor(species),
                  y=age_adjust,
                  fill=factor(species)),
              #position=position_dodge(width=0.9),
              width=0.075,
              outlier.shape = NA)+
  scale_fill_manual(values = xx) +
  xlab("Species") +
  ylab("Emergence as proportion of age at sexual maturation") +
  ggtitle("(C) Chest Beat")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
plot2c


play_d$species <- ifelse(play_d$population=="Bwindi","Mountain","Western")

plot2d <- ggplot(play_d, aes(x = factor(species), 
                               y = age_adjust, 
                               fill = factor(species))) +
  geom_boxplot(position=position_dodge(width=0.9),
               #width=1.1,
               trim=FALSE) +
  geom_jitter(data=play_d,
              aes(x=factor(species),
                  y=age_adjust,
                  fill=factor(species)),
              #position=position_dodge(width=0.9),
              width=0.075,
              outlier.shape = NA)+
  scale_fill_manual(values = xx) +
  xlab("Species") +
  ylab("Emergence as proportion of age at sexual maturation") +
  ggtitle("(D) Social Play")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_bw() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.position = "none")
windows()
plot2d

sex_d$species <- ifelse(sex_d$population=="Bwindi","Mountain","Western")

plot2e <- ggplot(sex_d, aes(x = factor(species), 
                             y = age_adjust, 
                             fill = factor(species))) +
  geom_boxplot(position=position_dodge(width=0.9),
               #width=1.1,
               trim=FALSE) +
  geom_jitter(data=sex_d,
              aes(x=factor(species),
                  y=age_adjust,
                  fill=factor(species)),
              #position=position_dodge(width=0.9),
              width=0.075,
              outlier.shape = NA)+
  scale_fill_manual(values = xx) +
  xlab("Species") +
  ylab("Emergence as proportion of age at sexual maturation") +
  ggtitle("(E) Sociosexual Play")+
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.position = "none")
windows()
plot2e

fig2 <- plot2a+plot2b+plot2c+plot2d+plot2e+plot_layout(nrow=2)
windows()
fig2

# 29 Sex differences - absolute ####

# data checks

sort(unique(m2d$population))
sort(unique(m2d$sex))
length(unique(m2d$id_code))
length(unique(m2d$milestone))
sort(unique(m2d$milestone))
m2d$complexity <- as.factor(as.numeric(m2d$complexity))
table(m2d$sex,m2d$milestone_category)
m2d$group_size <- as.numeric(as.character(m2d$group_size))
range(m2d$group_size)

hist(m2d$age)
hist(m2d$age_adjust)
m2d$zsize = (m2d$group_size - mean(m2d$group_size))/sd(m2d$group_size)
range(m2d$zsize)

mod = lm(age ~ zsize + 
           milestone_category + 
           complexity +
           sex, 
         data = m2d)
vif(mod)

# vifs a bit high with both in so let's just use complexity

mod = lm(age ~ zsize + 
           complexity +
           sex, 
         data = m2d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 11"
xx$variable <- rownames(xx)
xx <- xx[,c(2,3,1)]

vif_table <- rbind(vif_table,xx)

msexprior <- get_prior(age ~ zsize + 
                       sex*complexity +
                       (1+zsize|id_code)+
                       (1+zsize|milestone),
                     data = m2d, family = lognormal())

msexprior = c(prior(lognormal(0.5,0.5), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*complexity +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=msexprior)

start_time <- Sys.time()
msex <- brm(age ~ zsize + 
            sex*complexity +
            (1+zsize|id_code)+
            (1+zsize|milestone),
          data = m2d, 
          family = lognormal(),
          prior=msexprior,
          chains = 3, 
          cores = 3,
          backend = "cmdstanr",
          stan_model_args = list(stanc_options=list("O1")),
          iter = 5000, 
          warmup = 2500, 
          thin = 2,
          control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(msex)
summary(msex)
pppsex <- pp_check(msex)
windows()
pppsex
pppsex <- pppsex+labs(title="Model 21: Weakly regularising priors")

modsex_table <- as.data.frame(round(posterior_summary(msex, probs = c(0.05, 0.95)),3))
modsex_table <- modsex_table[c(1:7),]
modsex_table$variable <- rownames(modsex_table)
rownames(modsex_table) <- NULL
modsex_valid <- as.data.frame(summary(msex)$fixed)
modsex_valid$model <- "Model 22a"
modsex_valid <- modsex_valid[,c(6:8)]

msex_coef <- mcmc_plot(msex, variable = "^b_", regex=TRUE)
msex_coef <- msex_coef+labs(title="Model 21: Weakly regularising priors")

# prior sensitivity checks

msexprior2 = c(prior(lognormal(0.5,100), class=Intercept),
              prior(normal(0,100), class=b),
              prior(exponential(1), class=sd),
              prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*complexity +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=msexprior2)

start_time <- Sys.time()
msexuni <- brm(age ~ zsize + 
              sex*complexity +
              (1+zsize|id_code)+
              (1+zsize|milestone),
            data = m2d, 
            family = lognormal(),
            prior=msexprior2,
            chains = 3, 
            cores = 3,
            backend = "cmdstanr",
            stan_model_args = list(stanc_options=list("O1")),
            iter = 5000, 
            warmup = 2500, 
            thin = 2,
            control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(msexuni)
summary(msexuni)
pppsexuni <- pp_check(msexuni)
windows()
pppsexuni
pppsexuni <- pppsexuni+labs(title="Model 21: Uniform priors")

msexuni_coef <- mcmc_plot(msexuni, variable = "^b_", regex=TRUE)
msexuni_coef <- msexuni_coef+labs(title="Model 21: Uniform priors")

# strong

msexprior3 = c(prior(lognormal(0.25,0.25), class=Intercept),
               prior(normal(0,0.5), class=b),
               prior(exponential(1), class=sd),
               prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*complexity +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=msexprior3)

start_time <- Sys.time()
msexstr <- brm(age ~ zsize + 
                 sex*complexity +
                 (1+zsize|id_code)+
                 (1+zsize|milestone),
               data = m2d, 
               family = lognormal(),
               prior=msexprior3,
               chains = 3, 
               cores = 3,
               backend = "cmdstanr",
               stan_model_args = list(stanc_options=list("O1")),
               iter = 5000, 
               warmup = 2500, 
               thin = 2,
               control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(msexstr)
summary(msexstr)
pppsexstr <- pp_check(msexstr)
windows()
pppsexstr
pppsexstr <- pppsexstr+labs(title="Model 21: Strongly regularising priors")

msexstr_coef <- mcmc_plot(msexstr, variable = "^b_", regex=TRUE)
msexstr_coef <- msexuni_coef+labs(title="Model 21: Strongly regularising priors")

msex_fixed_comparisons <- msex_coef+
  msexuni_coef+
  msexstr_coef
windows()
msex_fixed_comparisons

msex_ppcs <- pppsex+pppsexuni+pppsexstr
windows()
msex_ppcs

# model passes all checks

# 30 Sex differences - relative ####

msexpriorb <- get_prior(age_adjust ~ zsize + 
                         sex*complexity +
                         (1+zsize|id_code)+
                         (1+zsize|milestone),
                       data = m2d, family = Beta())

msexpriorb = c(prior(normal(0,1), class=b),
              prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*complexity +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=msexpriorb)

start_time <- Sys.time()
msexb <- brm(age_adjust ~ zsize + 
              sex*complexity +
              (1+zsize|id_code)+
              (1+zsize|milestone),
            data = m2d, 
            family = Beta(),
            prior=msexpriorb,
            chains = 3, 
            cores = 3,
            backend = "cmdstanr",
            stan_model_args = list(stanc_options=list("O1")),
            iter = 5000, 
            warmup = 2500, 
            thin = 2,
            control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(msexb)
summary(msexb)
pppsexb <- pp_check(msexb)
windows()
pppsexb
pppsexb <- pppsexb+labs(title="Model 22: Weakly regularising priors")

modsexb_table <- as.data.frame(round(posterior_summary(msexb, probs = c(0.05, 0.95)),3))
modsexb_table <- modsexb_table[c(1:7),]
modsexb_table$variable <- rownames(modsexb_table)
rownames(modsexb_table) <- NULL
modsexb_valid <- as.data.frame(summary(msexb)$fixed)
modsexb_valid$model <- "Model 22b"
modsexb_valid <- modsexb_valid[,c(6:8)]

msexb_coef <- mcmc_plot(msexb, variable = "^b_", regex=TRUE)
msexb_coef <- msexb_coef+labs(title="Model 22: Weakly regularising priors")

# prior sensitivity checks

msexprior2b = c(prior(normal(0,100), class=b),
               prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*complexity +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=msexprior2b)

start_time <- Sys.time()
msexunib <- brm(age_adjust ~ zsize + 
                 sex*complexity +
                 (1+zsize|id_code)+
                 (1+zsize|milestone),
               data = m2d, 
               family = Beta(),
               prior=msexprior2b,
               chains = 3, 
               cores = 3,
               backend = "cmdstanr",
               stan_model_args = list(stanc_options=list("O1")),
               iter = 5000, 
               warmup = 2500, 
               thin = 2,
               control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(msexunib)
summary(msexunib)
pppsexunib <- pp_check(msexunib)
windows()
pppsexunib
pppsexunib <- pppsexunib+labs(title="Model 22: Uniform priors")

msexunib_coef <- mcmc_plot(msexunib, variable = "^b_", regex=TRUE)
msexunib_coef <- msexunib_coef+labs(title="Model 22: Uniform priors")

# strong

msexprior3b = c(prior(normal(0,0.5), class=b),
               prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*complexity +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=msexprior3)

start_time <- Sys.time()
msexstrb <- brm(age_adjust ~ zsize + 
                 sex*complexity +
                 (1+zsize|id_code)+
                 (1+zsize|milestone),
               data = m2d, 
               family = Beta(),
               prior=msexprior3b,
               chains = 3, 
               cores = 3,
               backend = "cmdstanr",
               stan_model_args = list(stanc_options=list("O1")),
               iter = 5000, 
               warmup = 2500, 
               thin = 2,
               control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(msexstrb)
summary(msexstrb)
pppsexstrb <- pp_check(msexstrb)
windows()
pppsexstrb
pppsexstrb <- pppsexstrb+labs(title="Model 22: Strongly regularising priors")

msexstrb_coef <- mcmc_plot(msexstrb, variable = "^b_", regex=TRUE)
msexstrb_coef <- msexunib_coef+labs(title="Model 22: Strongly regularising priors")

msexb_fixed_comparisons <- msexb_coef+
  msexunib_coef+
  msexstrb_coef
windows()
msexb_fixed_comparisons

msexb_ppcs <- pppsexb+pppsexunib+pppsexstrb
windows()
msexb_ppcs

# model passes all checks

# 31 Summarise Models 22 and 23 ####

# model 22

post <- posterior_samples(msex)
names(post)

round(sum(post$`b_zsize` < 0) / 
        length(post$`b_zsize`),3)
#0.550
round(sum(post$`b_sexmale` < 0) / 
        length(post$`b_sexmale`),3)
#0.851
round(sum(post$`b_complexity1` < 0) / 
        length(post$`b_complexity1`),3)
#0.995
round(sum(post$`b_complexity2` < 0) / 
        length(post$`b_complexity2`),3)
#0.881

round(sum(post$`b_sexmale:complexity1` > 0) / 
        length(post$`b_sexmale:complexity1`),3)
#0.517

round(sum(post$`b_sexmale:complexity2` > 0) / 
        length(post$`b_sexmale:complexity2`),3)
#0.979

# model 23
post <- posterior_samples(msexb)
names(post)

round(sum(post$`b_zsize` < 0) / 
        length(post$`b_zsize`),3)
#0.550
round(sum(post$`b_sexmale` < 0) / 
        length(post$`b_sexmale`),3)
#0.851
round(sum(post$`b_complexity1` < 0) / 
        length(post$`b_complexity1`),3)
#0.995
round(sum(post$`b_complexity2` < 0) / 
        length(post$`b_complexity2`),3)
#0.881

round(sum(post$`b_sexmale:complexity1` > 0) / 
        length(post$`b_sexmale:complexity1`),3)
#0.917

round(sum(post$`b_sexmale:complexity2` > 0) / 
        length(post$`b_sexmale:complexity2`),3)
#0.996

# 32 Figure 3 ####

fit_msex <- 
  fitted(msexb,probs=c(0.05,0.95)) %>%
  as_tibble() %>%
  bind_cols(m2d)

pred_vals = round(predict(msexb, summary = FALSE), 3)
pred_vals = as.data.frame(t(pred_vals))

# calculate mean of predictions

pred_vals$predict = rowMeans(pred_vals)

fit_msex$predict = pred_vals$predict
range(fit_msex$predict)

# generate a list of 100 random numbers between 1 and 15000

xx = floor(runif(100, 1, 3750))

# extract those from our predictions object

zzz = pred_vals[,xx]

colnames(zzz) = c("pred1", "pred2", "pred3", "pred4", "pred5", "pred6", "pred7", "pred8", "pred9", "pred10",
                  "pred11", "pred12", "pred13", "pred14", "pred15", "pred16", "pred17", "pred18", "pred19", "pred20",
                  "pred21", "pred22", "pred23", "pred24", "pred25", "pred26", "pred27", "pred28", "pred29", "pred30",
                  "pred31", "pred32", "pred33", "pred34", "pred35", "pred36", "pred37", "pred38", "pred39", "pred40",
                  "pred41", "pred42", "pred43", "pred44", "pred45", "pred46", "pred47", "pred48", "pred49", "pred50",
                  "pred51", "pred52", "pred53", "pred54", "pred55", "pred56", "pred57", "pred58", "pred59", "pred60",
                  "pred61", "pred62", "pred63", "pred64", "pred65", "pred66", "pred67", "pred68", "pred69", "pred70",
                  "pred71", "pred72", "pred73", "pred74", "pred75", "pred76", "pred77", "pred78", "pred79", "pred80",
                  "pred81", "pred82", "pred83", "pred84", "pred85", "pred86", "pred87", "pred88", "pred89", "pred90",
                  "pred91", "pred92", "pred93", "pred94", "pred95", "pred96", "pred97", "pred98", "pred99", "pred100")

predictions <- rowMeans(zzz)
predictions <- as.data.frame(predictions)

fit_msex <- as.data.frame(cbind(fit_msex,zzz,predictions))


xx = aggregate(fit_msex$Estimate, by = list(fit_msex$complexity,fit_msex$sex), FUN=mean)
colnames(xx) = c("complexity", "sex","mean_mean")
xx$mean_mean

fit_msex$mean_mean = NA

for(i in 1:nrow(fit_msex)){
  yy = which(xx$complexity == fit_msex$complexity[i]&
               xx$sex == fit_msex$sex[i])
  if(length(yy>0)){fit_msex$mean_mean[i] = xx$mean_mean[yy]}}

xx = aggregate(fit_msex$Estimate, 
               by = list(fit_msex$complexity,fit_msex$sex), 
               FUN=sd)
colnames(xx) = c("complexity", "sex","sd")


fit_msex$sd = NA

for(i in 1:nrow(fit_msex)){
  yy = which(xx$complexity == fit_msex$complexity[i]&
               xx$sex == fit_msex$sex[i])
  if(length(yy>0)){fit_msex$sd[i] = xx$sd[yy]}}

plot3_d <- fit_msex

sort(unique(plot3_d$sex))
plot3_d = plot3_d %>%
  mutate(sex = gsub("^female$", "Female", sex),
         sex = gsub("^male$", "Male", sex))

xx <- met.brewer(name="Troy", n=2, type="discrete")
plot3 <- ggplot(plot3_d, aes(x = factor(complexity), 
                             y = Estimate, 
                             fill = factor(sex))) +
  geom_violin(position=position_dodge(width=0.9),
              width=1.1,
              trim=FALSE) +
  geom_boxplot(data=plot3_d,
               aes(x=factor(complexity),
                   y=Estimate,
                   col=factor(sex)),
               fill="white",
               position=position_dodge(width=0.9),
               width=0.075,
               outlier.shape = NA)+
  geom_boxplot(position = position_dodge(width = 0.9), 
               alpha=0, width=0.075)+
  scale_fill_manual(values = xx) +
  xlab("Complexity") +
  ylab("Emergence as proportion of age at sexual maturation") +
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
plot3


# 33 Sex differences - with category + absolute ages ####

m23prior <- get_prior(age ~ zsize + 
                         sex*milestone_category +
                         (1+zsize|id_code)+
                         (1+zsize|milestone),
                       data = m2d, family = lognormal())

m23prior = c(prior(lognormal(0.5,0.5), class=Intercept),
              prior(normal(0,1), class=b),
              prior(exponential(1), class=sd),
              prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=msexprior)

start_time <- Sys.time()
m23 <- brm(age ~ zsize + 
              sex*milestone_category +
              (1+zsize|id_code)+
              (1+zsize|milestone),
            data = m2d, 
            family = lognormal(),
            prior=msexprior,
            chains = 3, 
            cores = 3,
            backend = "cmdstanr",
            stan_model_args = list(stanc_options=list("O1")),
            iter = 5000, 
            warmup = 2500, 
            thin = 2,
            control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(m23)
summary(m23)
ppp23 <- pp_check(m23)
windows()
ppp23
ppp23 <- ppp23+labs(title="Model 23: Weakly regularising priors")

mod23_table <- as.data.frame(round(posterior_summary(m23, probs = c(0.05, 0.95)),3))
mod23_table <- mod23_table[c(1:11),]
mod23_table$variable <- rownames(mod23_table)
rownames(mod23_table) <- NULL
mod23_valid <- as.data.frame(summary(m23)$fixed)
mod23_valid$model <- "Model 23"
mod23_valid <- mod23_valid[,c(6:8)]

m23_coef <- mcmc_plot(m23, variable = "^b_", regex=TRUE)
m23_coef <- m23_coef+labs(title="Model 23: Weakly regularising priors")

# prior sensitivity checks

m23prior2 = c(prior(lognormal(0.5,100), class=Intercept),
               prior(normal(0,100), class=b),
               prior(exponential(1), class=sd),
               prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=m23prior2)

start_time <- Sys.time()
m23uni <- brm(age ~ zsize + 
                 sex*milestone_category +
                 (1+zsize|id_code)+
                 (1+zsize|milestone),
               data = m2d, 
               family = lognormal(),
               prior=m23prior2,
               chains = 3, 
               cores = 3,
               backend = "cmdstanr",
               stan_model_args = list(stanc_options=list("O1")),
               iter = 5000, 
               warmup = 2500, 
               thin = 2,
               control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(m23uni)
summary(m23uni)
ppm23uni <- pp_check(m23uni)
windows()
ppm23uni
ppm23uni <- ppm23uni+labs(title="Model 23: Uniform priors")

m23uni_coef <- mcmc_plot(m23uni, variable = "^b_", regex=TRUE)
m23uni_coef <- m23uni_coef+labs(title="Model 23: Uniform priors")

# strong

m23prior3 = c(prior(lognormal(0.25,0.25), class=Intercept),
               prior(normal(0,0.5), class=b),
               prior(exponential(1), class=sd),
               prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=m23prior3)

start_time <- Sys.time()
m23str <- brm(age ~ zsize + 
                 sex*milestone_category +
                 (1+zsize|id_code)+
                 (1+zsize|milestone),
               data = m2d, 
               family = lognormal(),
               prior=m23prior3,
               chains = 3, 
               cores = 3,
               backend = "cmdstanr",
               stan_model_args = list(stanc_options=list("O1")),
               iter = 5000, 
               warmup = 2500, 
               thin = 2,
               control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(m23str)
summary(m23str)
ppm23str <- pp_check(m23str)
windows()
ppm23str
ppm23str <- ppm23str+labs(title="Model 23: Strongly regularising priors")

m23str_coef <- mcmc_plot(m23str, variable = "^b_", regex=TRUE)
m23str_coef <- m23str_coef+labs(title="Model 23: Strongly regularising priors")

m23_fixed_comparisons <- m23_coef+
  m23uni_coef+
  m23str_coef
windows()
m23_fixed_comparisons

m23_ppcs <- ppp23+ppm23uni+ppm23str
windows()
m23_ppcs

# model passes all checks

# 34 Sex differences - relative ####

m24prior <- get_prior(age_adjust ~ zsize + 
                        sex*milestone_category +
                        (1+zsize|id_code)+
                        (1+zsize|milestone),
                      data = m2d, family = Beta())

m24prior = c(prior(normal(0,1), class=b),
             prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=m24prior)

start_time <- Sys.time()
m24 <- brm(age_adjust ~ zsize + 
             sex*milestone_category +
             (1+zsize|id_code)+
             (1+zsize|milestone),
           data = m2d, 
           family = Beta(),
           prior=m24prior,
           chains = 3, 
           cores = 3,
           backend = "cmdstanr",
           stan_model_args = list(stanc_options=list("O1")),
           iter = 5000, 
           warmup = 2500, 
           thin = 2,
           control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(m24)
summary(m24)
ppp24 <- pp_check(m24)
windows()
ppp24
ppp24 <- ppp24+labs(title="Model 24: Weakly regularising priors")

mod24_table <- as.data.frame(round(posterior_summary(m24, probs = c(0.05, 0.95)),3))
mod24_table <- mod24_table[c(1:11),]
mod24_table$variable <- rownames(mod24_table)
rownames(mod24_table) <- NULL
mod24_valid <- as.data.frame(summary(m24)$fixed)
mod24_valid$model <- "Model 24"
mod24_valid <- mod24_valid[,c(6:8)]

m24_coef <- mcmc_plot(m24, variable = "^b_", regex=TRUE)
m24_coef <- m24_coef+labs(title="Model 24: Weakly regularising priors")

# prior sensitivity checks

m24prior2 = c(prior(normal(0,100), class=b),
              prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=m24prior2)

start_time <- Sys.time()
m24uni <- brm(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=m24prior2,
              chains = 3, 
              cores = 3,
              backend = "cmdstanr",
              stan_model_args = list(stanc_options=list("O1")),
              iter = 5000, 
              warmup = 2500, 
              thin = 2,
              control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(m24uni)
summary(m24uni)
ppm24uni <- pp_check(m24uni)
windows()
ppm24uni
ppm24uni <- ppm24uni+labs(title="Model 24: Uniform priors")

m24uni_coef <- mcmc_plot(m24uni, variable = "^b_", regex=TRUE)
m24uni_coef <- m24uni_coef+labs(title="Model 24: Uniform priors")

# strong

m24prior3 = c(prior(normal(0,0.5), class=b),
              prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=m24prior3)

start_time <- Sys.time()
m24str <- brm(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=m24prior3,
              chains = 3, 
              cores = 3,
              backend = "cmdstanr",
              stan_model_args = list(stanc_options=list("O1")),
              iter = 5000, 
              warmup = 2500, 
              thin = 2,
              control = list(adapt_delta = 0.975))
end_time <- Sys.time()
mod_time <- end_time - start_time
beep(8)
mod_time
#Time difference of 1.054398 mins

windows()
plot(m24str)
summary(m24str)
ppm24str <- pp_check(m24str)
windows()
ppm24str
ppm24str <- ppm24str+labs(title="Model 24: Strongly regularising priors")

m24str_coef <- mcmc_plot(m24str, variable = "^b_", regex=TRUE)
m24str_coef <- m24str_coef+labs(title="Model 24: Strongly regularising priors")

m24_fixed_comparisons <- m24_coef+
  m24uni_coef+
  m24str_coef
windows()
m24_fixed_comparisons

m24_ppcs <- ppp24+ppm24uni+ppm24str
windows()
m24_ppcs

# model passes all checks

# 35 Summarise Models  ####

# model 23

post <- posterior_samples(m23)
names(post)

round(sum(post$`b_zsize` < 0) / 
        length(post$`b_zsize`),3)
#0.550
round(sum(post$`b_sexmale` < 0) / 
        length(post$`b_sexmale`),3)
#0.851

round(sum(post$`b_milestone_categoryphysical` > 0) / 
        length(post$`b_milestone_categoryphysical`),3)
#0.910
round(sum(post$`b_milestone_categoryproximity` > 0) / 
        length(post$`b_milestone_categoryproximity`),3)
#0.997
round(sum(post$`b_milestone_categorysocial_affiliative` > 0) / 
        length(post$`b_milestone_categorysocial_affiliative`),3)
#0.899
round(sum(post$`b_milestone_categorysocial_aggressive` > 0) / 
        length(post$`b_milestone_categorysocial_aggressive`),3)
#0.968

round(sum(post$`b_sexmale:milestone_categoryphysical` < 0) / 
        length(post$`b_sexmale:milestone_categoryphysical`),3)
#0.508
round(sum(post$`b_sexmale:milestone_categoryproximity` > 0) / 
        length(post$`b_sexmale:milestone_categoryproximity`),3)
#0.526
round(sum(post$`b_sexmale:milestone_categorysocial_affiliative` > 0) / 
        length(post$`b_sexmale:milestone_categorysocial_affiliative`),3)
#0.948

round(sum(post$`b_sexmale:milestone_categorysocial_aggressive` > 0) / 
        length(post$`b_sexmale:milestone_categorysocial_aggressive`),3)
#0.982


# relative age model

post <- posterior_samples(m24)
names(post)

round(sum(post$`b_zsize` < 0) / 
        length(post$`b_zsize`),3)
#0.533
round(sum(post$`b_sexmale` < 0) / 
        length(post$`b_sexmale`),3)
#0.998

round(sum(post$`b_milestone_categoryphysical` > 0) / 
        length(post$`b_milestone_categoryphysical`),3)
#0.929
round(sum(post$`b_milestone_categoryproximity` > 0) / 
        length(post$`b_milestone_categoryproximity`),3)
#0.997
round(sum(post$`b_milestone_categorysocial_affiliative` > 0) / 
        length(post$`b_milestone_categorysocial_affiliative`),3)
#0.958
round(sum(post$`b_milestone_categorysocial_aggressive` > 0) / 
        length(post$`b_milestone_categorysocial_aggressive`),3)
#0.933

round(sum(post$`b_sexmale:milestone_categoryphysical` < 0) / 
        length(post$`b_sexmale:milestone_categoryphysical`),3)
#0.537
round(sum(post$`b_sexmale:milestone_categoryproximity` < 0) / 
        length(post$`b_sexmale:milestone_categoryproximity`),3)
#0.855
round(sum(post$`b_sexmale:milestone_categorysocial_affiliative` > 0) / 
        length(post$`b_sexmale:milestone_categorysocial_affiliative`),3)
#0.808

round(sum(post$`b_sexmale:milestone_categorysocial_aggressive` > 0) / 
        length(post$`b_sexmale:milestone_categorysocial_aggressive`),3)
#0.897

# 32 Figure 4 ####

fit_m23 <- 
  fitted(m23,probs=c(0.05,0.95)) %>%
  as_tibble() %>%
  bind_cols(m2d)

pred_vals = round(predict(m23, summary = FALSE), 3)
pred_vals = as.data.frame(t(pred_vals))

# calculate mean of predictions

pred_vals$predict = rowMeans(pred_vals)

fit_m23$predict = pred_vals$predict
range(fit_m23$predict)

# generate a list of 100 random numbers between 1 and 15000

xx = floor(runif(100, 1, 3750))

# extract those from our predictions object

zzz = pred_vals[,xx]

colnames(zzz) = c("pred1", "pred2", "pred3", "pred4", "pred5", "pred6", "pred7", "pred8", "pred9", "pred10",
                  "pred11", "pred12", "pred13", "pred14", "pred15", "pred16", "pred17", "pred18", "pred19", "pred20",
                  "pred21", "pred22", "pred23", "pred24", "pred25", "pred26", "pred27", "pred28", "pred29", "pred30",
                  "pred31", "pred32", "pred33", "pred34", "pred35", "pred36", "pred37", "pred38", "pred39", "pred40",
                  "pred41", "pred42", "pred43", "pred44", "pred45", "pred46", "pred47", "pred48", "pred49", "pred50",
                  "pred51", "pred52", "pred53", "pred54", "pred55", "pred56", "pred57", "pred58", "pred59", "pred60",
                  "pred61", "pred62", "pred63", "pred64", "pred65", "pred66", "pred67", "pred68", "pred69", "pred70",
                  "pred71", "pred72", "pred73", "pred74", "pred75", "pred76", "pred77", "pred78", "pred79", "pred80",
                  "pred81", "pred82", "pred83", "pred84", "pred85", "pred86", "pred87", "pred88", "pred89", "pred90",
                  "pred91", "pred92", "pred93", "pred94", "pred95", "pred96", "pred97", "pred98", "pred99", "pred100")

predictions <- rowMeans(zzz)
predictions <- as.data.frame(predictions)

fit_m23 <- as.data.frame(cbind(fit_m23,zzz,predictions))


xx = aggregate(fit_m23$Estimate, by = list(fit_m23$milestone_category,fit_m23$sex), FUN=mean)
colnames(xx) = c("milestone_category", "sex","mean_mean")
xx$mean_mean

fit_m23$mean_mean = NA

for(i in 1:nrow(fit_m23)){
  yy = which(xx$milestone_category == fit_m23$milestone_category[i]&
               xx$sex == fit_m23$sex[i])
  if(length(yy>0)){fit_m23$mean_mean[i] = xx$mean_mean[yy]}}

xx = aggregate(fit_m23$Estimate, 
               by = list(fit_m23$milestone_category,fit_m23$sex), 
               FUN=sd)
colnames(xx) = c("milestone_category", "sex","sd")


fit_m23$sd = NA

for(i in 1:nrow(fit_m23)){
  yy = which(xx$milestone_category == fit_m23$milestone_category[i]&
               xx$sex == fit_m23$sex[i])
  if(length(yy>0)){fit_m23$sd[i] = xx$sd[yy]}}

plot4_d <- fit_m23

sort(unique(plot4_d$sex))
plot4_d = plot4_d %>%
  mutate(sex = gsub("^female$", "Female", sex),
         sex = gsub("^male$", "Male", sex))

sort(unique(plot4_d$milestone_category))
plot4_d = plot4_d %>%
  mutate(milestone_category = gsub("^Carry$", "Carrying", milestone_category),
         milestone_category = gsub("^physical$", "Non-social physicality", milestone_category),
         milestone_category = gsub("^proximity$", "Proximity", milestone_category),
         milestone_category = gsub("^social_affiliative$", "Social-affiliative", milestone_category),
         milestone_category = gsub("^social_aggressive$", "Social-agonistic", milestone_category))

xx <- met.brewer(name="Troy", n=2, type="discrete")
plot4 <- ggplot(plot4_d, aes(x = factor(milestone_category), 
                             y = Estimate, 
                             fill = factor(sex))) +
  geom_violin(position=position_dodge(width=0.9),
              width=1.1,
              trim=FALSE) +
  geom_boxplot(data=plot4_d,
               aes(x=factor(milestone_category),
                   y=Estimate,
                   col=factor(sex)),
               fill="white",
               position=position_dodge(width=0.9),
               width=0.075,
               outlier.shape = NA)+
  geom_boxplot(position = position_dodge(width = 0.9), 
               alpha=0, width=0.075)+
  scale_fill_manual(values = xx) +
  xlab("Milestone category") +
  ylab("Emergence as proportion of age at sexual maturation") +
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
plot4


# 33 Chimp-gorilla comparison plot ####

str(chimp_gorilla)

xx <- met.brewer(name="Pillement", n=3, type="discrete")

sort(unique(chimp_gorilla$Species))
chimp_gorilla$Species2 <- factor(chimp_gorilla$Species, 
                                 levels=c('Mountain gorilla',
                                          'Western gorilla',
                                          'Western chimpanzees'))
levels(chimp_gorilla$Species2)

plot3a <- ggplot(chimp_gorilla,
                aes(x = reorder(Milestone,-Mean),
                    y = Mean)) + 
  geom_pointrange(aes(ymin = Mean - sd,
                      ymax = Mean + sd,
                      color=Species2,
                      shape=Species2),
                  position = position_dodge(0.5), 
                  size=1) +
  #scale_colour_discrete(limits=c('Western gorilla', 'Western chimpanzees', 'Mountain gorilla'))+
  scale_colour_manual(values=xx)+
  coord_flip()+
  theme_minimal()+
  xlab("Milestone")+
  ylab("Mean age of emergence (days)")+
  ggtitle("(A)")+
  theme(axis.title=element_text(face="bold",size=20),
        axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=20),
        legend.position="none",
        strip.text=element_text(face="bold", size=20))
#plot3 <- plot3+scale_colour_discrete(breaks=c('Western gorilla', 'Western chimpanzees', 'Mountain gorilla'))
#windows()  
#plot3a


plot3b <- ggplot(chimp_gorilla,
                 aes(x = reorder(Milestone,-Mean_adj),
                     y = Mean_adj)) + 
  geom_pointrange(aes(ymin = Mean_adj - sd_adj,
                      ymax = Mean_adj + sd_adj,
                      color=Species2,
                      shape=Species2),
                  position = position_dodge(0.5), 
                  size=1) +
  #scale_colour_discrete(limits=c('Western gorilla', 'Western chimpanzees', 'Mountain gorilla'))+
  scale_colour_manual(values=xx)+
  coord_flip()+
  theme_minimal()+
  xlab("Milestone")+
  ylab("Mean age of emergence as proportion of development")+
  ggtitle("(B)")+
  theme(axis.title.x=element_text(face="bold",size=20),
        axis.text=element_text(size=15),
        axis.title.y=element_blank(),
        plot.title=element_text(face="bold",size=20),
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        strip.text=element_text(face="bold", size=20))
#plot3 <- plot3+scale_colour_discrete(breaks=c('Western gorilla', 'Western chimpanzees', 'Mountain gorilla'))
#windows()  
#plot3b

plot4 <- plot3a+plot3b
windows()  
plot4

# 34 Check supplementary materials info ####

save.image("working.Rdata")
# age plot

gorilla_age_plot$duration <- gorilla_age_plot$age_last_datapoint -
  gorilla_age_plot$age_first
gorilla_age_plot$id_code <- with(gorilla_age_plot,
                                 factor(id_code, levels = id_code[order(-duration)]))
ggplot(gorilla_age_plot,
       aes(x = age_first,
           xend = age_last_datapoint,
           y = id_code,
           yend = id_code,
           colour = Species,
           linetype = Sex)) +
  scale_color_manual(values = xx) +
  geom_segment(linewidth = 1) +
  labs(x = "Age",
       y = "Individual ID",
       colour = "Species",
       linetype = "Sex") +
  #facet_wrap(~Sex)+
  theme_minimal()



vif_table
ppp1
windows()
m1_ppcs
windows()
m1b_ppcs
windows()
m2_pppcs
windows()
m2a_pppcs

# create prior table


prior_table <- as.data.frame(prior_summary(m1))
prior_table$model <- "Model 1-10: Weakly regularizing priors"

xx <- as.data.frame(prior_summary(m1str))
xx$model <- "Model 1-10: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1uni))
xx$model <- "Model 1-10: Uniform priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1b))
xx$model <- "Model 11-20: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1bstr))
xx$model <- "Model 11-20: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1buni))
xx$model <- "Model 11-20: Uniform priors"

prior_table <- rbind(prior_table,xx)

# now the sex difference models

xx <- as.data.frame(prior_summary(msex))
xx$model <- "Model 21: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(msexstr))
xx$model <- "Model 21: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(msexuni))
xx$model <- "Model 21: Uniform priors"

prior_table <- rbind(prior_table,xx)

#

xx <- as.data.frame(prior_summary(msexb))
xx$model <- "Model 22: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(msexstrb))
xx$model <- "Model 22: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(msexunib))
xx$model <- "Model 22: Uniform priors"

prior_table <- rbind(prior_table,xx)

#

xx <- as.data.frame(prior_summary(m23))
xx$model <- "Model 23: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m23str))
xx$model <- "Model 23: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m23uni))
xx$model <- "Model 23: Uniform priors"

prior_table <- rbind(prior_table,xx)

#

xx <- as.data.frame(prior_summary(m24))
xx$model <- "Model 24: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m24str))
xx$model <- "Model 24: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m24uni))
xx$model <- "Model 24: Uniform priors"

prior_table <- rbind(prior_table,xx)

# Save ###

save.image("revision_final.Rdata")


























