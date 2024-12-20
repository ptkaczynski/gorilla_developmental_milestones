# Gorilla milestones analysis ####

# 1 Description ####

# In model 1, we test whether: 
# (a) different categories of developmental milestones emerged at different ages of ontogeny within each species and  
# (b) whether each category of developmental milestones emerged at different ages in the two gorilla species 

# Model 1 needs to include only milestones that occur in each species.
# Model 1 is tested with object m1d
# We also run model 1 using ages adjusted for age at first reproduction for females, and at physical maturity for males

# In model 2, we test whether: there are sex differences in the timing of milestones using the entire list of milestones in Bwindi

# Model 1 needs to include only milestones that occur in each species.
# Model 1 is tested with object m1d
# We also run model 1 using ages adjusted for age at first reproduction for females, and at physical maturity for males

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

# 3 Model 1a ####

# do data checks

sort(unique(m1d$sex))
length(unique(m1d$id_code))
length(unique(m1d$milestone))
sort(unique(m1d$milestone))
table(m1d$population,m1d$milestone_category)
m1d$milestone_category <- factor(m1d$milestone_category, levels = c("Carry","proximity","physical","social_affiliative","social_aggressive"))
m1d$group_size <- as.numeric(as.character(m1d$group_size))
range(m1d$group_size)
hist(m1d$age)
hist(m1d$age_adjust)
m1d$zsize = (m1d$group_size - mean(m1d$group_size))/sd(m1d$group_size)
range(m1d$zsize)

covees = as.data.frame(
  cbind(
    m1d$zsize, 
    m1d$population, 
    m1d$milestone_category, 
    m1d$sex))
colnames(covees) = c("group_size",
                     "population",
                     "category",
                     "sex")

mod = lm(age ~ zsize + 
           population + 
           milestone_category +
           sex, 
         data = m1d)
vif(mod)

vif_table <- as.data.frame(vif(mod))
vif_table$model <- "Model 1"
vif_table$variable <- rownames(vif_table)
vif_table <- vif_table[,c(4,5,1,2,3)]

rm(mod)
windows()
ggpairs(covees)# nothing to look at as all categorical

# check priors

# can try a lognormal model
# let's check prior

m1prior <- get_prior(age ~ zsize + population*milestone_category +
                       sex*milestone_category +
                       (1+zsize|id_code)+
                       (1+zsize+sex|milestone),
                     data = m1d, family = lognormal())

m1prior = c(prior(lognormal(0.5,0.5), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population*milestone_category +
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = lognormal(),
              prior=m1prior)

start_time <- Sys.time()
prior_check <- brm(
  age ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = lognormal(),
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
  age ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = lognormal(),
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
ppp1 <- ppp1+labs(title="Model 1a: Weakly regularising priors")

mod1_table <- as.data.frame(round(posterior_summary(m1, probs = c(0.05, 0.95)),3))
mod1_table <- mod1_table[c(1:16),]
mod1_table$variable <- rownames(mod1_table)
rownames(mod1_table) <- NULL
mod1_valid <- as.data.frame(summary(m1)$fixed)
mod1_valid$model <- "Model 1"
mod1_valid <- mod1_valid[,c(6:8)]

m1_coef <- mcmc_plot(m1, variable = "^b_", regex=TRUE)
m1_coef <- m1_coef+labs(title="Model 1a: Weakly regularising priors")

# prior sensitivity checks

m1prior_uniform = c(prior(lognormal(0.5,100), class=Intercept),
            prior(normal(0,100), class=b),
            prior(exponential(1), class=sd),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population*milestone_category +
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = lognormal(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m1uni <- brm(
  age ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = lognormal(),
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
ppp1uni <- ppp1uni+labs(title="Model 1a: Uniform priors")

m1uni_coef <- mcmc_plot(m1uni, variable = "^b_", regex=TRUE)
m1uni_coef <- m1uni_coef+labs(title="Model 1a: Uniform priors")

##

m1prior_strong = c(prior(lognormal(0.5,0.25), class=Intercept),
                    prior(normal(0,0.5), class=b),
                    prior(exponential(1), class=sd),
                    prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + population*milestone_category +
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = lognormal(),
              prior=m1prior_strong)

start_time <- Sys.time()
m1str <- brm(
  age ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = lognormal(),
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
ppp1str <- ppp1str+labs(title="Model 1a: Strong priors")


m1str_coef <- mcmc_plot(m1str, variable = "^b_", regex=TRUE)
m1str_coef <- m1str_coef+labs(title="Model 1a: Strong priors")

m1_fixed_comparisons <- m1_coef+
  m1uni_coef+
  m1str_coef
windows()
m1_fixed_comparisons

m1_ppcs <- ppp1+ppp1uni+ppp1str
windows()
m1_ppcs

# model passes all checks

# 4 Model 1b ####

# do data checks

m1aprior <- get_prior(age_adjust ~ zsize + population*milestone_category +
                       sex*milestone_category +
                       (1+zsize|id_code)+
                       (1+zsize+sex|milestone),
                     data = m1d, family = Beta())

m1aprior = c(prior(normal(0,1.1), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + population*milestone_category +
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = Beta(),
              prior=m1aprior)

start_time <- Sys.time()
prior_check <- brm(
  age_adjust ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = Beta(),
  prior=m1aprior,
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
m1a <- brm(
  age_adjust ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = Beta(),
  prior=m1aprior,
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
plot(m1a)
summary(m1a)
windows()
pp_check(m1a) # fine
ppp1a <- pp_check(m1a)
ppp1a <- ppp1a+labs(title="Model 1b: Weakly regularising priors")

mod1a_table <- as.data.frame(round(posterior_summary(m1a, probs = c(0.05, 0.95)),3))
mod1a_table <- mod1a_table[c(1:16),]
mod1a_table$variable <- rownames(mod1a_table)
rownames(mod1a_table) <- NULL
mod1a_valid <- as.data.frame(summary(m1a)$fixed)
mod1a_valid$model <- "Model 1a"
mod1a_valid <- mod1a_valid[,c(6:8)]

m1a_coef <- mcmc_plot(m1a, variable = "^b_", regex=TRUE)
m1a_coef <- m1_coef+labs(title="Model 1b: Weakly regularising priors")

# prior sensitivity checks

m1aprior_uni = c(prior(normal(0,101), class=Intercept),
             prior(normal(0,1), class=b),
             prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + population*milestone_category +
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = Beta(),
              prior=m1aprior_uni)

start_time <- Sys.time()
m1a_uni <- brm(
  age_adjust ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = Beta(),
  prior=m1aprior_uni,
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
plot(m1a_uni)
summary(m1a_uni)
windows()
pp_check(m1a_uni) # fine
ppp1auni <- pp_check(m1a_uni)
ppp1auni <- ppp1auni+labs(title="Model 1b: Uniform priors")

m1auni_coef <- mcmc_plot(m1a_uni, variable = "^b_", regex=TRUE)
m1auni_coef <- m1auni_coef+labs(title="Model 1b: Uniform priors")

##

m1aprior_str = c(prior(normal(0,0.5), class=Intercept),
                 prior(normal(0,1), class=b),
                 prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + population*milestone_category +
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = Beta(),
              prior=m1aprior_str)

start_time <- Sys.time()
m1a_str <- brm(
  age_adjust ~ zsize + population*milestone_category +
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = Beta(),
  prior=m1aprior_str,
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
plot(m1a_str)
summary(m1a_str)
windows()
pp_check(m1a_str) # fine
ppp1astr <- pp_check(m1a_str)
ppp1astr <- ppp1astr+labs(title="Model 1b: Strongly regularising priors")

m1astr_coef <- mcmc_plot(m1a_str, variable = "^b_", regex=TRUE)
m1astr_coef <- m1astr_coef+labs(title="Model 1b: Strongly regularising priors")


m1a_fixed_comparisons <- m1a_coef+
  m1auni_coef+
  m1astr_coef
windows()
m1a_fixed_comparisons

m1a_ppcs <- ppp1a+
  ppp1auni+
  ppp1astr
windows()
m1a_ppcs

# model passes all checks

# 5 Summarise Models 1a and 1b ####

mod1_table
mod1a_table

# proportion support for m1 coefficients

postm1 <- posterior_samples(m1)
names(postm1)

round(sum(postm1$`b_zsize` > 0) / 
        length(postm1$`b_zsize`),3)

round(sum(postm1$`b_populationLoango` < 0) / 
        length(postm1$`b_populationLoango`),3)

round(sum(postm1$`b_milestone_categoryproximity` > 0) / 
        length(postm1$`b_milestone_categoryproximity`),3)

round(sum(postm1$`b_milestone_categoryphysical` > 0) / 
        length(postm1$`b_milestone_categoryphysical`),3)

round(sum(postm1$`b_milestone_categorysocial_affiliative` > 0) / 
        length(postm1$`b_milestone_categorysocial_affiliative`),3)

round(sum(postm1$`b_milestone_categorysocial_aggressive` > 0) / 
        length(postm1$`b_milestone_categorysocial_aggressive`),3)

round(sum(postm1$`b_sexmale` < 0) / 
        length(postm1$`b_sexmale`),3)

round(sum(postm1$`b_populationLoango:milestone_categoryproximity` < 0) / 
        length(postm1$`b_populationLoango:milestone_categoryproximity`),3)

round(sum(postm1$`b_populationLoango:milestone_categoryphysical` > 0) / 
        length(postm1$`b_populationLoango:milestone_categoryphysical`),3)

round(sum(postm1$`b_populationLoango:milestone_categorysocial_affiliative` < 0) / 
        length(postm1$`b_populationLoango:milestone_categorysocial_affiliative`),3)

round(sum(postm1$`b_populationLoango:milestone_categorysocial_aggressive` > 0) / 
        length(postm1$`b_populationLoango:milestone_categorysocial_aggressive`),3)

round(sum(postm1$`b_milestone_categoryproximity:sexmale` < 0) / 
        length(postm1$`b_milestone_categoryproximity:sexmale`),3)

round(sum(postm1$`b_milestone_categoryphysical:sexmale` > 0) / 
        length(postm1$`b_milestone_categoryphysical:sexmale`),3)

round(sum(postm1$`b_milestone_categorysocial_affiliative:sexmale` > 0) / 
        length(postm1$`b_milestone_categorysocial_affiliative:sexmale`),3)

round(sum(postm1$`b_milestone_categorysocial_aggressive:sexmale` > 0) / 
        length(postm1$`b_milestone_categorysocial_aggressive:sexmale`),3)


postm1a <- posterior_samples(m1a)
names(postm1a)

round(sum(postm1a$`b_zsize` > 0) / 
        length(postm1a$`b_zsize`),3)

round(sum(postm1a$`b_populationLoango` < 0) / 
        length(postm1a$`b_populationLoango`),3)

round(sum(postm1a$`b_milestone_categoryproximity` > 0) / 
        length(postm1a$`b_milestone_categoryproximity`),3)

round(sum(postm1a$`b_milestone_categoryphysical` > 0) / 
        length(postm1a$`b_milestone_categoryphysical`),3)

round(sum(postm1a$`b_milestone_categorysocial_affiliative` > 0) / 
        length(postm1a$`b_milestone_categorysocial_affiliative`),3)

round(sum(postm1a$`b_milestone_categorysocial_aggressive` > 0) / 
        length(postm1a$`b_milestone_categorysocial_aggressive`),3)

round(sum(postm1a$`b_sexmale` < 0) / 
        length(postm1a$`b_sexmale`),3)

round(sum(postm1a$`b_populationLoango:milestone_categoryproximity` < 0) / 
        length(postm1a$`b_populationLoango:milestone_categoryproximity`),3)

round(sum(postm1a$`b_populationLoango:milestone_categoryphysical` < 0) / 
        length(postm1a$`b_populationLoango:milestone_categoryphysical`),3)

round(sum(postm1a$`b_populationLoango:milestone_categorysocial_affiliative` < 0) / 
        length(postm1a$`b_populationLoango:milestone_categorysocial_affiliative`),3)

round(sum(postm1a$`b_populationLoango:milestone_categorysocial_aggressive` < 0) / 
        length(postm1a$`b_populationLoango:milestone_categorysocial_aggressive`),3)

round(sum(postm1a$`b_milestone_categoryproximity:sexmale` < 0) / 
        length(postm1a$`b_milestone_categoryproximity:sexmale`),3)

round(sum(postm1a$`b_milestone_categoryphysical:sexmale` > 0) / 
        length(postm1a$`b_milestone_categoryphysical:sexmale`),3)

round(sum(postm1a$`b_milestone_categorysocial_affiliative:sexmale` > 0) / 
        length(postm1a$`b_milestone_categorysocial_affiliative:sexmale`),3)

round(sum(postm1a$`b_milestone_categorysocial_aggressive:sexmale` > 0) / 
        length(postm1a$`b_milestone_categorysocial_aggressive:sexmale`),3)

# 6 Figure 1 ####

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
  #facet_wrap(~species)+
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
plot1

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
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        axis.text.x=element_blank(),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
plot2

fig1 <- plot1a+plot1b
windows()
fig1

# 7 Model 1 reporting ####

fit_m1 <- 
  fitted(m1,probs=c(0.05,0.95)) %>%
  as_tibble() %>%
  bind_cols(m1d)

pred_vals = round(predict(m1, summary = FALSE), 3)
pred_vals = as.data.frame(t(pred_vals))

# calculate mean of predictions

pred_vals$predict = rowMeans(pred_vals)

fit_m1$predict = pred_vals$predict
range(fit_m1$predict)

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

fit_m1 <- as.data.frame(cbind(fit_m1,zzz))

# differences between physical and carry in absolute ages

sort(unique(fit_m1$milestone_category))

diff <- round((median(fit_m1$Estimate[fit_m1$milestone_category=="physical"])-
         median(fit_m1$Estimate[fit_m1$milestone_category=="Carry"])),3)

#[1] 0.431
diffl <- round((median(fit_m1$Q5[fit_m1$milestone_category=="physical"])-
         median(fit_m1$Q5[fit_m1$milestone_category=="Carry"])),3)
#[1] 0.365
diffu <- round((median(fit_m1$Q95[fit_m1$milestone_category=="physical"])-
         median(fit_m1$Q95[fit_m1$milestone_category=="Carry"])),3)
#[1] 0.529

# differences between physical and social affiliation in absolute ages

diff <- round((median(fit_m1$Estimate[fit_m1$milestone_category=="social_affiliative"])-
                 median(fit_m1$Estimate[fit_m1$milestone_category=="physical"])),3)

#[1] 1.165
diffl <- round((median(fit_m1$Q5[fit_m1$milestone_category=="social_affiliative"])-
                  median(fit_m1$Q5[fit_m1$milestone_category=="physical"])),3)
#[1] 0.947
diffu <- round((median(fit_m1$Q95[fit_m1$milestone_category=="social_affiliative"])-
                  median(fit_m1$Q95[fit_m1$milestone_category=="physical"])),3)
#[1] 1.325

# 8 Model 1a reporting ####

fit_m1a <- 
  fitted(m1a,probs=c(0.05,0.95)) %>%
  as_tibble() %>%
  bind_cols(m1d)

pred_vals = round(predict(m1a, summary = FALSE), 3)
pred_vals = as.data.frame(t(pred_vals))

# calculate mean of predictions

pred_vals$predict = rowMeans(pred_vals)

fit_m1a$predict = pred_vals$predict
range(fit_m1a$predict)

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

fit_m1a <- as.data.frame(cbind(fit_m1a,zzz))

# proximity differences
sort(unique(fit_m1a$milestone_category))
diff <- round((median(fit_m1a$Estimate[fit_m1a$population=="Bwindi"&
                                         fit_m1a$milestone_category=="proximity"])-
                 median(fit_m1a$Estimate[fit_m1a$population=="Loango"&
                                           fit_m1a$milestone_category=="proximity"])),3)
diffu <- round((median(fit_m1a$Q5[fit_m1a$population=="Bwindi"&
                                    fit_m1a$milestone_category=="proximity"])-
                  median(fit_m1a$Q5[fit_m1a$population=="Loango"&
                                      fit_m1a$milestone_category=="proximity"])),3)
diffl <- round((median(fit_m1a$Q95[fit_m1a$population=="Bwindi"&
                                     fit_m1a$milestone_category=="proximity"])-
                  median(fit_m1a$Q95[fit_m1a$population=="Loango"&
                                       fit_m1a$milestone_category=="proximity"])),3)
total <- round((median(fit_m1a$Estimate[fit_m1a$population=="Bwindi"&
                                          fit_m1a$milestone_category=="proximity"])+
                  median(fit_m1a$Estimate[fit_m1a$population=="Loango"&
                                            fit_m1a$milestone_category=="proximity"])),3)
totalu <- round((median(fit_m1a$Q5[fit_m1a$population=="Bwindi"&
                                     fit_m1a$milestone_category=="proximity"])+
                   median(fit_m1a$Q5[fit_m1a$population=="Loango"&
                                       fit_m1a$milestone_category=="proximity"])),3)
totall <- round((median(fit_m1a$Q95[fit_m1a$population=="Bwindi"&
                                      fit_m1a$milestone_category=="proximity"])+
                   median(fit_m1a$Q95[fit_m1a$population=="Loango"&
                                        fit_m1a$milestone_category=="proximity"])),3)
diff/(total/2)*100
# 60.36585

diffu/(totalu/2)*100
# 62.73063

diffl/(totall/2)*100
# 57.94872

#social-affiliative differences

diff <- round((median(fit_m1a$Estimate[fit_m1a$population=="Bwindi"&
                                        fit_m1a$milestone_category=="social_affiliative"])-
                 median(fit_m1a$Estimate[fit_m1a$population=="Loango"&
                                          fit_m1a$milestone_category=="social_affiliative"])),3)
diffu <- round((median(fit_m1a$Q5[fit_m1a$population=="Bwindi"&
                                   fit_m1a$milestone_category=="social_affiliative"])-
                  median(fit_m1a$Q5[fit_m1a$population=="Loango"&
                                     fit_m1a$milestone_category=="social_affiliative"])),3)
diffl <- round((median(fit_m1a$Q95[fit_m1a$population=="Bwindi"&
                                    fit_m1a$milestone_category=="social_affiliative"])-
                  median(fit_m1a$Q95[fit_m1a$population=="Loango"&
                                      fit_m1a$milestone_category=="social_affiliative"])),3)
total <- round((median(fit_m1a$Estimate[fit_m1a$population=="Bwindi"&
                                         fit_m1a$milestone_category=="social_affiliative"])+
                  median(fit_m1a$Estimate[fit_m1a$population=="Loango"&
                                           fit_m1a$milestone_category=="social_affiliative"])),3)
totalu <- round((median(fit_m1a$Q5[fit_m1a$population=="Bwindi"&
                                    fit_m1a$milestone_category=="social_affiliative"])+
                   median(fit_m1a$Q5[fit_m1a$population=="Loango"&
                                      fit_m1a$milestone_category=="social_affiliative"])),3)
totall <- round((median(fit_m1a$Q95[fit_m1a$population=="Bwindi"&
                                     fit_m1a$milestone_category=="social_affiliative"])+
                   median(fit_m1a$Q95[fit_m1a$population=="Loango"&
                                       fit_m1a$milestone_category=="social_affiliative"])),3)

diff/(total/2)*100
# 38.37209

diffu/(totalu/2)*100
# 43.16547

diffl/(totall/2)*100
# 35.09615

# 9 Figure 2 ####

xx <- met.brewer(name="Veronese", n=2, type="discrete")
fig2 <- ggplot(plot2_d, aes(x = factor(milestone_category), 
                             y = Estimate, 
                             fill = factor(species))) +
  geom_violin(position=position_dodge(width=0.9),
              width=1.1,
              trim=FALSE) +
  geom_boxplot(data=plot2_d,
               aes(x=factor(milestone_category),
                   y=Estimate,
                   col=factor(species)),
               fill="white",
               position=position_dodge(width=0.9),
               width=0.075,
               outlier.shape = NA)+
  geom_boxplot(position = position_dodge(width = 0.9), 
               alpha=0, width=0.075)+
  scale_fill_manual(values = xx) +
  xlab("Milestone Category") +
  ylab("Milestone emergence as proportion of developmental phase") +
  #ggtitle("Distribution of MPG by Number of Cylinders") +
  theme_minimal() +
  theme(axis.title=element_text(face="bold",size=18),
        axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=18),
        legend.text=element_text(size=18),
        legend.title=element_blank())
windows()
fig2

# 10 Model 2a ####

# data checks

sort(unique(m2d$population))
sort(unique(m2d$sex))
length(unique(m2d$subject))
length(unique(m2d$milestone))
sort(unique(m2d$milestone))
table(m2d$sex,m2d$milestone_category)
m2d$group_size <- as.numeric(as.character(m2d$group_size))
range(m2d$group_size)

hist(m2d$age)
hist(m2d$age_adjust)
m2d$zsize = (m2d$group_size - mean(m2d$group_size))/sd(m2d$group_size)
range(m2d$zsize)

mod = lm(age ~ zsize + 
           milestone_category +
           sex, 
         data = m2d)
vif(mod)

xx <- as.data.frame(vif(mod))
xx$model <- "Model 2"
xx$variable <- rownames(xx)
xx <- xx[,c(4,5,1,2,3)]

vif_table <- rbind(vif_table,xx)

m2prior <- get_prior(age ~ zsize + 
                       sex*milestone_category +
                       (1+zsize|id_code)+
                       (1+zsize|milestone),
                     data = m2d, family = lognormal())

m2prior = c(prior(lognormal(0.5,0.5), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd),
            prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = lognormal(),
              prior=m2prior)

start_time <- Sys.time()
m2 <- brm(age ~ zsize + 
            sex*milestone_category +
            (1+zsize|id_code)+
            (1+zsize|milestone),
          data = m2d, 
          family = lognormal(),
          prior=m2prior,
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
plot(m2)
summary(m2)
ppp2 <- pp_check(m2)
windows()
ppp2
ppp2 <- ppp2+labs(title="Model 2a: Weakly regularising priors")

mod2_table <- as.data.frame(round(posterior_summary(m2, probs = c(0.05, 0.95)),3))
mod2_table <- mod2_table[c(1:11),]
mod2_table$variable <- rownames(mod2_table)
rownames(mod2_table) <- NULL
mod2_valid <- as.data.frame(summary(m2)$fixed)
mod2_valid$model <- "Model 2"
mod2_valid <- mod2_valid[,c(6:8)]

m2_coef <- mcmc_plot(m2, variable = "^b_", regex=TRUE)
m2_coef <- m2_coef+labs(title="Model 2a: Weakly regularising priors")

# prior sensitivity checks

m2prior_uniform = c(prior(lognormal(0.5,100), class=Intercept),
                    prior(normal(0,100), class=b),
                    prior(exponential(1), class=sd),
                    prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = lognormal(),
              prior=m1prior_uniform)

start_time <- Sys.time()
m2uni <- brm(
  age ~ zsize + 
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = lognormal(),
  prior=m2prior_uniform,
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
ppp2uni <- ppp2uni+labs(title="Model 2a: Uniform priors")

m2uni_coef <- mcmc_plot(m2uni, variable = "^b_", regex=TRUE)
m2uni_coef <- m2uni_coef+labs(title="Model 2a: Uniform priors")

##

m2prior_strong = c(prior(lognormal(0.5,0.25), class=Intercept),
                   prior(normal(0,0.5), class=b),
                   prior(exponential(1), class=sd),
                   prior(student_t(3, 0, 2.5), class=sigma))

make_stancode(age ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = lognormal(),
              prior=m1prior_strong)

start_time <- Sys.time()
m2str <- brm(
  age ~ zsize + 
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m2d, 
  family = lognormal(),
  prior=m2prior_strong,
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
ppp2str <- ppp2str+labs(title="Model 2a: Strong priors")


m2str_coef <- mcmc_plot(m2str, variable = "^b_", regex=TRUE)
m2str_coef <- m2str_coef+labs(title="Model 2a: Strong priors")

m2_fixed_comparisons <- m2_coef+
  m2uni_coef+
  m2str_coef
windows()
m2_fixed_comparisons

m2_pppcs <- ppp2+
  ppp2uni+
    ppp2str
windows()
m2_pppcs


# model passes all checks

# 11 Model 2b ####

m2aprior <- get_prior(age_adjust ~ zsize + 
                       sex*milestone_category +
                       (1+zsize|id_code)+
                       (1+zsize|milestone),
                     data = m2d, family = Beta())

m2aprior = c(prior(normal(0,1.1), class=Intercept),
            prior(normal(0,1), class=b),
            prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize|milestone),
              data = m2d, 
              family = Beta(),
              prior=m2aprior)

start_time <- Sys.time()
m2a <- brm(age_adjust ~ zsize + 
            sex*milestone_category +
            (1+zsize|id_code)+
            (1+zsize|milestone),
          data = m2d, 
          family = Beta(),
          prior=m2aprior,
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
plot(m2a)
summary(m2a)
ppp2a <- pp_check(m2a)
windows()
ppp2a
ppp2a <- ppp2a+labs(title="Model 2b: Weakly regularising priors")

mod2a_table <- as.data.frame(round(posterior_summary(m2a, probs = c(0.05, 0.95)),3))
mod2a_table <- mod2a_table[c(1:11),]
mod2a_table$variable <- rownames(mod2a_table)
rownames(mod2a_table) <- NULL
mod2a_valid <- as.data.frame(summary(m2a)$fixed)
mod2a_valid$model <- "Model 2a"
mod2a_valid <- mod2a_valid[,c(6:8)]

m2a_coef <- mcmc_plot(m2a, variable = "^b_", regex=TRUE)
m2a_coef <- m2a_coef+labs(title="Model 2b: Weakly regularising priors")

# prior sensitivity checks

m2aprior_uniform = c(prior(normal(0,100), class=Intercept),
                    prior(normal(0,100), class=b),
                    prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m2d, 
              family = Beta(),
              prior=m2aprior_uniform)

start_time <- Sys.time()
m2auni <- brm(
  age_adjust ~ zsize + 
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m1d, 
  family = Beta(),
  prior=m2aprior_uniform,
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
plot(m2auni)
summary(m2auni)
windows()
pp_check(m2auni) # fine
ppp2auni <- pp_check(m2auni)
ppp2auni <- ppp2auni+labs(title="Model 2b: Uniform priors")

m2auni_coef <- mcmc_plot(m2uni, variable = "^b_", regex=TRUE)
m2auni_coef <- m2uni_coef+labs(title="Model 2b: Uniform priors")

##

m2aprior_strong = c(prior(normal(0,0.5), class=Intercept),
                   prior(normal(0,0.5), class=b),
                   prior(exponential(1), class=sd))

make_stancode(age_adjust ~ zsize + 
                sex*milestone_category +
                (1+zsize|id_code)+
                (1+zsize+sex|milestone),
              data = m1d, 
              family = Beta(),
              prior=m2aprior_strong)

start_time <- Sys.time()
m2astr <- brm(
  age_adjust ~ zsize + 
    sex*milestone_category +
    (1+zsize|id_code)+
    (1+zsize+sex|milestone),
  data = m2d, 
  family = Beta(),
  prior=m2aprior_strong,
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
plot(m2astr)
summary(m2astr)
windows()
pp_check(m2astr) # fine
ppp2astr <- pp_check(m2astr)
ppp2astr <- ppp2astr+labs(title="Model 2b: Strong priors")


m2astr_coef <- mcmc_plot(m2astr, variable = "^b_", regex=TRUE)
m2astr_coef <- m2astr_coef+labs(title="Model 2b: Strong priors")

m2a_fixed_comparisons <- m2a_coef+
  m2auni_coef+
  m2astr_coef
windows()
m2a_fixed_comparisons

m2a_pppcs <- ppp2a+
  ppp2auni+
  ppp2astr
windows()
m2a_pppcs


# model passes all checks

# 12 Summarise Models 1a and 1b ####

mod2_table
mod2a_table

# proportion support for m1 coefficients

postm2 <- posterior_samples(m2)
names(postm2)

round(sum(postm2$`b_zsize` < 0) / 
        length(postm2$`b_zsize`),3)

round(sum(postm2$`b_sexmale` < 0) / 
        length(postm2$`b_sexmale`),3)

round(sum(postm2$`b_milestone_categoryproximity` > 0) / 
        length(postm2$`b_milestone_categoryproximity`),3)

round(sum(postm2$`b_milestone_categoryphysical` > 0) / 
        length(postm2$`b_milestone_categoryphysical`),3)

round(sum(postm2$`b_milestone_categorysocial_affiliative` > 0) / 
        length(postm2$`b_milestone_categorysocial_affiliative`),3)

round(sum(postm2$`b_milestone_categorysocial_aggressive` > 0) / 
        length(postm2$`b_milestone_categorysocial_aggressive`),3)

round(sum(postm2$`b_sexmale:milestone_categoryphysical` < 0) / 
        length(postm2$`b_sexmale:milestone_categoryphysical`),3)

round(sum(postm2$`b_sexmale:milestone_categoryproximity` > 0) / 
        length(postm2$`b_sexmale:milestone_categoryproximity`),3)

round(sum(postm2$`b_sexmale:milestone_categorysocial_affiliative` > 0) / 
        length(postm2$`b_sexmale:milestone_categorysocial_affiliative`),3)

round(sum(postm2$`b_sexmale:milestone_categorysocial_aggressive` > 0) / 
        length(postm2$`b_sexmale:milestone_categorysocial_aggressive`),3)


postm2a <- posterior_samples(m2a)
names(postm2a)

round(sum(postm2a$`b_zsize` < 0) / 
        length(postm2a$`b_zsize`),3)

round(sum(postm2a$`b_sexmale` < 0) / 
        length(postm2a$`b_sexmale`),3)

round(sum(postm2a$`b_milestone_categoryproximity` > 0) / 
        length(postm2a$`b_milestone_categoryproximity`),3)

round(sum(postm2a$`b_milestone_categoryphysical` > 0) / 
        length(postm2a$`b_milestone_categoryphysical`),3)

round(sum(postm2a$`b_milestone_categorysocial_affiliative` > 0) / 
        length(postm2a$`b_milestone_categorysocial_affiliative`),3)

round(sum(postm2a$`b_milestone_categorysocial_aggressive` > 0) / 
        length(postm2a$`b_milestone_categorysocial_aggressive`),3)

round(sum(postm2a$`b_sexmale:milestone_categoryphysical` < 0) / 
        length(postm2a$`b_sexmale:milestone_categoryphysical`),3)

round(sum(postm2a$`b_sexmale:milestone_categoryproximity` < 0) / 
        length(postm2a$`b_sexmale:milestone_categoryproximity`),3)

round(sum(postm2a$`b_sexmale:milestone_categorysocial_affiliative` > 0) / 
        length(postm2a$`b_sexmale:milestone_categorysocial_affiliative`),3)

round(sum(postm2a$`b_sexmale:milestone_categorysocial_aggressive` > 0) / 
        length(postm2a$`b_sexmale:milestone_categorysocial_aggressive`),3)

# 13 Chimp-gorilla comparison plot ####

str(chimp_gorilla)

xx <- met.brewer(name="Pillement", n=3, type="discrete")

sort(unique(chimp_gorilla$Species))
chimp_gorilla$Species2 <- factor(chimp_gorilla$Species, 
                                 levels=c('Western chimpanzees', 
                                          'Mountain gorilla', 
                                          'Western gorilla'))
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

plot3 <- plot3a+plot3b
windows()  
plot3

# 14 Check supplementary materials info ####

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

#m1a

prior_table <- as.data.frame(prior_summary(m1))
prior_table$model <- "Model 1a: Weakly regularizing priors"

xx <- as.data.frame(prior_summary(m1str))
xx$model <- "Model 1a: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1uni))
xx$model <- "Model 1a: Uniform priors"

prior_table <- rbind(prior_table,xx)

#m1b

xx <- as.data.frame(prior_summary(m1a))
xx$model <- "Model 1b: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1a_str))
xx$model <- "Model 1b: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m1a_uni))
xx$model <- "Model 1b: Uniform priors"

prior_table <- rbind(prior_table,xx)

#m2a

xx <- as.data.frame(prior_summary(m2))
xx$model <- "Model 2a: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m2str))
xx$model <- "Model 2a: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m2uni))
xx$model <- "Model 2a: Uniform priors"

prior_table <- rbind(prior_table,xx)

#m2b

xx <- as.data.frame(prior_summary(m2a))
xx$model <- "Model 2b: Weakly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m2astr))
xx$model <- "Model 2b: Strongly regularizing priors"

prior_table <- rbind(prior_table,xx)

xx <- as.data.frame(prior_summary(m2uni))
xx$model <- "Model 2b: Uniform priors"

prior_table <- rbind(prior_table,xx)

# coef comparisons

windows()
m1_fixed_comparisons

windows()
m1a_fixed_comparisons

windows()
m2_fixed_comparisons

windows()
m2a_fixed_comparisons

# 15 Save and backup ####

load("C:/Users/besptkac/OneDrive/Project_data/gorilla_milestones/gorilla_milestones_analysis.Rdata")

save.image("E:/OneDrive/Project_data/gorilla_milestones/gorilla_milestones_materials.Rdata")
































