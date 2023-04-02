# load packages
library(rstan)
# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(13)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
options(tidyverse.quiet = TRUE)
library(tidyverse)
options(dplyr.summarise.inform = FALSE)
library(bayesplot)

viabVr <- read.csv("https://www.dropbox.com/s/jfkgoxgv8o1fgqx/viability.csv?dl=1") #seed viability and germination


# seed viability
viabVr %>% 
  dplyr::select( plot, totS, yesMaybe, sr_f ) %>% 
  rename( SR        = sr_f,
          y_viab = yesMaybe,
          tot_seeds_viab = totS) %>% 
  dplyr::select(y_viab, tot_seeds_viab, SR ) %>% 
  na.omit->viab

# seed germination
viabVr %>% 
  dplyr::select( plot, germTot, germFail, sr_f ) %>% 
  rename( SR        = sr_f,
          y_germ    = germTot ) %>% 
  mutate(tot_seeds_germ = y_germ + germFail ) %>% 
  dplyr::select(y_germ, tot_seeds_germ, SR ) %>% 
  na.omit->germ

data_viab_germ <- list( # viability
                       n_v       = nrow(viab),
                       y_v       = viab$y_viab,
                       tot_seeds_v = viab$tot_seeds_viab,
                       SR_v        = viab$SR,
                       
                       # germination
                       n_m       = nrow(germ),
                       y_m       = germ$y_germ,
                       tot_seeds_m = germ$tot_seeds_germ,
                       SR_m        = germ$SR)

sim_pars <- list(
  warmup = 1000,
  iter = 4000,
  thin = 3,
  chains = 3
)

fit_viab_germ <- stan(
  file = "stan/poar_viab_germ.stan",
  data = data_viab_germ,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains)
mcmc_trace(fit_viab_germ,pars=c("v0","a_v","phi_v","m","phi_m"))

stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'
mod <- stan_model(model_code = stancode, verbose = TRUE)
