# ----------------------- #
### ---- MOD3 --------- ###
### --- BORAL --------- ###
### --- Poor Knight --- ###
# ----------------------- #

# First you need to install the following packages or programs: 
# JAGS from https://sourceforge.net/projects/mcmc-jags/files/
# install.packages("rjags")
# devtools::install_github('andyhoegh/UncertainOrd')


# setup -------------------------------------------------------------------
pacman::p_load(boral,
               corrplot,
               mvabund,
               readxl,
               rjags,
               UncertainOrd)


# prepare data ------------------------------------------------------------
dt_fish <- read_excel(path = "../../../Supplementary Informations/Anderson19Copula/ece34948-sup-0001-tables1.xlsx",
                      skip = 1)

# extract abundance data.frame from list 
Y <- dt_fish[,-c(1,2)]
X <- factor(dt_fish$Time, levels = c("Sep.98", "Mar.99", "Sep.99"))

# fit unconstrained Poisson and nb models 
fit_unconstrained_po <- boral(y = Y, 
                              family = "poisson", 
                              lv.control = list(num.lv = 2), 
                              row.eff = "fixed", 
                              save.model = T)
fit_unconstrained_nb <- boral(y = Y, 
                              family = "negative.binomial", 
                              lv.control = list(num.lv = 2), 
                              row.eff = "fixed", 
                              save.model = T)

# save model objects to file 
saveRDS(fit_unconstrained_po, "003_processed_data/boral/poor knight fish/boral_unconstrained_poisson.RDS")
saveRDS(fit_unconstrained_nb, "003_processed_data/boral/poor knight fish/boral_unconstrained_negbinom.RDS")

summary(fit_unconstrained_po)
plot(fit_unconstrained_po)
plot(fit_unconstrained_nb)
colors = c("green", "red", "blue")
lvsplot(fit_unconstrained_nb, biplot =F , col = colors[X])

# Fails uninformatively 
fit_constrained_nb <- boral(y = Y, X = X, family = "negative.binomial", lv.control = list(num.lv = 2), save.model = TRUE)
fit_constrained_nb <- boral(y = Y, X = X, family = "poisson", lv.control = list(num.lv = 2), save.model = TRUE)

# get MCMC samples from model 
samples <- get.mcmcsamples(fit_unconstrained_nb)
colnames(samples)
s1c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,1")
s2c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,2")
samples1 <- samples[,s1c]
samples2 <- samples[,s2c]
po0 = CredibleViz(coord1=samples1, coord2=samples2)
po1 <-CredibleViz(coord1=samples1, coord2=samples2,type = "scatter",items=c(2,13, 18))
po1