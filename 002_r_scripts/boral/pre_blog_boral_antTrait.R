# ----------------------- #
### ---- MOD3 --------- ###
### --- BORAL --------- ###
### --- AntTraits --- ###
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
data(antTraits)

# extract abundance data.frame from list 
Y <- antTraits$abund
X <- scale(antTraits$env)
traits = antTraits$traits
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
saveRDS(fit_unconstrained_po, "003_processed_data/boral/antTraits/boral_unconstrained_poisson.RDS")
saveRDS(fit_unconstrained_nb, "003_processed_data/boral/antTraits/boral_unconstrained_negbinom.RDS")

summary(fit_unconstrained_po)
plot(fit_unconstrained_po)
plot(fit_unconstrained_nb)
suppressMessages(lvsplot(fit_unconstrained_nb))

## uncertainty plots 
samples <- get.mcmcsamples(fit_unconstrained_nb)
colnames(samples)
s1c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,1")
s2c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,2")
samples1 <- samples[,s1c]
samples2 <- samples[,s2c]

options(warn=-1)

CredibleViz(coord1=samples1, coord2=samples2,type = "scatter",items=c(4, 21))



fit_constrained_nb <- boral(y = Y, X = X, family = "negative.binomial", lv.control = list(num.lv = 2), save.model = TRUE)
saveRDS(fit_constrained_nb, "003_processed_data/boral/antTraits/boral_constrained_negbinom.RDS")

par(mfrow=c(2,2))
plot(fit_constrained_nb)

envcors <- get.enviro.cor(fit_constrained_nb)
rescors <- get.residual.cor(fit_constrained_nb) 
par(mfrow = c(1,1))
corrplot(
        envcors$sig.cor,
        type = "lower",
        diag = FALSE,
        title =  "Correlations due to covariates", 
        mar = c(3,0.5,2,1), tl.srt = 45) 
corrplot(
        rescors$sig.cor,
        type = "lower",
        diag = FALSE,
        title =  "Residual correlations",
        mar = c(3, 0.5, 2, 1),
        tl.srt = 45
)
par(mfrow=c(1,1))
lvsplot(fit_constrained_nb)
par(mfrow=c(1,2))
lvsplot(fit_unconstrained_nb, main="Unconstrained ordination")
lvsplot(fit_constrained_nb, main="Constrained ordination")

summary(fit_unconstrained_nb)
# get MCMC samples from model 
samples <- get.mcmcsamples(fit_constrained_nb)
colnames(samples)
s1c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,1")
s2c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,2")
samples1 <- samples[,s1c]
samples2 <- samples[,s2c]

options(warn=-1)

CredibleViz(coord1=samples1, coord2=samples2,type = "scatter",items=c(4, 21))


# model with traits 
# full trait design 
# Non-numeric variables have to be removed
traits2 = traits[,-c(3,4)]

example_which_traits <- vector("list",ncol(X)+1)
for(i in 1:length(example_which_traits)) 
        example_which_traits[[i]] <- 1:ncol(traits2)




fit_constrained_nb_trait <- boral(y = Y, X = X, traits = traits2, which.traits = example_which_traits,family = "negative.binomial", lv.control = list(num.lv = 2), save.model = TRUE)
saveRDS(fit_constrained_nb_trait, "003_processed_data/boral/antTraits/boral_constrained_negbinom_trait.RDS")

corrplot(fit_constrained_nb_trait$geweke.diag$geweke.diag$traits.coefs, is.corr = FALSE)
