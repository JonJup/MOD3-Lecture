# ------------------ #
### ---- MOD3 --- ###
### --- BORAL --- ###
### --- BIRDS --- ###
# ----------------- #

# First you need to install the following packages or programs: 
# JAGS from https://sourceforge.net/projects/mcmc-jags/files/
# install.packages("rjags")
# devtools::install_github('andyhoegh/UncertainOrd')


# setup -------------------------------------------------------------------
pacman::p_load(boral,
               corrplot,
               data.table,
               mvabund,
               readxl,
               rjags,
               UncertainOrd)


# prepare data ------------------------------------------------------------
data = fread("001_raw_data/hmsc_birds/data/data.csv")
traits = fread("001_raw_data/hmsc_birds/data/traits.csv")

subset = sample(1:nrow(data), 300)

data = data[subset, ]

env = data[,1:9]
data = data[,-c(1:9)]
env[,c("Year", "Habitat") := .(factor(Year), factor(Habitat))]

# extract abundance data.frame from list 
Y = data
X =  env
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
saveRDS(fit_unconstrained_po, "003_processed_data/boral/finbird//boral_unconstrained_poisson.RDS")
saveRDS(fit_unconstrained_nb, "003_processed_data/boral/finbird/boral_unconstrained_negbinom.RDS")

fit_unconstrained_po = readRDS("003_processed_data/boral/finbird/boral_unconstrained_poisson.RDS")
fit_unconstrained_nb = readRDS("003_processed_data/boral/finbird/boral_unconstrained_negbinom.RDS")

summary(fit_unconstrained_po)
par(mfrow=c(2,2))
plot(fit_unconstrained_po)
plot(fit_unconstrained_nb)

par(mfrow=c(1,1))
lvsplot(fit_unconstrained_nb)
lvsplot(fit_unconstrained_nb, biplot=F)


# add arrows 
su_nb = summary(fit_unconstrained_nb)
cor(su_nb$lvs, X[,c(3,4,5,7,8,9)])

X
X_sub = X[,c(3,4,5,7,8,9)]
X_sub = scale(X_sub)
fit_constrained_nb <- boral(y = Y, X = X_sub, family = "negative.binomial", lv.control = list(num.lv = 2), save.model = TRUE)
saveRDS(fit_constrained_nb, "003_processed_data/boral/finbird/boral_constrained_negbinom.RDS")

par(mfrow=c(2,2))
plot(fit_constrained_nb)

envcors <- get.enviro.cor(fit_constrained_nb)
rescors <- get.residual.cor(fit_constrained_nb) 

saveRDS(envcors, "003_processed_data/boral/finbird/envcors.RDS")
saveRDS(rescors, "003_processed_data/boral/finbird/rescors.RDS")

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
lvsplot(fit_lvmp, main="Unconstrained ordination")
lvsplot(fit.Xnb, main="Constrained ordination")

# get MCMC samples from model 
samples <- get.mcmcsamples(fit_constrained_nb)
colnames(samples)
s1c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,1")
s2c <- grepl(x=colnames(samples), pattern="lvs+.+\\d,2")
samples1 <- samples[,s1c]
samples2 <- samples[,s2c]

options(warn=-1)

po1 <-CredibleViz(coord1=samples1, coord2=samples2,type = "scatter",items=c(40,80,245))


# model with traits 
# full trait design 
# Non-numeric variables have to be removed. This al
traits2 = data.frame(mass = traits$Mass)
rownames(traits2) = traits$Species
example_which_traits <- vector("list",ncol(X))
for(i in 1:length(example_which_traits)) 
        example_which_traits[[i]] <- 1:ncol(traits2)

# I got an error that suggested JAGS was unable to resolve the habitat parameter so i remove it 
fit_constrained_nb_trait <- boral(y = Y, X = X[,-6], traits = traits2, which.traits = example_which_traits,family = "negative.binomial", lv.control = list(num.lv = 2), save.model = TRUE)
saveRDS(fit_constrained_nb_trait, "003_processed_data/boral/finbird/boral_constrained_negbinom_trait.RDS")

corrplot(fit_constrained_nb_trait$geweke.diag$geweke.diag$traits.coefs, is.corr = FALSE)
