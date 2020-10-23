# ----------------- #
### ---- MOD3 --- ###
### --- BORAL --- ###
# ----------------- #

# First you need to install the following packages or programs: 
# JAGS from https://sourceforge.net/projects/mcmc-jags/files/
# install.packages("rjags")
# devtools::install_github('andyhoegh/UncertainOrd')


# setup -------------------------------------------------------------------
pacman::p_load(boral,
               corrplot,
               mvabund,
               rjags,
               UncertainOrd)


# prepare data ------------------------------------------------------------
#load spider data from mvabund package as example data set. 
data(spider)
# extract abundance data.frame from list 
Y <- spider$abund

# fit unconstrained poisson and nb models 
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
saveRDS(fit_unconstrained_po, "003_processed_data/boral_unconstrained_poisson.RDS")
saveRDS(fit_unconstrained_nb, "003_processed_data/boral_unconstrained_negbinom.RDS")

summary(fit_unconstrained_po)
plot(fit_unconstrained_po)
suppressMessages(lvsplot(fit_unconstrained_nb))

X <- scale(spider$x)
fit_constrained_nb <- boral(y = Y, X = X, family = "negative.binomial", lv.control = list(num.lv = 2), save.model = TRUE)
saveRDS(fit_constrained_nb, "003_processed_data/boral_constrained_negbinom.RDS")

par(mfrow=c(2,2))
plot(fit_constrained_nb)

envcors <- get.enviro.cor(fit_constrained_nb)
rescors <- get.residual.cor(fit_constrained_nb) 


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
lvsplot(fit.Xnb)
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

po1 <-CredibleViz(coord1=samples1, coord2=samples2,type = "point",items=c(1,20,26))
po2 <-CredibleViz(coord1=samples1, coord2=samples2,type = "scatter",items=c(1,20,26))
po3 <-CredibleViz(coord1=samples1, coord2=samples2,type = "circles",items=c(1,20,26))
po4 <-CredibleViz(coord1=samples1, coord2=samples2,type = "density",items=c(1,20,26))

uncertain_plot <- list(po1, po2, po3, po4)

saveRDS(uncertain_plot, "003_processed_data/boral_uncertainplot.RDS")


