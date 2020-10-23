# ---------------------------------- #
### ------------- MOD3 ----------- ###
### --- Visuallizg Uncertainty --- ###
### ------- in BORAL Models ------ ###
# ---------------------------------- #

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
# extract and scale environmental data from list 
X <- scale(spider$x)
# fit constrained negative binomial boral model  
fit_constrained_nb <- boral(y = Y, 
                            X = X, 
                            family = "negative.binomial", 
                            lv.control = list(num.lv = 2), 
                            save.model = TRUE)

# get MCMC samples from model 
samples <- get.mcmcsamples(fit_constrained_nb)
# extract ID of columns with relevant data. First for the first latent variable then for the second. 
s1c <- grepl(x=colnames(samples), 
             pattern="lvs+.+\\d,1")
s2c <- grepl(x=colnames(samples), 
             pattern="lvs+.+\\d,2")
# extract those columns from the samples
samples1 <- samples[,s1c]
samples2 <- samples[,s2c]

# Plot uncertanty plots 
#The type argument determines the type of plot The items argument determines for
#which items the uncertainty should be displayed. Not applicable to type="point"
#as this only shows the ordination plot.

po1 <-CredibleViz(coord1=samples1, 
                  coord2=samples2,
                  type = "point")

po2 <-CredibleViz(coord1=samples1, 
                  coord2=samples2,
                  type = "scatter",
                  items=c(1,20,26))

po3 <-CredibleViz(coord1=samples1,
                  coord2=samples2,
                  type = "circles",
                  items=c(1,20,26))

po4 <-CredibleViz(coord1=samples1,
                  coord2=samples2,
                  type = "density",
                  items=c(1,20,26))