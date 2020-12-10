# ---------------- #
### --- HMSC --- ### 
### --- SDM ---- ###
### --- ME  ---- ###
# ---------------- #

# date written: 21.10.20

# Fit a single species distribution model to Corvus monedula. This is similar to
# the example from Ovaskainen & Abrego (2020) Chapter 5.7


# setup -------------------------------------------------------------------
pacman::p_load(abind, data.table, ggplot2, here, Hmsc, magrittr)
setwd(here())

data.directory="001_raw_data/hmsc_birds/data"
model.path="003_processed_data/hmsc_mcmc/corvusmonedula/"

# load data ---------------------------------------------------------------
data=fread(file.path(data.directory, "data.csv"))

# carpet data -------------------------------------------------------------
data=data[Year==2014]
data%<>%droplevels()
XData=as.data.frame(data[,c("Habitat", "AprMay")])
names(XData)=c("hab","clim")
Y=as.matrix(data$Corvus_monedula)
colnames(Y)="Corvus monedula"
xy=as.matrix(data[,c("x","y")])
studyDesign=data.frame(route=factor(data$Route))
rownames(xy)=studyDesign[,1]
rL=HmscRandomLevel(sData = xy)
XFormula=~hab+poly(clim, degree=2, raw=TRUE)

# fit models --------------------------------------------------------------
m_full = Hmsc(Y=Y,
              XData=XData,
              XFormula=XFormula,
              distr = "lognormal poisson",
              studyDesign = studyDesign,
              ranLevels=list(route=rL))

m_env = Hmsc(Y=Y,
              XData=XData,
              XFormula=XFormula,
              distr = "lognormal poisson")

m_space = Hmsc(Y=Y,
               XData=XData,
               XFormula=~1,
               distr = "lognormal poisson",
               studyDesign = studyDesign,
               ranLevels=list(route=rL))

models=list(m_full, m_env, m_space)
# MCMC parameters 
nChains = 2
thin=100
nParallel = max(round(parallel::detectCores()/2), nChains)
samples = 1000
transient=500*thin
verbose=500*thin

for (i in 1:1) {
        models[[i]] = sampleMcmc(
                models[[i]],
                thin = thin,
                samples = samples,
                transient = transient,
                nChains = nChains,
                verbose = verbose,
                initPar = "fixed effects"
        )
}

#saveRDS(models, "003_processed_data/hmsc_mcmc/corvusmonedula/temp_models.RDS")
# saveRDS(models, file.path(model.path,"temp_models_100.RDS"))
# saveRDS(models, file.path(model.path, "temp_models_10.RDS"))
models_1=readRDS(file.path(model.path,"temp_models.RDS"))
models_10=readRDS(file.path(model.path,"temp_models_10.RDS"))
models_100=readRDS(file.path(model.path, "hmsc_cm_models100.RDS"))

# In spNamesNumbers we decide weather to to name species by Names(T) or
# Numbers(F). Same for covariables in covNamesNumbers.
mpost = convertToCodaObject(models_100[[1]],
                            spNamesNumbers = c(T, F),
                            covNamesNumbers = c(T, F))

plot(mpost$Beta)

# effective sample size - beta 
ess.beta = effectiveSize(mpost$Beta)
#potential scale reduction factor - beta 
psrf.beta = gelman.diag(mpost$Beta,
                        multivariate = FALSE)$psrf
# Alpha is the spatial scale factor 
# effective sample size - alpha 
ess.alpha = effectiveSize(mpost$Alpha[[1]])
# potential scale reduction factor - Alpha 
psrf.alpha = gelman.diag(mpost$Alpha[[1]],
                         multivariate = FALSE)$psrf
# histograms 


MF = list()
for (i in 1:3){
        preds = computePredictedValues(models[[i]], expected = FALSE)
        MF[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
}
round(head(models[[1]]$X),2)
groupnames = c("habitat", "climate")
group = c(1,1,1,1,1,2,2)
VP = list()
for (i in 1:2){
        VP[[i]] = computeVariancePartitioning(models[[i]],
                                              group = group, groupnames = groupnames)
}
round(summary(mpost$Beta, quantiles = c(0.025, 0.5, 0.975))
      [[2]],2)

# prediction  -------------------------------------------------------------

m = models_100[[1]]
par(mfrow = c(1,2))
Gradient = constructGradient(m, focalVariable = "clim",
                             non.focalVariables = list(hab = 1))
predY = predict(m, Gradient = Gradient, expected = TRUE)
plotGradient(m, Gradient, pred = predY, measure = "Y",
             index = 1,showData = TRUE)
Gradient = constructGradient(m, focalVariable = "clim",
                             non.focalVariables = list(hab = 2))
mpost = convertToCodaObject(models[[3]])
round(summary(mpost$Alpha[[1]], quantiles = c(0.025, 0.5,
                                              0.975))[[2]], 2)
partition = createPartition(models[[1]], nfolds = 2,
                            column = "route")
MF = list()
for (i in 1:3){
        preds = computePredictedValues(models[[i]],
                                       partition = partition)
        MF[[i]] = evaluateModelFit(hM = models[[i]], predY = preds)
}

# predict on grid 

m = models_100[[1]]
grid = read.csv(file.path(data.directory,
                          "grid_10000.csv"))
grid = droplevels(subset(grid, !(Habitat=="Ma")))
xy.grid = as.matrix(cbind(grid$x, grid$y))
XData.grid = data.frame(hab = grid$Habitat,
                        clim = grid$AprMay)
Gradient = prepareGradient(m, XDataNew = XData.grid,
                           sDataNew = list(route = xy.grid))
predY = predict(m, Gradient = Gradient)
EpredY = apply(abind(predY,along = 3), c(1,2), mean)
EpredO = apply(abind(predY,along = 3), c(1,2), FUN =
                       function(a) {mean(a > 0)})
rowSums(predY)
mapData=data.frame(xy.grid, EpredY,EpredO)
names(mapData)=c("xCoordinates", "yCoordinates", "PredictedAbundance", "PredictedOccurence")
spO <- ggplot(data = mapData, 
             aes(x= xCoordinates, 
                 y= yCoordinates, 
                 color=PredictedOccurence)
             ) +
  geom_point(size=2)
spC <- ggplot(data = mapData, 
              aes(x= xCoordinates, 
                  y= yCoordinates, 
                  color=PredictedAbundance)
) +
  geom_point(size=2)

spO + 
  ggtitle("Predicted Corvus monedula occurrence") +
  xlab("East coordinate (km)") + 
  ylab("North coordinate (km)") + 
  scale_color_gradient(low = "blue", 
                       high="red", 
                       name ="Occurrence probability")
spC + 
  ggtitle("Predicted Corvus monedula abundance") +
  xlab("East coordinate (km)") + 
  ylab("North coordinate (km)") + 
  scale_color_gradient(low = "blue", 
                       high="red", 
                       name ="Abundance")
                                                                                                                                             "blue", high = "red", name = "Occurrence probability")
