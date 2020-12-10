# ----------------- #
### --- HMSC ---- ### 
### --- JSDM ---- ###
### --- ME  ----- ###
# ----------------- #


# Taken almost one to one from Tikhonov et al 2019 SI2 

# setup ------------------------------------------------------------------
if (!require(pacman)) install.packages(pacman)
p_load(ape, colorspace, ggplot2, here, Hmsc, data.table, dplyr, vioplot)

setwd(here())

## dirs -- 
dir_data = "001_raw_data/hmsc_birds/data/"
dir_model = "003_processed_data/hmsc_mcmc/jsdm"

# SPECIES AND ENVIRONMENTAL DATA 
dt_y = fread(file.path(dir_data, "data.csv"))
# ENVIRONMENTAL COVARIATES 
dt_x = dt_y[,c(5,6,7,8,9)]
# PHYLOGENY 
ph_phylo <- read.tree(file.path(dir_data, "CTree.tre"))
#TRAITS 
dt_traits = fread(file.path(dir_data, "traits.csv")) 
dt_traits$LogMass = log(dt_traits$Mass)


# subset species to 10 most abundant species ----------------------------------------------------
vc_subset = sort(colSums(dt_y[,10:59]), decreasing = TRUE)
vc_subset = names(vc_subset)[1:9]
vc_subset = append(vc_subset, "Corvus_monedula")
vc_subset = which(names(dt_y) %in% vc_subset)

dt_y=dt_y[,append(1:9, vc_subset), with = F]
ma_y=as.matrix(dt_y[,10:19])
# trimm traits
dt_traits = dt_traits[Species %in% colnames(dt_y)]


# STUDY DESIGN 
ma_studydesign = matrix(NA,nrow(ma_y),2) 
ma_studydesign[,1] = sprintf('Route_%.3d',dt_y$Route) 
ma_studydesign[,2] = sprintf('Year_%.3d',dt_y$Year)
df_studydesign = as.data.frame(ma_studydesign) 
colnames(df_studydesign) = c("Route","Year") 
df_studydesign[,1]=as.factor(df_studydesign[,1]) 
df_studydesign[,2]=as.factor(df_studydesign[,2])

# RANDOM EFFECT STRUCTURE, HERE ROUTE AS A SPATIAL LATENT VARIABLE 
vc_routes = levels(df_studydesign[,1]) 
in_nroutes = length(vc_routes) 
ma_xy = matrix(0, nrow = in_nroutes, ncol = 2) 
for (i in 1:in_nroutes){ 
        rows=df_studydesign[,1]==vc_routes[[i]] 
        ma_xy[i,1] = mean(dt_y[rows,]$x) 
        ma_xy[i,2] = mean(dt_y[rows,]$y)
}
colnames(ma_xy) = c("x","y") 
ma_sRL = ma_xy 
rownames(ma_sRL) = vc_routes 
rL = HmscRandomLevel(sData=ma_sRL) 
rL$nfMin = 5 
rL$nfMax = 10

XFormula = ~ Habitat + poly(AprMay, degree = 2, raw = TRUE)
TrFormula = ~Migration + LogMass

m = Hmsc(Y = ma_y, 
         XData = as.data.frame(dt_x), 
         XFormula = XFormula, 
         TrData = dt_traits, 
         TrFormula = TrFormula, 
         phyloTree = ph_phylo, 
         distr = "lognormal poisson", 
         studyDesign = df_studydesign, 
         ranLevels = list(Route=rL))



thin = c(5)
samples = 1000
nChains = 4
set.seed(1)
computational.time=c() 
m_fit = list()
# for(i in seq_along(thin)){
#         ptm = proc.time()
#         
#         m_fit[[i]] = sampleMcmc(
#                 m,
#                 samples = samples,
#                 thin = thin[i],
#                 adaptNf = rep(ceiling(0.4 * samples * thin), 1),
#                 transient = ceiling(0.5 * samples * thin),
#                 nChains = nChains,
#                 nParallel = nChains,
#                 initPar = "fixed effects"
#         )
#         computational.time[i] = proc.time() - ptm
#         filename = file.path(dir_model,
#                              paste0("model_",
#                                     thin[i],
#                                     ".RDS"))
#         saveRDS(m_fit[[i]], file=filename)
#         rm(ptm);gc()
#         print(paste(i, Sys.time()))
# }
m_fit <- readRDS("003_processed_data/hmsc_mcmc/jsdm/model_5.RDS")
mpost = convertToCodaObject(m_fit,
                            spNamesNumbers = c(T, F),
                            covNamesNumbers = c(T, F))

es.beta = effectiveSize(mpost$Beta) 
ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
es.gamma = effectiveSize(mpost$Gamma) 
ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf

es.rho = effectiveSize(mpost$Rho) 
ge.rho = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf 
es.V = effectiveSize(mpost$V) 
ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf

es.omega = effectiveSize(mpost$Omega[[1]])
ge.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
mixing = list(
        es.beta = es.beta,
        ge.beta = ge.beta,
        es.gamma = es.gamma,
        ge.gamma = ge.gamma,
        es.rho = es.rho,
        ge.rho = ge.rho,
        es.V = es.V,
        ge.V = ge.V,
        es.omega = es.omega,
        ge.omega = ge.omega
)
par(mfrow=c(1,2))
hist(mixing$es.beta)
hist(mixing$ge.beta)
hist(mixing$es.gamma)
hist(mixing$ge.gamma)
hist(mixing$es.rho)
hist(mixing$ge.rho)
hist(mixing$es.V)
hist(mixing$ge.V)
hist(mixing$es.omega)
hist(mixing$ge.omega)

## modelfit 
# predY = computePredictedValues(hM=m_fit, expected=FALSE)
# MF = evaluateModelFit(hM=m_fit, predY=predY)
# saveRDS(predY, file=file.path(dir_model, "prediction.RDS"))
# saveRDS(MF, file=file.path(dir_model, "model_fit.RDS"))
predY = readRDS(file=file.path(dir_model, "prediction.RDS"))
MF    = readRDS(file=file.path(dir_model, "model_fit.RDS"))

# -- test 
plot_data = data.frame (species = colnames(predY), 
                        RMSE = MF$RMSE,
                        SR2 = MF$SR2,
                        O.AUC = MF$O.AUC,
                        O.TjurR2 = MF$O.TjurR2,
                        O.RMSE = MF$O.RMSE,
                        C.SR2 = MF$C.SR2,
                        C.RMSE = MF$C.RMSE
                        )

## -- variation partitioning 
VP = computeVariancePartitioning(m_fit)
vals = VP$vals
mycols = rainbow(nrow(VP$vals))
par(mfrow=c(1,1))
plotVariancePartitioning(hM=m_fit, VP=VP,cols = mycols, args.leg=list(bg="white",cex=0.7),
                         main = paste0("Proportion of explained variance, ",cex.main=0.8))

postBeta = getPostEstimate(m_fit, parName="Beta")
show.sp.names = (is.null(m_fit$phyloTree)) 
png(width=1024, height=768, filename = "004_plots/R_hmsc_beta.png")
plotBeta(m, post=postBeta, supportLevel = 0.95,param="Sign",
         plotTree = !is.null(m_fit$phyloTree),
         covNamesNumbers = c(TRUE,FALSE),
         spNamesNumbers=c(TRUE,FALSE),
         cex=c(1,1,1),
         split = 0.6)
dev.off()
mymain = paste0("BetaPlot, ",modelnames[[j]])

## -- Gamma plot 
postGamma = getPostEstimate(m_fit, parName="Gamma")
png(width=1024, height=768, filename = "004_plots/R_hmsc_gamma.png")
plotGamma(m_fit, post=postGamma, supportLevel = 0.9, param="Sign",
          covNamesNumbers = c(TRUE,FALSE),
          trNamesNumbers=c(TRUE,FALSE),
          cex=c(1,1,1))
dev.off()
## -- Omega plot 

OmegaCor = computeAssociations(m_fit)
supportLevel = 0.95

plotOrder = corrMatOrder(OmegaCor[[1]]$mean, order = "AOE")
toPlot = ((OmegaCor[[1]]$support > supportLevel) + (OmegaCor[[1]]$support <
                                                            (1 - supportLevel)) > 0) * sign(OmegaCor[[1]]$mean)
mymain = paste0("Associations: ",names(m_fit$ranLevels)[[1]])

corrplot(
        toPlot[plotOrder, plotOrder],
        method = "color",
        col = colorRampPalette(c("blue", "white", "red"))(3),
        mar = c(0, 0, 1, 0),
        main = mymain,
        cex.main = 0.8
)



## -- s8 predictions -- ## 

nm = 1
j= 1
m = m_fit
covariates = all.vars(m$XFormula)
ex.sp = which.max(colMeans(m$Y, na.rm = TRUE)) #most common species as example species
m$XData$Habitat = factor(m$XData$Habitat)
Gradient = constructGradient(m, focalVariable="AprMay", non.focalVariables=list(Habitat=list("3", "Urb")))
predY = predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew, ranLevels=Gradient$rLNew, expected=TRUE)
saveRDS(Gradient, )
saveRDS(predY, "003_processed_data/hmsc_mcmc/jsdm/prediction2.RDS")
plotGradient(m, Gradient, pred=predY, measure="S", las=1, showData = TRUE, main='Species richness (measure="S")')
plotGradient(m, Gradient, pred=predY, measure="Y", index=1, las=1, showData = TRUE, main='Focal species occurrence (measure="Y")')
plotGradient(m, Gradient, pred=predY, measure="T", index=3, las=1, showData = TRUE, main='Mean trait value (measure="T")')



        #predY = predict(m, Gradient = Gradient, expected = TRUE)
        predY2 = predict(m, Gradient = Gradient2, expected = TRUE)
        pl = plotGradient(
                m,
                Gradient2,
                pred = predY2,
                yshow = 0,
                measure = "S",
                showData = TRUE,
                main = paste0("summed response (marginal effect)")
        )
        if (inherits(pl, "ggplot")) {
                print(pl + labs(
                        title = paste0(
                                modelnames[[j]],
                                ": summed response (marginal effect)"
                        )
                ))
        }
        par(mfrow = c(2, 1))
        pl = plotGradient(
                m,
                Gradient,
                pred = predY,
                yshow = if (m$distr[1, 1] == 2) {
                        c(-0.1, 1.1)
                } else{
                        0
                },
                measure = "Y",
                index = ex.sp,
                showData = TRUE,
                main = paste0(modelnames[[j]], ": example species (total effect)")
        )
        if (inherits(pl, "ggplot")) {
                print(pl + labs(
                        title = paste0(modelnames[[j]], ": example species (total effect)")
                ))
        }
        pl = plotGradient(
                m,
                Gradient2,
                pred = predY2,
                yshow = if (m$distr[1, 1] == 2) {
                        c(-0.1, 1.1)
                } else{
                        0
                },
                measure = "Y",
                index = ex.sp,
                showData = TRUE,
                main = paste0(modelnames[[j]], ": example species (marginal effect)")
        )
        if (inherits(pl, "ggplot")) {
                print(pl + labs(
                        title = paste0(
                                modelnames[[j]],
                                ": example species (marginal effect)"
                        )
                ))
        }
        if (m$nt > 1) {
                for (l in 2:m$nt) {
                        par(mfrow = c(2, 1))
                        pl = plotGradient(
                                m,
                                Gradient,
                                pred = predY,
                                measure = "T",
                                index = l,
                                showData = TRUE,
                                yshow = 0,
                                main = paste0(
                                        modelnames[[j]],
                                        ": community weighted mean trait (total effect)"
                                )
                        )
                        if (inherits(pl, "ggplot")) {
                                print(pl + labs(
                                        title = paste0(
                                                modelnames[[j]],
                                                ": community weighted mean trait (total effect)"
                                        )
                                ))
                        }
                        pl = plotGradient(
                                m,
                                Gradient2,
                                pred = predY2,
                                measure = "T",
                                index = l,
                                showData = TRUE,
                                yshow = 0,
                                main = paste0(
                                        modelnames[[j]],
                                        ": community weighted mean trait (marginal effect)"
                                )
                        )
                        if (inherits(pl, "ggplot")) {
                                print(pl + labs(
                                        title = paste0(
                                                modelnames[[j]],
                                                ": community weighted mean trait (marginal effect)"
                                        )
                                ))
                        }
                }
        }
}
