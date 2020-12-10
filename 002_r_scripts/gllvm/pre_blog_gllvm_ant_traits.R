
save = FALSE

pacman::p_load(gllvm, corrplot,gclus)
data("antTraits")
y <- as.matrix(antTraits$abund)
X <- scale(as.matrix(antTraits$env))
TR <- antTraits$traits
fit_uo_po <- gllvm(y, family = poisson())
fit_uo_nb <- gllvm(y, family = "negative.binomial")

if (save){
        
   saveRDS(fit_uo_po, "003_processed_data/gllvm/anttraits/fit_model_uo_po.RDS")            
   saveRDS(fit_uo_nb, "003_processed_data/gllvm/fit_model_uo_nb.RDS")            
        
}

par(mfrow = c(1,2))
plot(fit_uo_nb, which = 1)
plot(fit_uo_po, which = 1)
plot(fit_uo_nb, which = 2)
plot(fit_uo_po, which = 2)
plot(fit_uo_nb, which = 3)
plot(fit_uo_po, which = 3)

ordiplot(fit_uo_nb,
         biplot = TRUE,
         ind.spp = 15,
         xlim = c(-3,3) ,
         ylim = c(-2, 1.6))


fit_co_nb2 <- gllvm(y, X, num.lv = 2,
                    formula = ~ Bare.ground + Shrub.cover +
                            Volume.lying.CWD,
                    family = "negative.binomial")
fit_co_nb3 <- gllvm(y, X, num.lv = 3,
                    formula = ~ Bare.ground + Shrub.cover +
                            Volume.lying.CWD,
                    family = "negative.binomial")

if (save){
        
        saveRDS(fit_co_nb2, "003_processed_data/gllvm/anttraits/fit_model_co_nb_2.RDS")            
        saveRDS(fit_co_nb3, "003_processed_data/gllvm/anttraits/fit_model_co_nb_3.RDS")            
        
}

AIC(fit_co_nb2, fit_co_nb3)


fit_4th <- gllvm(
        y = y,
        X = X,
        TR = TR,
        family = "negative.binomial",
        num.lv = 2,
        formula = y ~
                (Bare.ground + Shrub.cover + Volume.lying.CWD) +
                (Bare.ground + Shrub.cover + Volume.lying.CWD):(Pilosity + Polymorphism + Webers.length))

if (save){
        
        saveRDS(fit_4th, "003_processed_data/gllvm/anttraits/fit_model_trait_co_nb_r2.RDS")            
            
        
}

fourth <- fit_4th$fourth.corner
colort <- colorRampPalette(c("blue", "white", "red"))
a <- max( abs(fourth) )
plot.4th <- lattice::levelplot((as.matrix(fourth)), xlab = "Environmental Variables",
                               ylab = "Species traits", col.regions = colort(100), cex.lab =1.3,
                               at = seq(-a, a, length = 100), scales = list(x = list(rot = 45)))
plot.4th


fit_4th2 <- gllvm(y, X, TR, family = "negative.binomial", 
                  num.lv = 2, formula = y ~ (Bare.ground + Shrub.cover + Volume.lying.CWD))
trait_anoa = anova(fit_4th, fit_4th2)

if (save){
        saveRDS(fit_4th2, "003_processed_data/gllvm/anttraits/fit_model_traits_not_co_nb_2.RDS")
        saveRDS(trait_anoa, "003_processed_data/gllvm/anttraits/anova_traits.RDS")
}