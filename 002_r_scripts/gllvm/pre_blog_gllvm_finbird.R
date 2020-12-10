
save = TRUE

pacman::p_load(data.table, gllvm, corrplot,gclus, magrittr)

data = fread("001_raw_data/hmsc_birds/data/data.csv")
TR = fread("001_raw_data/hmsc_birds/data/traits.csv")

y <- as.matrix(data[,10:59])
x = data[,2:6]
x$Habitat %<>% factor()
fit_uo_po <- gllvm(y, family = poisson())
fit_uo_nb <- gllvm(y, family = "negative.binomial")

if (save){
        
   saveRDS(fit_uo_po, "003_processed_data/gllvm/finbird/fit_model_uo_po.RDS")            
   saveRDS(fit_uo_nb, "003_processed_data/gllvm/finbird/fit_model_uo_nb.RDS")            
        
}

par(mfrow = c(1,2))
plot(fit_uo_nb, which = 1)
plot(fit_uo_po, which = 1) #fan shaped
plot(fit_uo_nb, which = 2)
plot(fit_uo_po, which = 2) # very bad 
plot(fit_uo_nb, which = 3)
plot(fit_uo_po, which = 3)

ordiplot(fit_uo_nb,
         biplot = TRUE,
         ind.spp = 15,
         xlim = c(-3,3) ,
         ylim = c(-2, 1.6))


fit_co_nb2 <- gllvm(y, x, num.lv = 2,
                    formula = ~.,
                    family = "negative.binomial")
fit_co_nb3 <- gllvm(y, x, num.lv = 3,
                    formula = ~ .,
                    family = "negative.binomial")

if (save){
        
        saveRDS(fit_co_nb2, "003_processed_data/gllvm/finbird/fit_model_co_nb_2.RDS")            
        saveRDS(fit_co_nb3, "003_processed_data/gllvm/finbird/fit_model_co_nb_3.RDS")            
        
}

AIC(fit_co_nb2, fit_co_nb3) # 2d better than 3d 
AIC(fit_co_nb2, fit_uo_nb) # constrained better 

par(mfrow = c(1,2))
plot(fit_co_nb2, which = 1)
plot(fit_co_nb3, which = 1)
plot(fit_co_nb2, which = 2)
plot(fit_co_nb3, which = 2) 
plot(fit_co_nb2, which = 3)
plot(fit_co_nb3, which = 3)

AIC(fit_co_nb2, fit_co_nb3)

AIC(fit_co_nb2, fit_uo_nb)

ordiplot(fit_co_nb2, biplot = T)

coefplot(fit_co_nb2, cex.ylab = 0.7, # mar = c(4,9,2,1),
         #xlim.list = list(NULL, NULL, c(-4,4)),
         which.Xcoef = 1:4)
coefplot(fit_co_nb2, cex.ylab = 0.7, # mar = c(4,9,2,1),
         #xlim.list = list(NULL, NULL, c(-4,4)),
         which.Xcoef = 5)

par(mfrow = c(1,1))
cr <- getResidualCor(fit_co_nb2)
corrplot(cr[order.single(cr),
            order.single(cr)],
         diag = FALSE,
         type = "lower",
         method = "square",
         tl.cex = 0.8,
         tl.srt = 45,
         tl.col = "red")




# Trait model  ------------------------------------------------------------
TR = as.data.frame(TR)
row.names(TR) = TR$Species
tr = dplyr::select(TR, - "Species")
#tr = dplyr::select(tr, - "Migration")
tr$Migration = factor(tr$Migration)
fit_4th <- gllvm(
        y = y,
        X = x,
        TR = tr,
        family = "negative.binomial",
        num.lv = 2,
        formula = y ~
                (Year + x + y + Effort + Habitat) +
                (Year + x + y + Effort + Habitat):(Mass + Migration + Urb + Br + Co + Op + Ma + We))




if (save){
        
        saveRDS(fit_4th, "003_processed_data/gllvm/finbird/fit_model_trait_co_nb_r2.RDS")            
            
        
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