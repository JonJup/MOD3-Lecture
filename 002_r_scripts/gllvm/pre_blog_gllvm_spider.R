
save = TRUE

pacman::p_load(gllvm, corrplot,gclus, readxl)

data(spider)

# extract abundance data.frame from list 
y <- as.matrix(spider$abund)
X <- as.matrix(scale(spider$x))

fit_uo_po <- gllvm(y, family = poisson())
fit_uo_nb <- gllvm(y, family = "negative.binomial")

if (save){
        
        saveRDS(fit_uo_po, "003_processed_data/gllvm/hspider/fit_model_uo_po.RDS")            
        saveRDS(fit_uo_nb, "003_processed_data/gllvm/hspider/fit_model_uo_nb.RDS")            
        
}

par(mfrow = c(1,2))
plot(fit_uo_nb, which = 1) # upward trend 
plot(fit_uo_po, which = 1) 
plot(fit_uo_nb, which = 2)
plot(fit_uo_po, which = 2) # NB fits better
plot(fit_uo_nb, which = 3)
plot(fit_uo_po, which = 3) # both equally 

AIC(fit_uo_nb,fit_uo_po ) # nb has lower AIC 

ordiplot(fit_uo_nb,
         biplot = TRUE,
         ind.spp = 15,
         xlim = c(-3,3) ,
         ylim = c(-2, 1.6))


fit_co_nb2 <- gllvm(y, X, num.lv = 2,
                    formula = ~ .,
                    family = "negative.binomial")
fit_co_nb3 <- gllvm(y, X, num.lv = 3,
                    formula = ~.,
                    family = "negative.binomial")

if (save){
        
        saveRDS(fit_co_nb2, "003_processed_data/gllvm/hspider/fit_model_co_nb_2.RDS")            
        saveRDS(fit_co_nb3, "003_processed_data/gllvm/hspider/fit_model_co_nb_3.RDS")            
        
}

fit_co_nb2 = readRDS("003_processed_data/gllvm/hspider/fit_model_co_nb_2.RDS")

par(mfrow = c(1,2))
plot(fit_co_nb2, which = 1)
plot(fit_co_nb3, which = 1) # slight upward tild 
plot(fit_co_nb2, which = 2)
plot(fit_co_nb3, which = 2) # strong deviation in lower and higher qunaties 
plot(fit_co_nb2, which = 3)
plot(fit_co_nb3, which = 3)

AIC(fit_co_nb2, fit_co_nb3) # 2 is better

AIC(fit_co_nb2, fit_uo_nb) # CO is better

ordiplot(fit_co_nb2, biplot = T)

par(mfrow=c(2,3))
coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)), which.Xcoef = 1)
coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)), which.Xcoef = 2)
coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)), which.Xcoef = 3)
coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)), which.Xcoef = 4)
coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)), which.Xcoef = 5)
coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)), which.Xcoef = 6)

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
