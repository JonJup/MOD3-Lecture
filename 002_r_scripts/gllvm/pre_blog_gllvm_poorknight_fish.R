
save = TRUE

pacman::p_load(gllvm, corrplot,gclus, readxl)

dt_fish <- read_excel(path = "../../../Supplementary Informations/Anderson19Copula/ece34948-sup-0001-tables1.xlsx",
                      skip = 1)

# extract abundance data.frame from list 
y <- as.matrix(dt_fish[,-c(1,2)])
X <- factor(dt_fish$Time, levels = c("Sep.98", "Mar.99", "Sep.99"))
X <- as.matrix(X)
names(X) <- "season"


fit_uo_po <- gllvm(y, family = poisson())
fit_uo_nb <- gllvm(y, family = "negative.binomial")

if (save){
        
   saveRDS(fit_uo_po, "003_processed_data/gllvm/poor knigh fish//fit_model_uo_po.RDS")            
   saveRDS(fit_uo_nb, "003_processed_data/gllvm/poor knigh fish/fit_model_uo_nb.RDS")            
        
}

par(mfrow = c(1,2))
plot(fit_uo_nb, which = 1)
plot(fit_uo_po, which = 1) # fan shape in poisson 
plot(fit_uo_nb, which = 2)
plot(fit_uo_po, which = 2) # strong deviation in lower and higher qunaties 
plot(fit_uo_nb, which = 3)
plot(fit_uo_po, which = 3)

AIC(fit_uo_nb,fit_uo_po )

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
        
        saveRDS(fit_co_nb2, "003_processed_data/gllvm/poor knigh fish/fit_model_co_nb_2.RDS")            
        saveRDS(fit_co_nb3, "003_processed_data/gllvm/poor knigh fish/fit_model_co_nb_3.RDS")            
        
}

par(mfrow = c(1,2))
plot(fit_co_nb2, which = 1)
plot(fit_co_nb3, which = 1) # fan shape in poisson 
plot(fit_co_nb2, which = 2)
plot(fit_co_nb3, which = 2) # strong deviation in lower and higher qunaties 
plot(fit_co_nb2, which = 3)
plot(fit_co_nb3, which = 3)

AIC(fit_co_nb2, fit_co_nb3)

AIC(fit_co_nb2, fit_uo_nb)

ordiplot(fit_co_nb2, biplot = T)

coefplot(fit_co_nb2, cex.ylab = 0.7, mar = c(4,9,2,1),
         xlim.list = list(NULL, NULL, c(-4,4)))

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

