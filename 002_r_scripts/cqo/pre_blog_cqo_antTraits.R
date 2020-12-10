pacman::p_load(data.table, dplyr, ggplot2, purrr, VGAM, mvabund)
source("002_r_scripts/cqo/cqo_residual_plot.R")

data("antTraits")

bio = antTraits$abund
env = antTraits$env
env = scale(env)


# In Yee (2015, p.241), Thomas Yee suggests to keep S <= 10, n <=500 and p <=5. We will stick to this. 
# We have 41 species so I pick the 10 most abundand ones. 
colSums(bio) %>% 
        sort(decreasing=T) %>% 
        names %>% 
        .[11:10]

comb = cbind(bio,env)

# A df1.nl of 1.5 approximates CQO. By setting df1.nl to 1 we are conservative so that species will rather not look unimodal. 
set.seed = 1
cao_p = cao(formula = cbind(Iridomyrmex.rufoniger,Pheidole.sp..A,Rhytidoponera.metallica.sp..A,Pheidole.sp..E,
                            Iridomyrmex.bicknelli,Camponotus.consobrinus,Monomorium.leae,Iridomyrmex.mjobergi,Heteroponera.sp..A,Nylanderia.sp..A) ~
                    Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
            family = poissonff, data = comb, Rank = 1, df1.nl = 2, Bestof = 30 
)
saveRDS(cao_p, "003_processed_data/cao/antTraits/cao.RDS")

cao_p_summary = summary(cao_p)
sum(round(cao_p_summary@misc$deviance.Bestof,0) == round(min(cao_p_summary@misc$deviance.Bestof), 0))
x11()
par(mfrow=c(3,4))
plot(cao_p, lcol = "blue", lwd = 2, ylim = c(-5, 5), xlab = "", ylab = "")

# Most species are not unimodal so I kick out those an replace them with new ones 
cao_p = cao(formula = cbind(Iridomyrmex.rufoniger, Iridomyrmex.bicknelli,Camponotus.consobrinus,
                            Monomorium.rothsteini, Notoncus.ectatommoides, Tapinoma.sp..A, Rhytidoponera.sp..B,
                            Iridomyrmex.suchieri,Melophorus.sp..F,Aphaenogaster.longiceps) ~
                    Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
            family = poissonff, data = comb, Rank = 1, df1.nl = 2, Bestof = 30 
)
cao_p_summary = summary(cao_p)
sum(round(cao_p_summary@misc$deviance.Bestof,0) == round(min(cao_p_summary@misc$deviance.Bestof), 0))
x11()
par(mfrow=c(3,4))
plot(cao_p, lcol = "blue", lwd = 2, ylim = c(-5, 5), xlab = "", ylab = "")
saveRDS(cao_p, "003_processed_data/cao/antTraits/cao2.RDS")
# one more time 
cao_p = cao(formula = cbind(Iridomyrmex.rufoniger, Iridomyrmex.bicknelli,Camponotus.consobrinus,
                            Monomorium.rothsteini, Iridomyrmex.suchieri,Aphaenogaster.longiceps,
                            Iridomyrmex.purpureus,Camponotus.claripes,Tetramorium.sp..A,Meranoplus.sp..A,
                            Monomorium.sydneyense) ~
                    Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
            family = poissonff, data = comb, Rank = 1, df1.nl = 2, Bestof = 30 
)
cao_p_summary = summary(cao_p)
sum(round(cao_p_summary@misc$deviance.Bestof,0) == round(min(cao_p_summary@misc$deviance.Bestof), 0))
x11()
par(mfrow=c(3,4))
plot(cao_p, lcol = "blue", lwd = 2, ylim = c(-5, 5), xlab = "", ylab = "")
saveRDS(cao_p, "003_processed_data/cao/antTraits/cao3.RDS")


cqo_p_ut = cqo(formula = cbind(Iridomyrmex.rufoniger, Iridomyrmex.bicknelli,Camponotus.consobrinus,
                            Monomorium.rothsteini, Aphaenogaster.longiceps,
                            Iridomyrmex.purpureus,Tetramorium.sp..A,Meranoplus.sp..A
                            ) ~
                    Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
            family = poissonff, data = comb, Rank = 1,  Bestof = 20, eq.toler = F
)
cqo_p_summary = summary(cqo_p)
sum(round(cqo_p_summary@misc$deviance.Bestof,0) == round(min(cqo_p_summary@misc$deviance.Bestof), 0))
saveRDS(cqo_p_ut, "003_processed_data/cqo/antTraits/cqo_p_ut.RDS")

Tol(cqo_p)[1, 1, ]
is.bell(cqo_p)

cqo2_p_ut = cqo(formula = cbind(Iridomyrmex.rufoniger, Iridomyrmex.bicknelli,Camponotus.consobrinus,
                               Monomorium.rothsteini, Aphaenogaster.longiceps,
                               Iridomyrmex.purpureus,Tetramorium.sp..A,Meranoplus.sp..A
) ~
        Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
family = poissonff, data = comb, Rank = 2,  Bestof = 30, eq.toler = F
)
cqo2_p_ut_summary = summary(cqo2_p_ut)
sum(round(cqo2_p_ut_summary@misc$deviance.Bestof,0) == round(min(cqo2_p_ut_summary@misc$deviance.Bestof), 0))
saveRDS(cqo2_p_ut, "003_processed_data/cqo/antTraits/cqo2_p_ut.RDS")

Tol(cqo2_p_ut)
is.bell(cqo2_p_ut)

AIC(cqo_p_ut)
AIC(cqo2_p_ut)


cqo_resid_plot(list(cqo_p_ut,
                    cqo2_p_ut))

cqo_nb_ut = cqo(formula = cbind(Iridomyrmex.rufoniger, Iridomyrmex.bicknelli,Camponotus.consobrinus,
                               Monomorium.rothsteini, Aphaenogaster.longiceps,
                               Iridomyrmex.purpureus,Tetramorium.sp..A,Meranoplus.sp..A
) ~
        Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
family = negbinomial, data = comb, Rank = 1,  Bestof = 40, eq.toler = F
)
cqo_nb_ut_summary = summary(cqo_nb_ut)
sum(round(cqo_nb_ut_summary@misc$deviance.Bestof,0) == round(min(cqo_nb_ut_summary@misc$deviance.Bestof), 0))
saveRDS(cqo_nb_ut, "003_processed_data/cqo/antTraits/cqo_nb_ut.RDS")

Tol(cqo_nb_ut_summary)
is.bell(cqo_nb_ut_summary)

cqo_resid_plot(list(cqo_nb_ut, cqo_p_ut))


# lets look at the plots 
par(mfrow=c(1,1))
lvplot(cqo_p_ut, label = T)
persp(cqo_p_ut, label = T)
lvplot(cqo2_p_ut, label = T)
persp(cqo2_p_ut, label = T)

# other tolerance assumptions 
cqo2_p_et = cqo(formula = cbind(Iridomyrmex.rufoniger, Iridomyrmex.bicknelli,Camponotus.consobrinus,
                                Monomorium.rothsteini, Aphaenogaster.longiceps,
                                Iridomyrmex.purpureus,Tetramorium.sp..A,Meranoplus.sp..A
) ~
        Bare.ground + Canopy.cover + Shrub.cover + Volume.lying.CWD + Feral.mammal.dung,
family = poissonff, data = comb, Rank = 2,  Bestof = 50, eq.toler = T
)
cqo2_p_et_summary = summary(cqo2_p_et)
sum(round(cqo2_p_et_summary@misc$deviance.Bestof,0) == round(min(cqo2_p_et_summary@misc$deviance.Bestof), 0))
saveRDS(cqo2_p_et, "003_processed_data/cqo/antTraits/cqo2_p_et.RDS")

# trivial 
is.bell(cqo2_p_et)

cqo_resid_plot(list(cqo2_p_et, cqo2_p_ut))
lvplot(cqo2_p_et, label = T)
persp(cqo2_p_et, label = T)

concoef(cqo2_p_et)
