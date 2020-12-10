# pre blog cqo spider 

pacman::p_load(data.table, dplyr, ggplot2, purrr, VGAM)
source("002_r_scripts/cqo/cqo_residual_plot.R")
data(hspider)
# standardize data 

hspider[, 1:6] <-
        scale(hspider[, 1:6]) # Standardized environmental variables

# first check unimodality 
set.seed = 1


# In Yee (2015, p.241), Thomas Yee suggests to keep S <= 10, n <=500 and p <=5. We will stick to this. 
# We have 12 species and will remove the two rarest ones. 
colSums(hspider[,7:18]) # Arctlute and Arctperi
# We will also remove a predictor variable. From the manyglms we know that CoveMoss recieved the least support on the community level. 
# A df1.nl of 1.5 approximates CQO. By setting df1.nl to 1 we are conservative so that species will rather not look unimodal. 

cao_p = cao(formula = cbind(Alopacce,Alopcune,Alopfabr,Auloalbi,
                    Pardlugu,Pardmont,Pardnigr,Pardpull,Trocterr,Zoraspin) ~
            WaterCon + BareSand + FallTwig + CoveHerb + ReflLux,
    family = poissonff, data = hspider, Rank = 1, df1.nl = 2, Bestof = 50 
)


# how often was the lowest value reached? Should at least be 3 to ensure the minimum is not local
# If it is lower than 3 increase Bestof. 
cao_p_summary = summary(cao_p)
sum(round(cao_p_summary@misc$deviance.Bestof,0) == round(min(cao_p_summary@misc$deviance.Bestof), 0))

par(mfrow=c(3,3))
plot(cao_p, lcol = "blue", lwd = 2, ylim = c(-5, 5), xlab = "", ylab = "")

# Alopfabr and Pardlugu dont seem unimodel so we wont model them from here on 

# First we try a unequal tolerance model. 
cqo_p = cqo(formula = cbind(Alopacce,Alopcune,Auloalbi,Pardmont,Pardnigr,Pardpull,Trocterr,Zoraspin) ~
                    WaterCon + BareSand + FallTwig + CoveHerb + ReflLux,
            family = poissonff, data = hspider, Rank = 1,  Bestof = 10, eq.toler = F
)
cqo_p_summary = summary(cqo_p)
sum(round(cqo_p_summary@misc$deviance.Bestof,0) == round(min(cqo_p_summary@misc$deviance.Bestof), 0))
# negative scores here would suggest that fitting a unimodal response has failed
Tol(cqo_p)[1, 1, ]
is.bell(cqo_p)
# next we try an additional latent variable 
cqo_p2 = cqo(formula = cbind(Alopacce,Alopcune,Auloalbi,Pardmont,Pardnigr,Pardpull,Trocterr,Zoraspin) ~
                    WaterCon + BareSand + FallTwig + CoveHerb + ReflLux,
            family = poissonff, data = hspider, Rank = 2,  Bestof = 20, eq.toler = F
)
cqo_p2_summary = summary(cqo_p2)
sum(round(cqo_p2_summary@misc$deviance.Bestof,0) == round(min(cqo_p2_summary@misc$deviance.Bestof), 0))

# compare rank 1 with rank 2
AIC(cqo_p)
AIC(cqo_p2)
# rank to fits better 
# lets have a look at the residuals 
par(mfrow = c(2,5))
plot(cqo_p)
cqo_resid_plot(list(cqo_p, cqo_p2), legend = T)
is.bell(cqo_p2)
Tol(cqo_p2)

# They are not bell shaped if there are negative entries on the diagonal 
par(mfrow=c(1,1))
lvplot(cqo_p)
lvplot(cqo_p2, label = T)
persp(cqo_p2, label = T)

# lets see how an equal tolerance model performs 
cqo_p = cqo(formula = cbind(Alopacce,Alopcune,Auloalbi,Pardmont,Pardnigr,Pardpull,Trocterr,Zoraspin) ~
                    WaterCon + BareSand + FallTwig + CoveHerb + ReflLux,
            family = poissonff, data = hspider, Rank = 1,  Bestof = 10, eq.toler = T, I.Toler = F
)
cqo_p_summary = summary(cqo_p)
sum(round(cqo_p_summary@misc$deviance.Bestof,0) == round(min(cqo_p_summary@misc$deviance.Bestof), 0))
# negative scores here would suggest that fitting a unimodal response has failed
Tol(cqo_p)[1, 1, ]
is.bell(cqo_p)
# next we try an additional latent variable 
cqo_p2 = cqo(formula = cbind(Alopacce,Alopcune,Auloalbi,Pardmont,Pardnigr,Pardpull,Trocterr,Zoraspin) ~
                     WaterCon + BareSand + FallTwig + CoveHerb + ReflLux,
             family = poissonff, data = hspider, Rank = 2,  Bestof = 50, eq.toler = T, I.Toler = F
)
cqo_p2_summary = summary(cqo_p2)
sum(round(cqo_p2_summary@misc$deviance.Bestof,0) == round(min(cqo_p2_summary@misc$deviance.Bestof), 0))
Tol(cqo_p2)
is.bell(cqo_p)

AIC(cqo_p)
AIC(cqo_p2)

cqo_resid_plot(list(cqo_p, cqo_p2), legend = T)

lvplot(cqo_p)
persp(cqo_p)
lvplot(cqo_p2, label = T)
persp(cqo_p2, label = T)


# I toler models 

cqo_p2 = cqo(formula = cbind(Alopacce,Alopcune,Auloalbi,Pardmont,Pardnigr,Pardpull,Trocterr,Zoraspin) ~
                     WaterCon + BareSand + FallTwig + CoveHerb + ReflLux,
             family = poissonff, data = hspider, Rank = 2,  Bestof = 50, eq.toler = T, I.Toler = T
)

# also doesn't converge 

coef(cqo_p)
