pacman::p_load(data.table, dplyr, ggplot2, purrr, VGAM, mvabund)
source("002_r_scripts/cqo/cqo_residual_plot.R")

#load data 
data = fread("001_raw_data/hmsc_birds/data/data.csv")
env = data[,1:9]
data = data[,-c(1:9)]
env[,c("Year", "Habitat") := .(factor(Year), factor(Habitat))]

# choose species 
colSums(data) %>% 
        sort(decreasing=T) %>% 
        names %>% 
        .[1:11]

env = env[, c("x","y","Effort","JunJul","AprMay")]
env = scale(env)

comb = cbind(data,env)
nrow(comb)
sample_sites = sample(1:nrow(comb), 500)
comb_use = comb[sample_sites, ]
set.seed = 1
cao_p = cao(
        formula = cbind(Phylloscopus_trochilus, Fringilla_coelebs, Anthus_trivialis, Erithacus_rubecula,
                        Carduelis_spinus, Turdus_philomelos, Parus_major, Turdus_pilaris, Turdus_iliacus,
                        Emberiza_citrinella) ~
                x + y + Effort + JunJul + AprMay,
        family = poissonff,
        data = comb_use,
        Rank = 1,
        df1.nl = 2,
        Bestof = 10
)
cao_p_summary = summary(cao_p)
sum(round(cao_p_summary@misc$deviance.Bestof,0) == round(min(cao_p_summary@misc$deviance.Bestof), 0))
x11()
par(mfrow=c(3,4))
plot(cao_p, lcol = "blue", lwd = 2, ylim = c(-5, 5), xlab = "", ylab = "")
saveRDS(cao_p, "003_processed_data/cao/finbird/cao1.RDS")
set.seed = 1
cao_p = cao(
        formula = cbind(Phylloscopus_trochilus, Fringilla_coelebs, Anthus_trivialis, Erithacus_rubecula,
                        Carduelis_spinus, Turdus_philomelos, Parus_major, Turdus_pilaris, Columba_palumbus,
                        Emberiza_citrinella) ~
                x + y + Effort + JunJul + AprMay,
        family = poissonff,
        data = comb_use,
        Rank = 1,
        df1.nl = 2,
        Bestof = 10
)
cao_p_summary = summary(cao_p)
sum(round(cao_p_summary@misc$deviance.Bestof,0) == round(min(cao_p_summary@misc$deviance.Bestof), 0))
x11()
par(mfrow=c(3,4))
plot(cao_p, lcol = "blue", lwd = 2, ylim = c(-5, 5), xlab = "", ylab = "")
saveRDS(cao_p, "003_processed_data/cao/finbird/cao2.RDS")

# that ok 
cqo_p_ut = cqo(formula = cbind(Phylloscopus_trochilus, Fringilla_coelebs, Anthus_trivialis, Erithacus_rubecula,
                            Carduelis_spinus, Turdus_philomelos, Parus_major, Turdus_pilaris, Columba_palumbus,
                            Emberiza_citrinella) ~
                    x + y + Effort + JunJul + AprMay,
            family = poissonff, data = comb_use, Rank = 1, Bestof = 50, eq.toler = F
)
cqo_p_ut_sum = summary(cqo_p_ut)
sum(round(cqo_p_ut_sum@misc$deviance.Bestof,0) == round(min(cqo_p_ut_sum@misc$deviance.Bestof), 0))
saveRDS(cqo_p_ut, "003_processed_data/cqo/finbird/cqo_p_ut.RDS")
Tol(cqo_p_ut)[1, 1, ]
is.bell(cqo_p_ut)

concoef(cqo_p_ut)
lvplot(cqo_p_ut, label = T)
persp(cqo_p_ut, label = T)
cqo_resid_plot(cqo_p_ut)

cqo2_p_et = cqo(formula = cbind(Phylloscopus_trochilus, Fringilla_coelebs, Anthus_trivialis, Erithacus_rubecula,
                               Carduelis_spinus, Turdus_philomelos, Parus_major, Turdus_pilaris, Columba_palumbus,
                               Emberiza_citrinella) ~
                       x + y + Effort + JunJul + AprMay,
               family = poissonff, data = comb_use, Rank = 2,  Bestof = 100, eq.toler = F
)
cqo2_p_ut_sum = summary(cqo2_p_ut)
sum(round(cqo_p_ut_sum@misc$deviance.Bestof,0) == round(min(cqo_p_ut_sum@misc$deviance.Bestof), 0))
saveRDS(cqo_p_ut, "003_processed_data/cqo/finbird/cqo2_p_ut.RDS")
Tol(cqo2_p_ut)[1, 1, ]
is.bell(cqo2_p_ut)

lvplot(cqo2_p_ut, label = T)
persp(cqo2_p_ut, label = T)
cqo_resid_plot(cqo2_p_ut)

