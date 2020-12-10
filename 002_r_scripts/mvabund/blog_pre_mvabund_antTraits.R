# blog pre mvabund ant traits 

pacman::p_load(data.table, dplyr, gllvm, ggplot2, lattice, magrittr, mvabund, stringr, readxl)

data(antTraits)

antTraits$abund
antTraits$env
antTraits$trait


ant_mv <- mvabund(antTraits$abund)

plot(ant_mv)

# quadratic 
meanvar.plot(dune_mvabund,xlab=expression(mu), ylab = expression(sigma))

model_list = list()
cou = 1
for (i in 1:2){
        family = switch(i,"poisson","negative.binomial")
        for (k in 1:2){
                cor.var = switch(k, "I", "shrink")
                model_list[[cou]] = manyglm(ant_mv ~ ., data = antTraits$env, family = family, cor.type = cor.var)
                names(model_list)[cou] = paste0(family,"_",cor.var)
                cou = cou + 1
        }
}

names(model_list)

# first we compare nb vs poisson than cor.type
# Residual plots 
plot.manyglm(model_list[[1]], which = 1) # slight fan shape 
plot.manyglm(model_list[[3]], which = 1) # no fan shape 
plot.manyglm(model_list[[1]], which = 2) # bad 
plot.manyglm(model_list[[3]], which = 2) # ok 
plot.manyglm(model_list[[1]], which = 3) # slopes up with increasing linear predictor  
plot.manyglm(model_list[[3]], which = 3) # ok
# -> all suggest nb model 
# now cor.types 
plot.manyglm(model_list[[3]], which = 1)
plot.manyglm(model_list[[4]], which = 1)
plot.manyglm(model_list[[3]], which = 2)
plot.manyglm(model_list[[4]], which = 2)
plot.manyglm(model_list[[3]], which = 3)
plot.manyglm(model_list[[4]], which = 3)
# -> very similar, lets compare results 

lapply(X = model_list, FUN = AIC) %>% 
        unlist -> aic_list
aic_table = data.table(model = names(aic_list), 
                       aic = aic_list)
aic_table[, family := unlist(lapply(str_split(model, pattern = "_"), function(x)x[1]))]
aic_table[, cor.type := unlist(lapply(str_split(model, pattern = "_"), function(x)x[2]))]
aic_table[, cor.type := str_remove(cor.type, pattern = "[0-9].*$")]
aic_table[, taxon := rep(names(antTraits$abund), times = 4)]

aic_table %>% ggplot(aes(x = family, y = aic, col = taxon, group = taxon)) + 
        geom_point() + 
        geom_line() + 
        geom_text(data = filter(aic_table, family == "poisson"), aes(label=taxon), nudge_x = .2) + 
        theme(legend.position = "none")

# For most species NB AIC is lower than AIC 
# 2 sec for 10 -> 999 = 200 sec
test_mod = anova(model_list[[3]], nBoot = 10)

# We will also remove a predictor 
anova_list = list()
c = 1
for (j in 1:3) {
        test_st = switch(j, "LR", "wald", "score")
        for (l in 1:4) {
                resamp_method = switch(l,
                                       "case",
                                       "perm.resid",
                                       "montecarlo",
                                       "pit.trap")
                if (test_st == "score" &
                    resamp_method == "case")
                        next()
                anova_list[[c]] = anova.manyglm(
                        object = model_list[[3]],
                        p.uni = "adjusted",
                        test = test_st,
                        resamp = resamp_method
                )
                names(anova_list)[[c]] = paste(test_st,
                                               resamp_method,
                                               sep = "_")
                print(c)
                c = c + 1
        }
        
}

saveRDS(anova_list,  "003_processed_data/mvabund/anova/antTrait//anova_list.RDS")
anova_list = readRDS("003_processed_data/mvabund/anova/antTrait/anova_list.RDS")

plot_list = list()
c = 1
for (i in seq_along(anova_list)){
        loop_set = anova_list[[i]]
        split_name = str_split_fixed(names(anova_list)[i], "_", n = 2)
        loop_data = data.frame(test = split_name[,1],
                               resamp = split_name[,2],
                               p.value = loop_set$table[2:6,4],
                               variable = rownames(loop_set$table)[2:6])
        plot_list[[c]] = loop_data
        c = c + 1
}
plot_dt = rbindlist(plot_list)
x11()
ggplot(plot_dt , aes(x = resamp, y = p.value, col = test)) + 
        geom_jitter(size=2) + 
        facet_wrap(.~variable) + 
        geom_hline(yintercept = 0.05)

# canopy cover receives the lowest p-values across the board. Test statistics are similar, wald tends towards lower values. 
# Monte carlo generally receives the lowest p-values  

plot_data_species = anova_list$LR_pit.trap$uni.p[-1,]
plot_data_species = data.frame(plot_data_species)
plot_data_species$variable = rownames(plot_data_species)
plot_data_species %<>% tidyr::pivot_longer(cols = names(plot_data_species)[-which(names(plot_data_species) == "variable")])
plot_data_species %>% 
        ggplot(aes(x = variable, y = value)) + 
        geom_point(aes(col = name)) + 
        geom_hline(yintercept = 0.05) 

# most species seem to be distributed independently of the measured variables

# some final plots 
# these plots do not change with cor.type so only I is shown 
a <- max(abs(coef(model_list[[3]])))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(model_list[[3]]))), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
x11()
print(plot.tas)

# remove dominant variable 
a <- max(abs(coef(model_list[[1]])[-5,]))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(model_list[[1]])[-5,])), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)

# Modeling traits 

trait_model1=traitglm(antTraits$abund,antTraits$env,antTraits$traits,method="glm1path", family = "negative.binomial")
trait_model1$fourth.corner
a        = max( abs(trait_model1$fourth.corner) )
colort   = colorRampPalette(c("blue","white","red")) 
plot.4th = levelplot(t(as.matrix(trait_model1$fourth.corner)), xlab="Environmental Variables",
                     ylab="Species traits", col.regions=colort(100), at=seq(-a, a, length=100),
                     scales = list( x= list(rot = 45)))
print(plot.4th)
