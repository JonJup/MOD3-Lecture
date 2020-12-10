# blog pre mvabund fish 

pacman::p_load(data.table, dplyr, ggplot2, lattice, magrittr, mvabund, stringr, readxl)

data = read_xlsx("001_raw_data/anderson copula/ece34948-sup-0001-tables1.xlsx", skip = 1)

pkf_mva = mvabund(data[,-c(1,2)])

season = factor(pull(data[,2]), levels = c("Sep.98", "Mar.99", "Sep.99"))

plot(pkf_mva~season)

# quadratic 
meanvar.plot(pkf_mva,xlab=expression(mu), ylab = expression(sigma))
meanvar.plot(pkf_mva~season,xlab=expression(mu), ylab = expression(sigma))

model_list = list()
cou = 1
for (i in 1:2){
        family = switch(i,"poisson","negative.binomial")
        for (k in 1:2){
                cor.var = switch(k, "I", "shrink")
                model_list[[cou]] = manyglm(pkf_mva ~ season, family = family, cor.type = cor.var)
                names(model_list)[cou] = paste0(family,"_",cor.var)
                cou = cou + 1
        }
}

names(model_list)

# first we compare nb vs poisson than cor.type
# Residual plots 
plot.manyglm(model_list[[1]], which = 1) # fan shape 
plot.manyglm(model_list[[3]], which = 1) # no fan shape 
plot.manyglm(model_list[[1]], which = 2) # bad 
plot.manyglm(model_list[[3]], which = 2) # ok 
plot.manyglm(model_list[[1]], which = 3) # slopes up with increasing linear predictor  
plot.manyglm(model_list[[3]], which = 3) # god   
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
aic_table[, taxon := rep(names(data)[-c(1,2)], times = 4)]

aic_table %>% ggplot(aes(x = family, y = aic, col = taxon, group = taxon)) + 
        geom_point() + 
        geom_line() + 
        geom_text(data = filter(aic_table, family == "poisson"), aes(label=taxon), nudge_x = .2) + 
        theme(legend.position = "none")

# For most species NB AIC is lower than Poisson

# 5 seconds for 100 bootstraps 
test_mod1 = anova.manyglm(model_list[[4]], nBoot = 100) 

anova_list_I = list()
anova_list_shrink = list()
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
                anova_list_I[[c]] = anova.manyglm(
                        object = model_list[[3]],
                        p.uni = "adjusted",
                        test = test_st,
                        resamp = resamp_method
                )
                anova_list_shrink[[c]] = anova.manyglm(
                        object = model_list[[4]],
                        p.uni = "adjusted",
                        test = test_st,
                        resamp = resamp_method
                )
                names(anova_list_I)[[c]] = paste(test_st,
                                               resamp_method,
                                               sep = "_")
                names(anova_list_shrink)[[c]] = paste(test_st,
                                               resamp_method,
                                               sep = "_")
                print(c)
                c = c + 1
        }
        
}

saveRDS(anova_list_I,  "003_processed_data/mvabund/anova/poorknight/anova_list_I.RDS")
saveRDS(anova_list_shrink,  "003_processed_data/mvabund/anova/poorknight/anova_list_shrink.RDS")
anova_list = readRDS("003_processed_data/mvabund/anova/poorknight/anova_list.RDS")


plot_list_I = list()
c = 1
for (i in seq_along(anova_list_I)){
        loop_set = anova_list_I[[i]]
        split_name = str_split_fixed(names(anova_list_I)[i], "_", n = 2)
        loop_data = data.frame(
                               test = split_name[,1],
                               resamp = split_name[,2],
                               p.value = loop_set$table[2,4],
                               variable = rownames(loop_set$table)[2])
        plot_list_I[[c]] = loop_data
        c = c + 1
}
plot_dt_I = rbindlist(plot_list_I)

plot_list_shrink = list()
c = 1
for (i in seq_along(anova_list_shrink)){
        loop_set = anova_list_shrink[[i]]
        split_name = str_split_fixed(names(anova_list_shrink)[i], "_", n = 2)
        loop_data = data.frame(
                test = split_name[,1],
                resamp = split_name[,2],
                p.value = loop_set$table[2,4],
                variable = rownames(loop_set$table)[2])
        plot_list_shrink[[c]] = loop_data
        c = c + 1
}
plot_dt_shrink = rbindlist(plot_list_shrink)

ggplot(plot_dt_I , aes(x = resamp, y = p.value, col = test)) + 
        geom_point(size=2) + 
        facet_wrap(.~variable) + 
        geom_hline(yintercept = 0.05)
ggplot(plot_dt_shrink , aes(x = resamp, y = p.value, col = test)) + 
        geom_point(size=2) + 
        facet_wrap(.~variable) + 
        geom_hline(yintercept = 0.05)


p_uni_fun = function(x){
        plot_data_species = x 
        plot_data_species = data.frame(plot_data_species)
        plot_data_species$variable = rownames(plot_data_species)
        plot_data_species %<>% tidyr::pivot_longer(cols = "plot_data_species")
        plot_data_species$name = "season"
        plot_data_species %>% 
                ggplot(aes(y = variable, x = value)) + 
                geom_point(aes(col = name)) + 
                geom_vline(xintercept = 0.05) 
}

p_uni_fun(anova_list_I$wald_pit.trap$uni.p[-1,])
p_uni_fun(anova_list_shrink$wald_pit.trap$uni.p[-1,])
p_uni_fun(anova_list_I$LR_pit.trap$uni.p[-1,])
p_uni_fun(anova_list_shrink$LR_pit.trap$uni.p[-1,])

p_uni_fun(anova_list$LR_pit.trap$uni.p[-1,])
p_uni_fun(anova_list$score_pit.trap$uni.p[-1,])


# some final plots 
# these plots do not change with cor.type so only I is shown 
a <- max(abs(coef(model_list[[3]])))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(model_list[[3]]))), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)

a <- max(abs(coef(model_list[[4]])))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(model_list[[4]]))), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)

