### --- MOD3
### --- BLOG
### --- MvAbund 
### --- Tasmania Data 


## -- Overview 
# Copepode simple 
# Nematode simple 
# Copepod complex 
# Nematode complex 


# setup -------------------------------------------------------------------
library(mvabund)
library(lattice)
library(stringr)
library(data.table)
library(ggplot2)
library(dplyr)

data("Tasmania")
attach(Tasmania)

Tasmania$tr.block

head(copepods[,1:3])
cop_mvabund <- mvabund(copepods)
nem_mvabund <- mvabund(nematodes)

plot(cop_mvabund~treatment, col = as.numeric(block))
plot(nem_mvabund~treatment, col = as.numeric(block))

meanvar.plot(cop_mvabund~tr.block, col = as.numeric (treatment))
meanvar.plot(nem_mvabund~tr.block, col = as.numeric (treatment))

# lets fit gaussian, poisson and negative binomial data to all 
models_cop = list()
models_nem = list()
cou = 1
for (i in 1:2){
        family = switch(i,"poisson","negative.binomial")
        for (k in 1:2){
               cor.var = switch(k, "I", "shrink")
               models_cop[[cou]] = manyglm(cop_mvabund ~ block*treatment, family = family, cor.type = cor.var)
               models_nem[[cou]] = manyglm(nem_mvabund ~ block*treatment, family = family, cor.type = cor.var)
               names(models_cop)[cou] = paste0(family,"_",cor.var)
               names(models_nem)[cou] = paste0(family,"_",cor.var)
               cou = cou + 1
        }
}

names(models_cop)

# first we compare nb vs poisson than cor.type
# Residual plots for copepods 
par(mfrow = c(3,2))
plot.manyglm(models_cop[[1]], which = 1) # fan shape
plot.manyglm(models_cop[[3]], which = 1)
plot.manyglm(models_cop[[1]], which = 2)
plot.manyglm(models_cop[[3]], which = 2) # better fit 
plot.manyglm(models_cop[[1]], which = 3)
plot.manyglm(models_cop[[3]], which = 3) # better fit 
# -> all suggest nb model 
# now cor.types 
plot.manyglm(models_cop[[3]], which = 1)
plot.manyglm(models_cop[[4]], which = 1)
plot.manyglm(models_cop[[3]], which = 2)
plot.manyglm(models_cop[[4]], which = 2)
plot.manyglm(models_cop[[3]], which = 3)
plot.manyglm(models_cop[[4]], which = 3)

# -> very similar, lets compare results 
plot.manyglm(models_nem[[1]], which = 1) # fan shape
plot.manyglm(models_nem[[3]], which = 1)
plot.manyglm(models_nem[[1]], which = 2)
plot.manyglm(models_nem[[3]], which = 2) # better fit 
plot.manyglm(models_nem[[1]], which = 3)
plot.manyglm(models_nem[[3]], which = 3) # better fit 
# -> all suggest nb model 
# now cor.types 
plot.manyglm(models_nem[[3]], which = 1)
plot.manyglm(models_nem[[4]], which = 1)
plot.manyglm(models_nem[[3]], which = 2)
plot.manyglm(models_nem[[4]], which = 2)
plot.manyglm(models_nem[[3]], which = 3)
plot.manyglm(models_nem[[4]], which = 3)
anova_list = list()
c = 1

# AIC 
lapply(X = model_list, FUN = AIC) %>% 
        unlist -> aic_list
aic_table = data.table(model = names(aic_list), 
                       aic = aic_list)
aic_table[, family := unlist(lapply(str_split(model, pattern = "_"), function(x)x[1]))]
aic_table[, cor.type := unlist(lapply(str_split(model, pattern = "_"), function(x)x[2]))]
aic_table[, cor.type := str_remove(cor.type, pattern = "[0-9].*$")]
aic_table[, taxon := rep(names(dune), times = 4)]

aic_table %>% ggplot(aes(x = family, y = aic, col = taxon, group = taxon)) + 
        geom_point() + 
        geom_line() + 
        geom_text(data = filter(aic_table, family == "poisson"), aes(label=taxon), nudge_x = .2) + 
        theme(legend.position = "none")



for (i in 1:2){
        data = switch(i,"cop", "nem")
        for (k in 1:2) {
                cor = switch(k, 3,4)
                test_object = get(paste0("models_",data))
                test_object = test_object[[cor]]
                cor_save_name = ifelse(cor == 3, "I", "shrink")
                for(j in 1:3){
                        test_st = switch(j, "LR", "wald", "score")
                        if (cor == 4 & test_st == "LR") next()
                        for (l in 1:4){
                                resamp_method = switch(l, "case", "perm.resid", "montecarlo", "pit.trap")
                                if (test_st == "score" & resamp_method == "case") next()
                                anova_list[[c]] = anova.manyglm(object = test_object,
                                                                p.uni = "adjusted",
                                                                test = test_st,
                                                                resamp = resamp_method)
                                names(anova_list)[[c]] = paste(data, cor_save_name,test_st,resamp_method, sep = "_")
                                c = c + 1
                        }

                }

                anova_list[[c]] = anova.manyglm(object = test_object)
                c = c + 1
        }
}
# saveRDS(anova_list, "003_processed_data/mvabund/anova/tasmania/anova_list.RDS")
anova_list = readRDS("003_processed_data/mvabund/anova/tasmania/anova_list.RDS")
# show community p-values for nematodes for different methods 
nem_id = which(str_detect(names(anova_list), "nem"))
cop_id = which(str_detect(names(anova_list), "cop"))

plot_list = list()
c = 1
for (i in nem_id){
        loop_set = anova_list[[i]]
        split_name = str_split_fixed(names(anova_list)[i], "_", n = 4)
        loop_data = data.frame(cor = split_name[,2],
                               test = split_name[,3],
                               resamp = split_name[,4],
                               p.value = loop_set$table[2:4,4],
                               variable = c("block", "treatment", "interaction"))
        plot_list[[c]] = loop_data
        c = c + 1
}
plot_dt = rbindlist(plot_list)
ggplot(plot_dt , aes(x = resamp, y = p.value, col = test, shape = cor)) + geom_jitter(size=2) + facet_wrap(.~variable) + geom_hline(yintercept = 0.05)

plot_list = list()
c = 1
for (i in cop_id){
        loop_set = anova_list[[i]]
        split_name = str_split_fixed(names(anova_list)[i], "_", n = 4)
        loop_data = data.frame(cor = split_name[,2],
                               test = split_name[,3],
                               resamp = split_name[,4],
                               p.value = loop_set$table[2:4,4],
                               variable = c("block", "treatment", "interaction"))
        plot_list[[c]] = loop_data
        c = c + 1
}
plot_dt = rbindlist(plot_list)
ggplot(plot_dt , aes(x = resamp, y = p.value, col = test, shape = cor)) + geom_jitter(size=2) + facet_wrap(.~variable) + geom_hline(yintercept = 0.05)

# some final plots 

# these plots do not change with cor.type so only I is shown 

a <- max(abs(coef(models_cop[[3]])))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(models_cop[[3]]))), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)

a <- max(abs(coef(models_nem[[3]])))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(models_nem[[3]]))), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)


cop_pred = predict(models_cop[[3]], type = "response")
matplot(t(cop_pred[c(1,3),]), type = "l", xaxt = "n", log = "y", ylab = "Mean abundance [log]")
axis(1, at =1:12, labels = colnames(copepods), las = 3)
legend("topright", legend = levels(treatment), col = 1:2, lty = 1:2)



