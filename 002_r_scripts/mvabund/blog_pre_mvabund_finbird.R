# blog pre mvabund finish birds 

pacman::p_load(data.table, dplyr, ggplot2, lattice, magrittr, mvabund, stringr, readxl)

data = fread("001_raw_data/hmsc_birds/data/data.csv")
traits = fread("001_raw_data/hmsc_birds/data/traits.csv")
env = data[,1:9]
data = data[,-c(1:9)]

env[,c("Year", "Habitat") := .(factor(Year), factor(Habitat))]

data_mv <- mvabund(data)

plot(data_mv)
plot(data_mv~env$Year)
plot(data_mv~env$Habitat)

# quadratic 
meanvar.plot(data_mv,xlab=expression(mu), ylab = expression(sigma))
meanvar.plot(data_mv~env$Year,xlab=expression(mu), ylab = expression(sigma))
meanvar.plot(data_mv~env$Habitat,xlab=expression(mu), ylab = expression(sigma))

model_list = list()
cou = 1
for (i in 1:2){
        family = switch(i,"poisson","negative.binomial")
        for (k in 1:2){
                cor.var = switch(k, "I", "shrink")
                model_list[[cou]] = manyglm(data_mv ~ ., data = env, family = family, cor.type = cor.var)
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
aic_table[, taxon := rep(names(data), times = 4)]

aic_table %>% ggplot(aes(x = family, y = aic, col = taxon, group = taxon)) + 
        geom_point() + 
        geom_line() + 
        geom_text(data = filter(aic_table, family == "poisson"), aes(label=taxon), nudge_x = .2) + 
        theme(legend.position = "none")

# For most species NB AIC is lower than AIC


#29 sec per bootstrap 
anova_obj = anova.manyglm(model_list[[3]], nBoot = 200, p.uni = "adjusted") 

saveRDS(anova_list,  "003_processed_data/mvabund/anova/finbird/anova_list.RDS")
anova_obj = readRDS("003_processed_data/mvabund/anova/finbird/anova_list.RDS")

anova_obj

my_color_palette <- c("#7fc97f","#d95f02","#1b9e77","#666666","#bf5b17","#5f64ff","#ff9a14","#dcce00","#03eaff","#e6ab02","#66a61e","#e7298a","#7570b3","#ff00bf","#00fe04","#a6cee3","#a6761d","#386cb0","#fdc086","#beaed4")
plot_data_species = anova_obj$uni.p
plot_data_species = data.frame(plot_data_species)
plot_data_species$variable = rownames(plot_data_species)
plot_data_species %<>% tidyr::pivot_longer(cols = names(plot_data_species)[-which(names(plot_data_species) == "variable")])
plot_data_species$variable %<>% factor
plot_data_species %>% 
        ggplot(aes(x = value, y = name)) + 
        geom_point(aes(col = variable)) + 
        geom_vline(xintercept = 0.05) + 
        scale_color_manual(values = my_color_palette)
plot_data_species %>% 
        ggplot(aes(y = value, x = variable)) +
        geom_boxplot()


# some final plots 
# these plots do not change with cor.type so only I is shown 
a <- max(abs(coef(model_list[[3]])))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(model_list[[3]]))), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)

a <- max(abs(coef(model_list[[3]])[-1,]))
colort <- colorRampPalette(c("blue", "white", "red"))
plot.tas <- levelplot(t(as.matrix(coef(model_list[[3]])[-1,])), ylab = "", xlab  = "", col.regions = colort(100), at=seq(-a,a,length = 100), scales = list( x= list(rot = 45)))
print(plot.tas)


