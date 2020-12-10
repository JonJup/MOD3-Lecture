# copula with finbirds

pacman::p_load(corrplot, ecoCopula, here, magrittr, readxl, data.table,dplyr, stringr, vegan, mvabund)
# load functions 
source("../../../Supplementary Informations/Anderson19Copula/RCode1.R")
source("../../../Supplementary Informations/Anderson19Copula/RCode2.R")
source("../../../Supplementary Informations/Anderson19Copula/RCode3_functions.R")
data = fread("001_raw_data/hmsc_birds/data/data.csv")

x = data[,1:9]
y = data[,-c(1:9)]
x[,c("Year", "Habitat") := .(factor(Year), factor(Habitat))]
rm(data)
levels(x$Habitat)
table(x$Habitat)
# Br: Broad leaved forest 
# Co: Coniferous forest 
# Op: Open habitat 
# Urb: Urban habitat 
# We: Wetlands 
# separate data based on land use and convert to matrix. The second step is necessary for the chooseDistr() function

ls_y = list("y_op" = y[x$Habitat == "Op",] %>% as.matrix(), 
            "y_ub" = y[x$Habitat == "Urb", ] %>% as.matrix(), 
            "y_br" = y[x$Habitat == "Br",] %>% as.matrix()
            )

# any species missing completely from one of the data sets? 
any(unlist(lapply(ls_y, colSums))==0)

# find marginal distributions 
ls_mar = lapply(ls_y, function(z)chooseDistr(z, X = x[])$marginals)
ls_mar = list()
ls_mar$y_op = chooseDistr(Y = ls_y$y_op, X = x[Habitat == "Op", x])
ls_mar$y_ub = chooseDistr(Y = ls_y$y_ub, X = x[Habitat == "Urb", x])
ls_mar$y_br = chooseDistr(Y = ls_y$y_br, X = x[Habitat == "Br", x])

saveRDS(ls_mar, "003_processed_data/copula/findbird/ls_mar_x.RDS")

# now we want to find associations between species. There are multiple ways to do this and 
# below I will look at three different ones: Index of associations - the original methods from Anderson et al 2019; 
# gaussian graphical copula models (Popovic et al 2019) and gllvm residual correlations. 

# first of index of association 
ls_pwa = lapply(ls_y, function(x)
                pairWise(
                        x,
                        nperm = 99999,
                        alpha_type = "PCER",
                        graphic = FALSE,
                        sig_level = 0.01
                ))

saveRDS(ls_pwa, "003_processed_data/copula/findbird/ls_pwa_use.RDS")
ls_pwa = readRDS("003_processed_data/copula/findbird/ls_pwa_use.RDS")

# correlation plots 
par(mfrow=c(1,1))
corrplot(ls_pwa$y_op$IoA.obs, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_ub$IoA.obs, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_br$IoA.obs, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_op$IoA.shrunk, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_ub$IoA.shrunk, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_br$IoA.shrunk, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_op$IoA.subset, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_ub$IoA.subset, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
corrplot(ls_pwa$y_br$IoA.subset, diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")

## -- I wont run the gllvm again here but simply refer to their results: 

# Ecocopula ---------------------------------------------------------------
ls_mvabund = ls()
ls_mvabund = lapply(ls_y, function(x)mvabund(x))
ls_mvabund$y_hp_mod_pois = manyglm(ls_y$y_hp ~ ., data = x[x$Use == "Haypastu",], family = "poisson")
par(mfrow=c(2,2))
for (i in 1:3) plot(y_0_model_p, which = i)

ls_mvabund$y_pa_mod_pois = manyglm(ls_y$y_hp ~ ., data = x[x$Use == "Pasture",], family = "poisson")
ls_mvabund$y_hf_mod_pois = manyglm(ls_y$y_hp ~ ., data = x[x$Use == "Hayfield",], family = "poisson")



y_0_model_p = manyglm(y_0 ~ ., data = x[x$sand_factor == 0,], family = "negbinom")
y0_gcgm_2 <- cgr(y_0_model_p,lambda=.2)
plot(y0_gcgm_2)
y0_ec_accepted = colSums(y0_gcgm_2$best_graph$graph) > 1

y1_mva = mvabund(y_1)
y1_model_p = manyglm(y_1 ~ ., data = x[x$sand_factor == 1,], family = "poisson")
for (i in 1:3) plot(y1_model_p, which = i)
y1_model_nb = manyglm(y_1 ~ ., data = x[x$sand_factor == 1,], family = "negbinom")
for (i in 1:3) plot(y1_model_nb, which = i)
y1_gcgm_2 <- cgr(y1_model_nb,lambda=.2)
plot(y1_gcgm_2)
y1_ec_accepted = colSums(y1_gcgm_2$best_graph$graph) > 1


y0_ec_accepted = names(y0_ec_accepted)[which(y0_ec_accepted)]
y1_ec_accepted = names(y1_ec_accepted)[which(y1_ec_accepted)]


# add species column. Later functions require this column.
y_0_mar$Species = rownames(y_0_mar)
y_1_mar$Species = rownames(y_1_mar)

# ids of significant pairwise associated taxa
y_op_pwa_id <- which(colnames(ls_y$y_op) %in%  ls_pwa$y_op$associated)
y_ub_pwa_id <- which(colnames(ls_y$y_ub) %in%  ls_pwa$y_ub$associated)
y_br_pwa_id <- which(colnames(ls_y$y_br) %in%  ls_pwa$y_br$associated)



y0_eco_id <- which(colnames(y_0) %in%  y0_ec_accepted)
y1_eco_id <- which(colnames(y_1) %in%  y1_ec_accepted)

# new data sets with associated species 
y0_assoc_pwa = y_0[,y0_pwa_id]
y0_assoc_eco = y_0[,y0_eco_id]
y1_assoc_pwa = y_1[,y1_pwa_id]
y1_assoc_eco = y_1[,y1_eco_id]

ls_y$y_op_pwa_assoc = ls_y$y_op[,y_op_pwa_id]

y0_assoc_pwa_marginal = y_0_mar[y0_pwa_id,]
y0_assoc_eco_marginal = y_0_mar[y0_eco_id,]
y1_assoc_pwa_marginal = y_1_mar[y1_pwa_id,]
y1_assoc_eco_marginal = y_1_mar[y1_eco_id,]

ls_copula = list()
ls_copula$pwa_op = estimate_copula(data = y0_assoc_pwa, marginal_details = y0_assoc_pwa_marginal)
y0_copula_pwa <- estimate_copula(data = y0_assoc_pwa, marginal_details = y0_assoc_pwa_marginal)
y1_copula_pwa <- estimate_copula(data = y1_assoc_pwa, marginal_details = y1_assoc_pwa_marginal)
y1_copula_pwa <- estimate_copula(data = y1_assoc_pwa, marginal_details = y1_assoc_pwa_marginal)
y0_copula_eco <- estimate_copula(data = y0_assoc_eco, marginal_details = y0_assoc_eco_marginal)
y1_copula_eco <- estimate_copula(data = y1_assoc_eco, marginal_details = y1_assoc_eco_marginal)

saveRDS(y0_copula_pwa, "003_processed_data/copula/spider/copula_sand0_pwa.RDS")
saveRDS(y1_copula_pwa, "003_processed_data/copula/spider/copula_sand1_pwa.RDS")
saveRDS(y0_copula_eco, "003_processed_data/copula/spider/copula_sand0_eco.RDS")
saveRDS(y1_copula_eco, "003_processed_data/copula/spider/copula_sand1_eco.RDS")

# Extract the MLE correlation matrix from the mcem_result list.
y0pwa_corr_mcem <- y0_copula_pwa[["cov_final"]]
y1pwa_corr_mcem <- y1_copula_pwa[["cov_final"]]
y0eco_corr_mcem <- y0_copula_eco[["cov_final"]]
y1eco_corr_mcem <- y1_copula_eco[["cov_final"]]

# We need to add the unassociated species back into the covariance matrix
add_unass <- function(data, spe_names) {
        new_names <- setdiff(spe_names, row.names(data))
        n_new     <- length(new_names)
        n_old     <- ncol(data)
        mt_add1   <- matrix(0, ncol = ncol(data), nrow = n_new)
        mt_add2   <- matrix(0, ncol = n_new, nrow = nrow(data))
        data <- rbind(data, mt_add1)
        row.names(data)[(n_old+1):(n_old+n_new)] <- new_names
        mt_add2   <- matrix(0, ncol = n_new, nrow = nrow(data))
        data <- cbind(data, mt_add2)
        colnames(data)[(n_old+1):(n_old+n_new)] <- new_names
        diag(data) <- 1
        return(data)
}

y0pwa_corr_mcem_add <- add_unass(data = y0pwa_corr_mcem, spe_names = colnames(y_0))
y1pwa_corr_mcem_add <- add_unass(data = y1pwa_corr_mcem, spe_names = colnames(y_1))
y0eco_corr_mcem_add <- add_unass(data = y0eco_corr_mcem, spe_names = colnames(y_0))
y1eco_corr_mcem_add <- add_unass(data = y1eco_corr_mcem, spe_names = colnames(y_1))

# evaluate quality of copula models 
par(mfrow=c(1,2))
N = 100
simulated_data_y0_pwa <- generate_copula_data(N, marginal_details = y_0_mar, cov = y0pwa_corr_mcem_add)
simulated_data_y0_eco <- generate_copula_data(N, marginal_details = y_0_mar, cov = y0eco_corr_mcem_add)
species_indices <- sample(x = colnames(y_0), 2)
pch <- c(rep(0, nrow(y_0)), rep(20, N))
plot(rbind(y_0[ , species_indices], simulated_data_y0_pwa[["observed"]][ , species_indices]), pch = pch, col = c(rep("red", nrow(y_0)), rep("black", N)))
plot(rbind(y_0[ , species_indices], simulated_data_y0_eco[["observed"]][ , species_indices]), pch = pch, col = c(rep("red", nrow(y_0)), rep("black", N)))
legend(x = "topright", pch = c(0, 20), legend = c("observed", "simulated"))



# 
n_sim = 100
n_g1 = nrow(y_0)
n_g2 = nrow(y_1)
group_var = append(rep("0", n_g1), rep("1", n_g2))
dt_base = data.table(Sample = 1:28, sand = group_var)

sim_list = list()
sim_list[[1]] = setDT(spider$abund)
sim_list[[1]][,null_model := FALSE]
for (i in 1:(2*n_sim)){
        
        ld_0 <-
                generate_copula_data(n_g1, marginal_details = y_0_mar, cov = y0pwa_corr_mcem_add)
        ld_1 <-
                generate_copula_data(n_g2, marginal_details = y_1_mar, cov = y1pwa_corr_mcem_add)
        
        mat_new = rbindlist(list(
                as.data.frame(ld_0[["observed"]]),
                as.data.frame(ld_1[["observed"]])
        ), fill = T)
        
        for (j in seq_len(ncol(mat_new)))
                set(mat_new, which(is.na(mat_new[[j]])), j, 0)
        
        mat_new = cbind(dt_base, mat_new)
        mat_new[, sand := paste0(sand, "_", i)]
        mat_new[, null_model := FALSE]
        if (i <= n_sim){
                sim_list[[i+1]] = mat_new
        } else{
                shuffle = sample(1:nrow(mat_new),nrow(mat_new))
                mat_new$sand = mat_new$sand[shuffle]
                mat_new[, null_model := TRUE]
                mat_new[, sand := paste0(sand, "_n", i-n_sim,"")]
                sim_list[[i+1]] = mat_new
        }
        
        rm(mat_new)
}

sim_mat = rbindlist(sim_list, use.names = TRUE, fill = TRUE)
sim_mat[, sand := factor(sand)]
saveRDS(sim_mat, "003_processed_data/copula/spider/y_super.RDS")
sim_mat = readRDS("003_processed_data/copula/anderson19/poor_knight_fish/y_super.RDS")

d_super = parallelDist(x = as.matrix(sim_mat[,-c(12:15)]), method = "bray", threads = 2)
saveRDS(d_super, "003_processed_data/copula/spider/d_super.RDS")
d_super = readRDS("003_processed_data/copula/anderson19/poor_knight_fish/d_super.RDS")

centroids = betadisper(d_super, sim_mat$sand, type = "centroid")
saveRDS(centroids, "003_processed_data/copula/spider/centroids.RDS")
centroids = readRDS("003_processed_data/copula/spider/centroids.RDS")

group_vec = rownames(centroids$centroids)
o_id     = which(!str_detect(group_vec, "_"))
null_id  = which(str_detect(group_vec, "n")) 
group_vec[null_id] = "null model"

group_vec %<>% str_remove("_.*")
plot_data = data.frame(axis1 = centroids$centroids[,1],
                       axis2 = centroids$centroids[,2],
                       axis3 = centroids$centroids[,3],
                       sand = factor(group_vec), 
                       og = FALSE)

plot_data$og[o_id] = TRUE
p12 = ggplot(data= filter(plot_data, sand != "null model"), aes(x = axis1, y = axis2, col = sand)) + 
        stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) +
        geom_point(aes(shape = sand)) + 
        #geom_point(data = filter(plot_data, og == TRUE), col = "black", size = 0.5) + 
        theme_classic()
p13 = ggplot(data= filter(plot_data, sand != "null model"), aes(x = axis1, y = axis3, col = sand)) + 
        stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
        geom_point(aes(shape = sand)) + 
        #geom_point(data = filter(plot_data, og == TRUE), col = "black", size = 0.5) + 
        theme_classic() +
        theme(legend.position = "none")
p23 = ggplot(data= filter(plot_data, sand != "null model"), aes(x = axis2, y = axis3, col = sand)) + 
        stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
        geom_point(aes(shape = sand)) + 
        #geom_point(data = filter(plot_data, og == TRUE), col = "black", size = 0.5) + 
        theme_classic()+
        theme(legend.position = "none")
p12_leg = cowplot::get_legend(p12)
cowplot::plot_grid(p12 + theme(legend.position = "none"),
                   p13,
                   p23#,
                   #p12_leg
)

p12 = ggplot(data= plot_data, aes(x = axis1, y = axis2, col = sand)) + 
        stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) +
        geom_point(aes(shape = sand)) + 
        #geom_point(data = filter(plot_data, og == TRUE), col = "black", size = 0.5) + 
        theme_classic()
p13 = ggplot(data= plot_data, aes(x = axis1, y = axis3, col = sand)) + 
        stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
        geom_point(aes(shape = sand)) + 
        #geom_point(data = filter(plot_data, og == TRUE), col = "black", size = 0.5) + 
        theme_classic() +
        theme(legend.position = "none")
p23 = ggplot(data= plot_data, aes(x = axis2, y = axis3, col = sand)) + 
        stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.5) + 
        geom_point(aes(shape = sand)) + 
        #geom_point(data = filter(plot_data, og == TRUE), col = "black", size = 0.5) + 
        theme_classic()+
        theme(legend.position = "none")
p12_leg = cowplot::get_legend(p12)
cowplot::plot_grid(p12 + theme(legend.position = "none"),
                   p13,
                   p23,
                   p12_leg)




colors3d = c('royalblue1', 'darkcyan', 'oldlace', "gray")
plot_data$color <- colors3d[ as.numeric(plot_data$season) ]
plot3d(x = plot_data$axis1, y = plot_data$axis2, z = plot_data$axis3, col = plot_data$color, type = 's', 
       radius = .005)