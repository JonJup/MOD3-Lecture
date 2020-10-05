### Copula example ### 



pacman::p_load(corrplot, here, readxl, data.table,dplyr, stringr)
setwd(here())
# load data  --------------------------------------------------------------

## -- fish data -- ##

dt_fish <- read_excel(path = "../../../Supplementary Informations/Anderson19Copula/ece34948-sup-0001-tables1.xlsx",
                      skip = 1)

## -- function -- ## Here I load the function provided in the supplementary
#materials of Anderson et al. 2019. The original third code file also contains
#an example use case. I split the file into one containing the functions which
#we load here and another one which contains the example use case. The
#delineation is made clear in the RCode3.R file.
source("../../../Supplementary Informations/Anderson19Copula/RCode1.R")
source("../../../Supplementary Informations/Anderson19Copula/RCode2.R")



# carpet data ------------------------------------------------------------- 
# What unique seasons are in the data? The select verb from dplyr is explicitly
# loaded from dplyr through the :: sign because the MASS package which will be
# loaded inside the chooseDistr() function also contains a select() function
# which will mask dplyr's.
vec_unique_seasons <-
        dt_fish %>%
        dplyr::select(Time) %>%
        pull() %>%
        unique

# Define a function to: 1. Subset data to a certain value of the variable
# "Time", 2. remove the two Variables Sample and Time, and 3. convert the result
# to a matrix. The last step is done because the chooseDistr() function requires
# a matrix as input. This is done in a function to adhere to the DRY principle:
# Don't repeat yourself. Which mean that in coding you should avoid writing the
# same lines of code over and over again to apply them to different objects.
# This would make it very cumbersome to fix errors and also to rerun the
# analysis.

tbl_to_mat <- function(sub){
        out <- filter(dt_fish, Time == sub) %>%
                dplyr::select(-c(Sample, Time)) %>%
                as.matrix()
}
#The original season names have dots (.) in them. In coding in general you
#should avoid dots in the names of objects. In R this is not really a problem
#but in other languages they can cause errors, so its a good habit to use
#underscores instead.
for (i in vec_unique_seasons){
        cr_new_i <- i %>%
                str_replace(pattern = "\\.", replacement = "_") %>%
                str_to_lower()
        assign(x = paste0("mt_",cr_new_i),
               value = tbl_to_mat(sub=i))
        rm(i, cr_new_i)
}

rm(tbl_to_mat, dt_fish, vec_unique_seasons)


# In this loop we first remove species that do not occur in the subsets (i.e.
# where the sum of the column is zero). They would cause errors in the
# chooseDistr() function. Afterwards we call the function on each matrix.
for(i in ls()[grepl(pattern="mt", x = ls())]){

        mt_loop                         <- get(i)
        rm_id                           <- which(colSums(mt_loop) == 0)
        if (length(rm_id) != 0) mt_loop <- mt_loop[,-rm_id]
        l_cd_out                        <- chooseDistr(Y = mt_loop)
        df_cd_out                       <- l_cd_out$marginals
        cr_save_name_part               <- paste(str_split_fixed(string=i,pattern="_",n=3) %>% .[,2:3], collapse = "_")
        cr_save_name                    <- paste0("df_", cr_save_name_part)
        assign(x=i,           value=mt_loop)
        assign(x=cr_save_name,value=df_cd_out)
        rm(i, mt_loop,l_cd_out,df_cd_out,cr_save_name_part,cr_save_name,rm_id);gc()
}


# Now we now the optimal marginal distributions for each species and can turn to
# the significant pair-wise associations. The function is not well written and
# occaisonally returns an error. To avoid having to run in manually over and
# over again I wrote a while loop.

done=FALSE
safeapply <- purrr::safely(lapply)
counter=1

while(done == FALSE) {
        li_pw_all_out = safeapply(
                list(mt_sep_98,
                     mt_mar_99,
                     mt_sep_99),
                pairWise,
                nperm = 999,
                graphic = FALSE
        )
        if (length(li_pw_all_out$error) == 0)
                done = TRUE
        
        print(counter)
        
        counter <- counter + 1
}
rm(counter, done, safeapply)
# extract results from safely list 
li_pw_all_out <- li_pw_all_out$result

# get single seasons 
pwa_sep_98 <- li_pw_all_out[[1]]
pwa_mar_99 <- li_pw_all_out[[2]]
pwa_sep_99 <- li_pw_all_out[[3]]

# We can construct correlation plots just like in the publication
# For Index of the observed association (IoA): 
pwa_sep_98$IoA.obs %>% 
        corrplot(diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
# ... and the shrunken IoA:
pwa_sep_98$IoA.shrunk %>% 
        corrplot(diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")
# ... and the subset of statistically significant IoAs:
pwa_sep_98$IoA.subset %>% 
        corrplot(diag = FALSE, type = "lower", tl.cex = 0.8, tl.srt = 45, tl.col = "black")

# add species column. Later functions require this column.
df_sep_98$Species <- rownames(df_sep_98)
df_mar_99$Species <- rownames(df_mar_99)
df_sep_99$Species <- rownames(df_sep_99)

# ids of significant pairwise associated taxa
s98mt_id <- which(colnames(mt_sep_98) %in%  pwa_sep_98$associated)
s98df_id <- which(rownames(df_sep_98) %in%  pwa_sep_98$associated)
m98mt_id <- which(colnames(mt_mar_99) %in%  pwa_mar_99$associated)
m98df_id <- which(rownames(df_mar_99) %in%  pwa_mar_99$associated)
s99mt_id <- which(colnames(mt_sep_99) %in%  pwa_sep_99$associated)
s99df_id <- which(rownames(df_sep_99) %in%  pwa_sep_99$associated)

# new data sets that only hold these subsets 
mt_sep_98_assoc <- mt_sep_98[,s98mt_id]
mt_mar_99_assoc <- mt_mar_99[,s99mt_id]
mt_sep_99_assoc <- mt_sep_99[,s99mt_id]
df_sep_98_assoc <- df_sep_98[s98df_id,]
df_mar_99_assoc <- df_mar_99[s99df_id,]
df_sep_99_assoc <- df_sep_99[s99df_id,]





source("../../../Supplementary Informations/Anderson19Copula/RCode3_functions.R")
li_copula_sep_98 <- estimate_copula(data = mt_sep_98_assoc,marginal_details = df_sep_98_assoc)
li_copula_mar_99 <- estimate_copula(data = mt_mar_99_assoc,marginal_details = df_mar_99_assoc)
li_copula_sep_99 <- estimate_copula(data = mt_sep_99_assoc,marginal_details = df_sep_99_assoc)


# Display a summary of the MCEM result
summary(li_copula_sep_98)

# Extract the MLE correlation matrix from the mcem_result list.
corr_mcem_sep_98 <- li_copula_sep_98[["cov_final"]]
corr_mcem_mar_99 <- li_copula_mar_99[["cov_final"]]
corr_mcem_sep_99 <- li_copula_sep_99[["cov_final"]]

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

corr_mcem_sep_98_add <- add_unass(data = corr_mcem_sep_98, spe_names = rownames(df_sep_98))
corr_mcem_mar_99_add <- add_unass(data = corr_mcem_mar_99, spe_names = rownames(df_mar_99))
corr_mcem_sep_99_add <- add_unass(data = corr_mcem_sep_99, spe_names = rownames(df_sep_99))

N = 1000

# Use the result of the MCEM calculation to generate a new set of data.
simulated_data <- generate_copula_data(N, marginal_details = df_sep_98, cov = corr_mcem_sep_98_add)
simulated_data <- generate_copula_data(N, marginal_details = df_mar_99, cov = corr_mcem_mar_99_add)
simulated_data <- generate_copula_data(N, marginal_details = df_sep_99, cov = corr_mcem_sep_99_add)

library(vegan)
test <- adonis(simulated_data)


# Produce some plots of the original set of data and the data generated using the correlation
# matrix found by the MCEM algorithm.
pch <- c(rep(0, N), rep(20, N))
species_indices <- c(1,2)
plot(rbind(full_data[["observed"]][ , species_indices], simulated_data[["observed"]][ , species_indices]), pch = pch)
legend(x = "topright", pch = c(0, 20), legend = c("observed", "simulated"))

species_indices <- c(3,4)
plot(rbind(full_data[["observed"]][ , species_indices], simulated_data[["observed"]][ , species_indices]), pch = pch)
legend(x = "topright", pch = c(0, 20), legend = c("observed", "simulated"))


####### End of EXAMPLE Usage #########



### ---- old shizzle ---- ##### 
