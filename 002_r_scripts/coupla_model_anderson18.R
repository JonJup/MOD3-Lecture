### Copula example ### 



pacman::p_load(here, readxl, data.table,dplyr)
setwd(here())
# load data  --------------------------------------------------------------

## -- fish data -- ##

dt_fish <- read_excel(path = "../../../Supplementary Informations/Anderson19Copula/ece34948-sup-0001-tables1.xlsx",
                      skip = 1)

## -- function -- ## 
# load choose distr
source("../../../Supplementary Informations/Anderson19Copula/RCode1.R")


# carpet data -------------------------------------------------------------
# names of seasons 
dt_fish %>% select(Time) %>% unique

tbl_to_mat <- function(x, sub){
        out <- filter(dt_fish, Time == sub) %>% select(-c(Sample, Time)) %>% as.matrix()
}
dt_fish_s98 <- tbl_to_mat(sub = "Sep.98")
dt_fish_m99 <- tbl_to_mat(sub = "Mar.99") 
dt_fish_s99 <- tbl_to_mat(sub = "Sep.99") 

# check if all rows are included in subsets 
nrow(dt_fish_m99) + nrow(dt_fish_s98) + nrow(dt_fish_s99) == nrow(dt_fish)

chooseDistr(Y = dt_fish_s98)
chooseDistr(Y = dt_fish_m99[,3:ncol(dt_fish_m99)])
chooseDistr(Y = dt_fish_s99[,3:ncol(dt_fish_s99)])

chooseDistr()