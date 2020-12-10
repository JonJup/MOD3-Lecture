### --- MOD3
### --- FIT caos
### --- HSPIDER data

library(VGAM)
library(purrr)
data(hspider)
hspider[, 1:6] <-
        scale(hspider[, 1:6]) # Standardized environmental variables

ls_cao = list()
counter = 1
counter2 = 0
safe_cao = safely(cao)

for (loop_family in c("poissonff", "negbinomial")) {
        counter2 = 0
        model_fit = FALSE
        while (counter2 < 5 & !model_fit) {
                ls_cao[[counter]] = safe_cao(
                        cbind(
                                Alopacce,
                                Alopcune,
                                Alopfabr,
                                Arctlute,
                                Arctperi,
                                Auloalbi,
                                Pardlugu,
                                Pardmont,
                                Pardnigr,
                                Pardpull,
                                Trocterr,
                                Zoraspin
                        ) ~
                                WaterCon + BareSand + FallTwig +
                                CoveMoss + CoveHerb + ReflLux,
                        family = loop_family,
                        data = hspider,
                        trace = FALSE
                )
                counter2 = counter2 + 1
                if (is.null(ls_cao[counter]$error))
                        model_fit = TRUE
        }
        if (is.null(ls_cao[[counter]]$error)) {
                ls_cao[[counter]] = ls_cao[[counter]]$result
        } else {
                next()
        }
        model_name = paste0(ifelse(loop_family == "poissonff", "p", "nb"))
        
        names(ls_cao)[counter] = model_name
        counter = counter + 1
        
}



saveRDS(object = ls_cao, "003_processed_data/cao/hspider/all_models.RDS")

