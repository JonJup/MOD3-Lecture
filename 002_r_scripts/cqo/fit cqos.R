### --- MOD3
### --- FIT CQOs
### --- HSPIDER data

library(VGAM)
library(purrr)
data(hspider)
hspider[, 1:6] <-
        scale(hspider[, 1:6]) # Standardized environmental variables

ls_cqo = list()
counter = 1
counter2 = 0
safe_cqo = safely(cqo)
for (loop_rank in 1:2) {
        for (loop_family in c("poissonff", "negbinomial")) {
                for (loop_eqtoler in c(TRUE, FALSE)) {
                        for (loop_itoler in c(TRUE, FALSE)) {
                                if (!loop_eqtoler & loop_itoler)
                                        next()
                                counter2 = 0
                                model_fit = FALSE
                                while (counter2 < 5 & !model_fit) {
                                        ls_cqo[[counter]] = safe_cqo(
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
                                                eq.toler = loop_eqtoler,
                                                I.tolerances = loop_itoler,
                                                trace = FALSE,
                                                Rank = loop_rank
                                        )
                                        counter2 = counter2 + 1
                                        if (is.null(ls_cqo[counter]$error))
                                                model_fit = TRUE
                                }
                                if (is.null(ls_cqo[[counter]]$error)){
                                        ls_cqo[[counter]] = ls_cqo[[counter]]$result
                                } else {
                                        next()
                                }
                                tolerance_name = ifelse(
                                        loop_itoler,
                                        "i",
                                        ifelse(loop_eqtoler, "e", "u")
                                )
                                
                                model_name = paste0(
                                        "r",
                                        loop_rank,
                                        "_",
                                        ifelse(loop_family == "poissonff","p","nb"),
                                        "_",
                                        tolerance_name
                                )
                                
                                names(ls_cqo)[counter] = model_name
                                counter = counter + 1
                                
                        }
                }
        }
        
}

saveRDS(object = ls_cqo, "003_processed_data/cqo/hspider/all_models.RDS")

