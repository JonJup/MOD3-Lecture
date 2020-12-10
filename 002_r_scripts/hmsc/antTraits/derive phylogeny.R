# .. HMSC on AntTrits data.. # 

pacman::p_load(ape, Hmsc, data.table, magrittr, dplyr, gllvm, stringr, rotl)

data(antTraits)

# derive phylogeny for ants 
taxa <- colnames(antTraits$abund)
taxa %<>% str_replace_all("\\.", "\\ ") %>%
        str_replace_all("sp\\ .*", "") %>%
        str_trim()
names(antTraits$abund) = taxa
taxa = taxa[str_detect(taxa, " ")]
antTraits$abund = antTraits$abund[, names(antTraits$abund) %in% taxa]
resolved_names  <- rotl::tnrs_match_names(names = taxa)

names(antTraits$abund)[which(names(antTraits$abund) == "Cardiocondyla nuda atalanta")] = "Cardiocondyla nuda"
names(antTraits$abund)[which(names(antTraits$abund) == "Myrmecia pilosula complex")] = "Myrmecia pilosula"

resolved_names  <- rotl::tnrs_match_names(names = names(antTraits$abund))
any(duplicated(resolved_names$unique_name))

antTraits$abund$`Iridomyrmex suchieri` = antTraits$abund$`Iridomyrmex suchieri` + antTraits$abund$`Iridomyrmex suchieroides`
antTraits$abund = antTraits$abund[, !names(antTraits$abund) %in% "Iridomyrmex suchieroides"]

resolved_names  <- rotl::tnrs_match_names(names = names(antTraits$abund))
any(duplicated(resolved_names$unique_name))

ott_ids <- ott_id(resolved_names)

species_tree <- rotl::tol_induced_subtree(ott_ids, 
                            label_format = "name")

species_tree$edge
rm(resolved_names, taxa, ott_ids)
        
# adjust traits to new data 
dt_traits = copy(antTraits$traits)
dt_traits %<>% setDT
dt_traits[, taxon := rownames(antTraits$traits) %>% str_replace_all("\\.", "\\ ") %>%
                  str_replace_all("sp\\ .*", "") %>%
                  str_trim()]
dt_traits = dt_traits[str_detect(taxon, " ")]
dt_traits[taxon == "Cardiocondyla nuda atalanta", taxon := "Cardiocondyla nuda"]
dt_traits[taxon == "Myrmecia pilosula complex", taxon := "Myrmecia pilosula"]
dt_traits = dt_traits[taxon != "Iridomyrmex suchieroides"]

df_traits = data.frame(dt_traits)
rownames(df_traits) = dt_traits$taxon
df_traits %<>% select(-taxon)
names(df_traits)[2] = "spines"

df_studydesign = data.frame(site = paste(1:30))
rL = HmscRandomLevel(N = 30)
rL$nfMin = 5
rL$nfMax = 10
# define hmsc model  ------------------------------------------------------

m_full = Hmsc(
        Y = as.matrix(antTraits$abund),
        XData = antTraits$env,
        XFormula =  ~ .,
        TrData = df_traits,
        TrFormula = ~ Femur.length + Pilosity + spines + Webers.length,
        #phyloTree = species_tree,
        distr = "lognormal poisson", 
        studyDesign = df_studydesign,
        ranLevels = list(site = rL)
        )

nChains = 2
thin = 1000
nParallel = max(round(parallel::detectCores() / 2), nChains)
samples = 1000
transient = 500 * thin
verbose = 500 * thin

m_full_fit = sampleMcmc(
        m_full,
        thin = thin,
        samples = samples,
        transient = transient,
        nChains = nChains,
        verbose = verbose,
        initPar = "fixed effects"
)

saveRDS(m_full_fit, "003_processed_data/hmsc/spider/fit_thin1000.RDS")

mpost = convertToCodaObject(m_full_fit,
                            spNamesNumbers = c(T, F),
                            covNamesNumbers = c(T, F))
plot(mpost$Beta)


# check convergence statistics 
es.beta  = effectiveSize(mpost$Beta)
ge.beta  = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
es.gamma = effectiveSize(mpost$Gamma)
ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
es.V = effectiveSize(mpost$V)
ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
es.omega = effectiveSize(mpost$Omega[[1]])
ge.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
ge = list(
        beta = ge.beta,
        gamma = ge.gamma,
        V      = ge.V,
        Omega = ge.omega
        
)
es =  list(
        beta = es.beta,
        gamma = es.gamma,
        V = es.V,
        Omega = es.omega
)
par(mfrow=c(2,2))
for (i in seq_along(es))vioplot::vioplot(es[[i]], main = names(es)[i])

for (i in seq_along(ge)) vioplot::vioplot(ge[[i]][,1], main = names(ge)[i])

