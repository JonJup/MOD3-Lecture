# Code S1.R
# R code function for estimating marginal distributions and their parameters
# using information criteria.


# chooseDistr () 
# Distributional choice for individual variables (marginals) in a multivariate count dataset.
#    Pois = Poisson,
#    NB = negative binmomial,
#    ZIP = zero-inflated Poisson,
#    ZINB = zero-inflated negative binomial.
#           (Note: zero-inflated models are fitted with just a single extra parameter:
#                  pi = the probability of an excess zero)
# The four models are fitted to each variable individually using ML, and the models are compared using information criteria.
#
# The function takes as input:
#    Y : (n x p) matrix of count data with n = no. samples and p = no. variables.
#    X : (n x q) option to input a matrix of factors or predictor variables of full rank
#                      (If not provided a simple intercept is used).
#    rel_cutoff: relative cutoff; e.g., default is 0.05: omit species that occur in less than 5% of samples
#    abs_cutoff: absolute cutoff; e.g., default is 1: omit species that occur in no more than 1 sample
#    criterion: which criterion should be used for ranking the models: "AIC", "AICc" or "BIC".
#               The default is "AICc".
#    imp_ord: whether or not to provide output of species in order of Importance, which is given
#             by their Frequency of occurrence. The default is FALSE.
#
# Output is a list with two items:
#    marginals: a data frame with the marginal distribution and the estimated parameters for each variable
#    diagnostics: a data frame summarising a host of diagnostics on each variable,
#                 including the chosen distribution and associated parameters.
#    detail: a list of length p containing a matrix of results for each species with results
#            obtained for each of the four models, including log-likelihoods and information criteria.
#
   chooseDistr = function(Y, X = 0, criterion = "AICc", rel_cutoff = 0.05, abs_cutoff = 1, imp_ord = FALSE) {
      require(pscl)
      require(MASS)
      N = dim(Y)[1]
      p = dim(Y)[2]
      if(X==0) {
         q=1
         } else {
         q = dim(model.matrix(~X))[2]
         }
      irank <- function (criterion) {
                 switch(criterion,
                    AIC = rank(AIC),
                    AICc = rank(AICc),
                    BIC = rank(BIC) )
               }
      detail = as.list(1:p); names(detail) <- colnames(Y)
   # Silence warnings, but then be sure to turn them back on at the end of the function!
      oldw <- getOption("warn")
      options(warn = -1)
   # Some initial diagnostics on the matrix of variables.
      Y.pa = Y>0
   # Calculate the Frequency of occurrence for each species and also the Percentage of samples in which they occur.
   # Importance is a ranking based on this Percentage.
      Freq = apply(Y.pa, MARGIN = 2, FUN = sum)
      Perc = Freq/N
      Importance = rank(100-Perc, ties.method = "average")
      Original.order = 1:p
      Max.Abund = apply(Y, MARGIN = 2, FUN = max)
      Tot.Abund = apply(Y, MARGIN = 2, FUN = sum)
      mu = apply(Y, MARGIN = 2, FUN = mean)
      SD = apply(Y, MARGIN = 2, FUN = sum)
      Var = apply(Y, MARGIN = 2, FUN = var)
      D = Var/mu # Dispersion index
      cut.off = max(c(abs_cutoff/N, rel_cutoff)) # allows for situation where 1/N is larger than 5%
      Rare = as.numeric(Perc <= cut.off)
      Singleton = as.numeric(Freq == 1)
      diagnostics = cbind(Original.order,Freq,Perc,Importance,Max.Abund,Tot.Abund,SD,Var,D,Rare,Singleton)
   # Consider each model for a given species
      Distr = rep(0,p)
      theta = rep(NA,p)
      pi = rep(NA,p)
      beta = matrix(rep(0,p*q),nrow=p,ncol=q)
      if(X==0) { colnames(beta) <- "(Intercept)"
         } else {
         colnames(beta) <- colnames(model.matrix(~X))
         }
      for (k in 1:p) {   # Begin p-species loop
      #   print(k)
         y = Y[,k]
      # Fit models and pull out relevant information from each...
         models = c("Pois","NB","ZIP","ZINB")
         nparms = c( q,q+1,q+1,q+2 )
         logl=rep(0,4)
      # The order of fit of the models is in increasing order of complexity
         if(X==0) {
         m1 = glm(y ~ 1, family = poisson(link = "log"))            #Poisson
         m2 = glm.nb(y ~ 1)                                         #NB
         if (min(y)==0) m3 = zeroinfl(y ~ 1 | 1, dist = "poisson", link = "logit") #ZIP
         if (min(y)==0) m4 = zeroinfl(y ~ 1 | 1, dist = "negbin", link = "logit")  #ZINB
         } else {
         m1 = glm(y ~ X, family = poisson(link = "log"))            #Poisson
         m2 = glm.nb(y ~ X)                                         #NB
         if (min(y)==0) m3 = zeroinfl(y ~ X | 1, dist = "poisson", link = "logit") #ZIP
         if (min(y)==0) m4 = zeroinfl(y ~ X | 1, dist = "negbin", link = "logit")  #ZINB
         }
      # Put these models in a list
         if (min(y)==0) {m = list(m1,m2,m3,m4)
                         names(m) <- models
                         n.m <- 4
                         } else {
                         m = list(m1,m2)
                         names(m) <- models[1:2]
                         n.m <- 2
                         logl <- c(0,0); nparms = c( q,q+1 )
                         }
         for (i in 1:n.m) { logl[i] = logLik(m[[i]]) }
         AIC = 2*nparms - 2*logl
         AICc = AIC + 2*nparms*(nparms+1)/(N-nparms-1)
         BIC = nparms*log(N) - 2*logl
         min.AIC = min(AIC)
         weights = exp((min.AIC - AIC)/2)
         results = cbind(nparms,logl,AIC,AICc,BIC,weights,irank(criterion))[1:n.m,]
         rownames(results) <- models[1:n.m]
         colnames(results) = c("nparms","logl","AIC","AICc","BIC","weight",paste("rank",criterion,sep="."))
         results.ord = results[order(irank(criterion)),]
         detail[[k]] = results.ord
      # And now the distribution of choice based on the chosen criterion should be named in the first position, hence
         Distr[k] = rownames(results.ord)[1] 
      # We can also pull out estimated values for parameters, including theta and pi, as may be appropriate
         if(Distr[k] == "Pois") { beta[k,]=m1$coef } 
         if(Distr[k] == "NB") { theta[k] = m2$theta ; beta[k,]=m2$coef }
         if(Distr[k] == "ZIP") { temp =  m3$coefficients$zero; pi[k] = 1/(1+exp(-temp)) ; beta[k,]=m3$coef$count}
         if(Distr[k] == "ZINB") { theta[k] = m4$theta; temp = m4$coefficients$zero
                                 pi[k] = 1/(1+exp(-temp)); beta[k,]=m4$coef$count }
      } # end of p-species loop
      marginals = data.frame(Distr, mu, theta, pi, beta)
      colnames(marginals) <- c("Distr","mu","theta","pi",colnames(beta))
      diagnostics = cbind(diagnostics, marginals)
      if (imp_ord == TRUE) { marginals = marginals[order(Importance), ]
                              diagnostics = diagnostics[order(Importance), ]
                              detail = detail[[order(Importance)]] }
      output = list(marginals,diagnostics,detail); names(output) <- c("marginals","diagnostics","detail")
      return(output)
   # Restore warnings to the original level.
      options(warn = oldw)
   } # end of function

# Usage:
#   my.output = chooseDistr(Y, X) # analysis with group structure
#   my.output = chooseDistr(Y) # analysis without group structure
#   my.output$marginals
#   my.output$diagnostics
#   my.output$detail
