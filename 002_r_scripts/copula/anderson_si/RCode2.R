# Code S2.R
# R code function for calculating significant pair-wise associations
# using an adjusted family-wise error rate.

#
# pairWise ()
#
# Input a matrix of multivariate data and work out which pairs of variables
# are statistically significantly associated using permutation tests
# on the index of association under chosen criteria.
# Experimentwise error can be controlled using Family-wise error rates (FWER) and the
# distribution of minimum p-values, derived from the full set of permutation distributions
# across all tests (Wheldon et al. 2007).
# Note: eta = the index of association (IoA).
#
#
# The function takes as input:
#    Y : (N x p) matrix of count data with n = no. samples (rows) and p.full = no. variables (columns).
#    sig_level: alpha-level for per-comparison error rate (= 0.05 by default)
#    nperm: no. of permutations for each pairwise test (= 99999 by deafult)
#    alpha_type: type of pairwise comparison to use (= "FWER by deafult);
#                this can be one of: ("PCER","FDR","FWER","BF"), where:
#                   PCER = per-comparison error rate (alpha[1] = sig.level)
#                   FDR = false discovery rate
#                   FWER = exact family-wise error-rate
#                   BF = Bonferroni error rate (alpha[4] = sig.level/no.tests) 
#    graphic: logical to choose whether or not to output the graphic (default = TRUE)
#    rel_cutoff: relative cutoff; e.g., default is 0.05: omit species that occur in less than 5% of samples
#    abs_cutoff: absolute cutoff; e.g., default is 1: omit species that occur in no more than 1 sample
#
# Output is a list with the following items: 
#        alpha: adjusted significance level (alpha.prime) used for each type of multiple-test procedure
#        associated: names of associated species
#        IoA.obs: (p x p) full matrix of observed values for the index of association (IoA)
#        sig.pairs: matrix of (1,0) values identifying pairs of species having sig. association.
#        IoA.subset: matrix of IoA values for the subset of species that showed any sig. association.
#        IoA.shrunk: matrix of IoA values for this subset with only the sig. IoA values shown (zeros elsewhere).
#        results.df: detailed results matrix (e.g., for graphics, etc.)
#
   pairWise = function(Y, sig_level = 0.05, nperm = 99999, alpha_type = "FWER",
                       graphic = TRUE, rel_cutoff = 0.05, abs_cutoff = 1) {

      N = dim(Y)[1] # total number of sample units 
      p = dim(Y)[2] # total number of species

# Calculate the Frequency of occurrence for each species and also the Percentage of samples in which they occur.
# Importance is a ranking based on this Percentage.
      cut.off = max(c(abs_cutoff/N, rel_cutoff)) # allows for situation where 1/N is larger than 5%
      Y.pa = Y>0
      Freq = apply(Y.pa, MARGIN = 2, FUN = sum)
      Perc = Freq/N
      Importance = rank(100-Perc, ties.method = "average")
      Original.order = 1:p
      Rare = as.numeric(Perc <= cut.off)

# Replace Y with Y.important, which has omitted the rare species and is given
# in decreasing order of importance (based on freq of occurrence).
# also replace p with p.common
      Y.common = Y[,Rare!=1.0] # omit rare species
      p.common = dim(Y.common) [2]
      Y.important = Y.common[,order(Importance[Rare!=1.0])]
      Y = Y.important; p = p.common

# Calculate eta, the measure of association between every pair of variables.
# The Index of Association
   eta.pair = function(Y1,Y2) { 1 - (sum(abs(Y1/sum(Y1)-Y2/sum(Y2)))/sum(Y1/sum(Y1)+Y2/sum(Y2))) }
# Get the full matrix of indices (all pairs)
   eta = function (Y) {
      N=dim(Y)[1] # total number of sample units 
      p=dim(Y)[2]     
      d=matrix(0,nrow=p,ncol=p) # Note: we will be ignoring the diagonal, so it doesn't matter that it is zero here instead of the desired value of 1.0
      for(i in 1:(p-1)) {
         for(j in (i+1):p) {
            d[i,j] = eta.pair(Y[,i],Y[,j])
         }
      }
      d = d+t(d)
      return(as.dist(d)) # treat it as a distance matrix object, to easily pull out just the off-diagonal elements as a vector.
      }
   IoA = eta(Y)
   IoA.obs=as.matrix(IoA)
   colnames(IoA.obs) <- rownames(IoA.obs) <- colnames(Y)

# Now calculate the permutation distribution for each pair
   nset = p*(p-1)/2 # number of pairs = number of individual tests to be done.
#   length(as.vector(IoA)); nset # check that this is what we expect it to be, i.e., equal to nset.
   nperm = nperm # number of permutations
# We will need the following:
# Iperm = the distribution of the test statistic for each dataset under permutation.
# Ival = the value of the test statistic for each dataset.
# pval = the p-value obtained for each dataset.
# psort = p-values sorted in ascending order (to obtain FDR error rate)
# pstar = the p-values that would have been obtained for each Iperm, if it had been the observed one
# pmin = the minimum p-value obtained across all datasets (of length nperm)
# Icrit = the critical values of F associated with each of the percentiles of the distribution.

   Ival <- pval <- psort <- rep(0,nset)
   name1 <- name2 <- rep(0,nset)
   Iperm <- pstar <- matrix(rep(0,(nperm+1)*nset),nrow=(nperm+1),ncol=nset)
   alpha <- rep(0,4); crit.hi <- crit.low <- list(rep(0,4))        # Note: "crit" is a function
   Icrit.hi <- Icrit.low <- matrix(rep(0,nset*4),nrow=nset,ncol=4) # Critical values are held in Icrit

   iset = 0
   for (i in 1:(p-1)) {
      y1 = Y[,i]
      for (j in (i+1):p) {
         iset = iset+1
         y2 = Y[,j]
         Ival[iset] = eta.pair(y1,y2)
         name1[iset] = colnames(Y)[i]
         name2[iset] = colnames(Y)[j]
         Iperm[1,iset] = Ival[iset] # Note that Iperm includes the observed value
         for (iperm in 2:(nperm+1)) {
            Iperm[iperm,iset] = eta.pair(y1,y2[sample(N)])
         } # end permutation loop
      } # end loop over second variable
   } # end loop over first variable

# A function for calculating the empirical cdf (i.e. individual p-values for each member of a null empirical distribution)
# Note: p.up = prob of getting a value >= observed.
#       p.low = prob of getting a value <= observed.   # Note: p.up should also be: (rank(x)-1)/length(x)
   px = function(x) {p.low <- rank(x)/length(x); p.up <- (1.0 - p.low);
                     mat <- rbind(p.low,p.up)
#                     2*apply(mat,MARGIN = 2, FUN = min) # This would make the analysis a two-tailed test (p-values are uniform on [0,1.0]   
                     apply(mat,MARGIN = 2, FUN = min) # We can look at the probabilities in either tail
                                                      # The dist of these p-values under a true null is uniform on [0,0.5]
                                                      # then the use of alpha/2 and 1 - alpha/2 to choose critical values should match these p-values.
                     }
   pstar = apply(Iperm, MARGIN = 2, FUN = px)
   pval = pstar[1,] # raw p-values for each dataset.
# Now get the minimum p-values across all datasets for each permutation:
   pmin = apply(pstar, MARGIN = 1, FUN = min)

# OK, Let's get the corrected alpha levels for each of the four methods:
# alpha[1] = Per-comparison error rate (PCER) (e.g., = 0.05)
# alpha[2] = False discovery rate (FDR)
# alpha[3] = Family-wise error rate (FWER)
# alpha[4] = Bonferroni error rate (BF)
   names(alpha) <- c("PCER","FDR","FWER","BF")
   alpha[1] = sig_level
   alpha[4] = sig_level/nset
# Get the FWER
# To figure out what adjusted significance level should be used to maintain
# family-wise error rate at alpha (e.g., = 0.05), then we need to find the
# alpha-quantile (e.g. the 5th percentile) of the pmin distribution.
   alpha[3] = quantile(pmin, probs = sig_level)
# Finally, for the false discovery rate...
# Sort the original p-values from smallest to largest...
   psort = sort(pval, decreasing = FALSE)
# Now, starting with the largest p-value and working backwards,
# determine the largest i for which psort[i] <= (i/nset)*sig_level
   index=1:length(psort)
   ind = psort<= (index/nset)*sig_level  
   k = max(index[ind])
   alpha[2]=psort[k]

# If you never got a situation where p-values were small enough then use Bonferroni.
   if(as.numeric(sum(ind)==0) | (psort[k] < alpha[4]))  { alpha[2]=alpha[4] } 

# Calculate critical values for each Iperm distribution.
# i.e., get the quantile corresponding to each alpha-level for each permutation distribution
   for (i in 1:4) {
      crit.low[[i]] <- function (x) {quantile(x,prob = (alpha[i]/2) )}
      crit.hi[[i]] <- function (x) {quantile(x,prob = (1 - alpha[i]/2) )}
      Icrit.low[,i] <- apply(Iperm, MARGIN=2, FUN = crit.low[[i]])
      Icrit.hi[,i] <- apply(Iperm, MARGIN=2, FUN = crit.hi[[i]])
      }
#
# Identify "significance", using the chosen criterion   
   itype = switch(alpha_type, PCER=1, FDR=2, FWER=3, BF=4)
   sig.type = paste("sig.",alpha_type, sep="")
   IcritL = Icrit.low[,itype]
   IcritH = Icrit.hi[,itype]

   sig = pval <= alpha[itype]/2

# Note: when we use IoA, there will be an IoA value of zero between two species that never occurred
# together anywhere at all. We may decide that we should not try to ascertain the relationship
# between any two species that never occur in any sites together at all. This is directly analogous to the
# notion that we can never ascertain the similarity between two sites that have no species in common - they are
# just always going to be 100% dissimilar.
# 
# Information in the detailed results output are:
# Ival, pval, and Icrit for each of the 4 methods in turn for each pair-wise test

   results.df=data.frame(name1,name2,Ival,pval,sig,Icrit.low,Icrit.hi)
   colnames(results.df) <- c("name1","name2","IoA","p.value",sig.type,
                          "IcritL.PCER","IcritL.FDR","IcritL.FWER","IcritL.BF",
                          "IcritH.PCER","IcritH.FDR","IcritH.FWER","IcritH.BF" )

# Draw the results...
   if(graphic == TRUE) {
   x=1:length(pval)

   plot(x,results.df$IoA, type = "n", xlab="Pair of variables", ylab = "Index of Association",ylim=c(0,1.0),las=1)
   points(x,results.df$IoA,col = "black", pch = 19, cex=0.5)
   points(x[sig],results.df$IoA[sig],col = "red", cex = 1.5)  
   lines(x,Icrit.low[,itype], col = "blue")
   lines(x,Icrit.hi[,itype], col = "blue")
   }

# "Rewind" the matrix into a correlation format...
   sig.mat = diag(rep(NA,p))
   IoA2 = diag(rep(0.5,p))
   iset=0
   for (i in 1:(p-1)) {
      for (j in (i+1):p) {
         iset = iset+1
         IoA2[i,j] = Ival[iset] 
         sig.mat[i,j] = sig[iset]
      } # end loop over second variable
   } # end loop over first variable
   IoA2 = IoA2 + t(IoA2)
   sig.mat = sig.mat + t(sig.mat)
   if(sum(sig) != 0 ) {
# Subset of variables for which modeling of correlation structures is worthwhile:
   pull = apply(sig.mat, MARGIN = 2, FUN = sum, na.rm = TRUE)
# species numbers you want to retain.
   sub = (1:p)[pull!=0]
   sig.pairs = sig.mat[sub,sub];
   rownames(sig.pairs) <- colnames(sig.pairs) <- colnames(Y)[sub] 
   IoA.subset = IoA2[sub,sub]
   rownames(IoA.subset) <- colnames(IoA.subset) <- colnames(Y)[sub] 
   IoA.shrunk = IoA.subset * sig.pairs
   for (i in 1:length(sub)) {IoA.shrunk[i,i] = 1.0}
   associated <-colnames(Y)[sub]
   }
   else {
   associated = "none"
   }

# Kick out the output...
   output = list(alpha,associated,IoA.obs,sig.pairs,IoA.subset,IoA.shrunk,results.df)
   names(output) = c("alpha.adj","associated","IoA.obs","sig.pairs","IoA.subset","IoA.shrunk","results.df")
   return(output)

   } # end

# Usage
#   my.pairs = pairWise(Y, sig_level = 0.05, nperm = 99999, alpha_type = "FWER",
#                       graphic = TRUE, rel_cutoff = 0.05, abs_cutoff = 1)
#   my.pairs$results.df

 
