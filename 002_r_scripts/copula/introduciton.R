
library(copula)
library(ggplot2)
library(GoFKernel)
library(magrittr)
n = 500
d = 2

normal <- normalCopula(param = 0.8, dim = d)
x_cop_mvn = mvdc(copula = normal,margins = c("norm", "norm"),paramMargins = list(list(mean = 0, sd=1), list(mean = 0, sd=1)))
x_cop_exp = mvdc(copula = normal,margins = c("exp", "exp"), paramMargins = list(list(rate = 1), list(rate = 1)))

n_samples <- rMvdc(n, x_cop_mvn) %>% data.frame()
x_samples <- rMvdc(n, x_cop_exp) %>% data.frame()

names(n_samples) = c("x1", "x2")
names(x_samples) = c("y1", "y2")

n_dens = n_samples %>% tidyr::pivot_longer(cols = c("x1", "x2"))
x_dens = x_samples %>% tidyr::pivot_longer(cols = c("y1", "y2"))

ggplot(data = n_samples, aes(x=x1,y=x2)) + geom_point() + ggtitle(paste("Correlation:", round(cor(n_samples[,1], n_samples[,2]),2)))
ggplot(data = x_samples, aes(x=y1,y=y2)) + geom_point() + ggtitle(paste("Correlation:", round(cor(x_samples[,1], x_samples[,2]),2)))

ggplot(data = n_dens, aes(x = value, fill = name, col = name)) + geom_density(alpha = 0.5)
ggplot(data = x_dens, aes(x = value, fill = name, col = name)) + geom_density(alpha = 0.5)

mvn_cdf = ecdf(as.matrix(n_samples))
exp_cdf = ecdf(as.matrix(x_samples))

n_pit = matrix(0, nrow = n, ncol = d)
x_pit = matrix(0, nrow = n, ncol = d)

n_pit[,1] = mvn_cdf(n_samples[,1])
n_pit[,2] = mvn_cdf(n_samples[,2])
x_pit[,2] = exp_cdf(x_samples[,2])
x_pit[,1] = exp_cdf(x_samples[,1])

ggplot(data = data.frame(n_pit), aes(x=X1,y=X2)) + geom_point() + ggtitle(paste("Correlation:", round(cor(n_pit[,1], n_pit[,2]),2)))
ggplot(data = data.frame(x_pit), aes(x=X1,y=X2)) + geom_point() + ggtitle(paste("Correlation:", round(cor(x_pit[,1], x_pit[,2]),2)))

## -- Quantile transfrom -- ## 

#mvn_icdf = inverse(f = function(x) pnorm(x, mean = 0, sd = 1))
n_qit    = qnorm(p = n_pit) %>% data.frame 
x_qit    = qnorm(p = x_pit) %>% data.frame 

n_qit <- n_qit[!is.infinite(rowSums(n_qit)),]
x_qit <- x_qit[!is.infinite(rowSums(x_qit)),]

names(n_qit) = c("x1", "x2")
names(x_qit) = c("y1", "y2")

ggplot(data = n_qit, aes(x=x1,y=x2)) + geom_point() + ggtitle(paste("Correlation:", round(cor(n_qit[,1], n_qit[,2]),2)))
ggplot(data = x_qit, aes(x=y1,y=y2)) + geom_point() + ggtitle(paste("Correlation:", round(cor(x_qit[,1], x_qit[,2]),2)))

both_dens = cbind(n_qit, x_qit)
both_dens %<>% tidyr::pivot_longer(cols = c("x1", "x2", "y1", "y2"))

ggplot(data = both_dens, aes(x = value, fill = name, col = name)) + geom_density(alpha = 0.5)



# Copula package
library(copula)
# Fancy 3D plain scatterplots
library(scatterplot3d)
# ggplot2
library(ggplot2)
# Useful package to set ggplot plots one next to the other
library(grid)
set.seed(235)

# Generate a bivariate normal copula with rho = 0.7
normal <- normalCopula(param = 0.7, dim = 2)
# Generate a bivariate t-copula with rho = 0.8 and df = 2
stc <- tCopula(param = 0.8, dim = 2, df = 2)

frank <- frankCopula(dim = 2, param = 8)
gumbel <- gumbelCopula(dim = 3, param = 5.6)
clayton <- claytonCopula(dim = 4, param = 19)

# Print information on the Frank copula
print(frank)

cp <- claytonCopula(param = c(3.4), dim = 2)

# Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
multivariate_dist <- mvdc(copula = cp,
                          margins = c("norm", "t"),
                          paramMargins = list(list(mean = 2, sd=3),
                                              list(df = 2)) )

print(multivariate_dist)


fr <- rCopula(2000, frank)
gu <- rCopula(2000, gumbel)
cl <- rCopula(2000, clayton)

p1 <- qplot(fr[,1], fr[,2], colour = fr[,1], main="Frank copula random samples theta = 8", xlab = "u", ylab = "v")
p2 <- qplot(gu[,1], gu[,2], colour = gu[,1], main="Gumbel copula random samples theta = 5.6", xlab = "u", ylab = "v") 
p3 <- qplot(cl[,1], cl[,2], colour = cl[,1], main="Clayton copula random samples theta = 19", xlab = "u", ylab = "v")
pushViewport(viewport(layout = grid.layout(1, 3)))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))

samples <- rMvdc(2000, multivariate_dist)
scatterplot3d(samples[,1], samples[,2], color = "blue",pch = ".")

# Generate the normal copula and sample some observations
coef_ <- 0.8
mycopula <- normalCopula(coef_, dim = 2)
u <- rCopula(2000, mycopula)

# Compute the density
pdf_ <- dCopula(u, mycopula)

# Compute the CDF
cdf <- pCopula(u, mycopula)

# Generate random sample observations from the multivariate distribution
v <- rMvdc(2000, multivariate_dist)

# Compute the density
pdf_mvd <- dMvdc(v, multivariate_dist)

# Compute the CDF
cdf_mvd <- pMvdc(v, multivariate_dist)


par(mfrow = c(1, 3))
# 3D plain scatterplot of the density, plot of the density and contour plot
scatterplot3d(u[,1], u[,2], pdf_, color="red", main="Density", xlab ="u1", ylab="u2", zlab="dCopula", pch=".")
persp(mycopula, dCopula, main ="Density")
contour(mycopula, dCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")

par(mfrow = c(1, 3))
# 3D plain scatterplot of the CDF, plot of the CDF and contour plot
scatterplot3d(u[,1], u[,2], cdf, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pCopula",pch=".")
persp(mycopula, pCopula, main = "CDF")
contour(mycopula, pCopula, xlim = c(0, 1), ylim=c(0, 1), main = "Contour plot")

# 3D plain scatterplot of the multivariate distribution
par(mfrow = c(1, 2))
scatterplot3d(v[,1],v[,2], pdf_mvd, color="red", main="Density", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")
persp(multivariate_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Density")
contour(multivariate_dist, dMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")
persp(multivariate_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "CDF")
contour(multivariate_dist, pMvdc, xlim = c(-4, 4), ylim=c(0, 2), main = "Contour plot")