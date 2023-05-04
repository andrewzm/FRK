library("DHARMa")

## Define a grid of BAUs
x = y = seq(0, 1, length.out = 100)
BAUs <- expand.grid(x = x, y = y)
coordinates(BAUs) = ~ x + y
gridded(BAUs) <- TRUE
N <- length(BAUs)
BAUs_df <- coordinates(BAUs) %>% as.data.frame 

## Construct some random effects: 5 vertical strips of the domain
set.seed(1)
ng  <- 5
fct <- cut(BAUs_df$x, ng)
sigma2gamma <- 2^2 # variance of random effects
gamma <- rnorm(ng, sd = sqrt(sigma2gamma))
names(gamma) <- levels(fct)

## Simulate the process, Y, over the BAUs 
sigma_xi <- 0.05
xi <- rnorm(N, sd = sigma_xi)
alpha <- c(3, 0.3)
BAUs_df$Y <- alpha[1] + alpha[2] * BAUs_df$x + xi + gamma[fct]
# BAUs_df$Y <- alpha[2] * BAUs_df$x + xi + gamma[fct]
# BAUs_df$Y <- alpha[2] * BAUs_df$x + xi + gamma[fct]
# BAUs_df$Y <- xi + gamma[fct] 

## Subsample some BAUs to act as observation locations, and simulate data
sigma_epsilon <- 0.05
n <- 1000
# obs_idx <- sample(1:length(BAUs), n) # random observations
obs_idx <- sample(which(BAUs_df$x < 0.29 | BAUs_df$x > 0.61), n) # missing strip missing so that one level is unobserved
zdf <- BAUs_df[obs_idx, ] %>% mutate(Z = rnorm(n, Y, sigma_epsilon))
coordinates(zdf) <- ~ x + y
plot_spatial_or_ST(zdf, column_names = "Z")

## Scalar fine-scale covariance matrix and add covariates
BAUs$fs   <- rep(1, N) 
BAUs$xcoordinate <- BAUs_df$x
BAUs$fct  <- fct

## Dummy Basis functions (tiny scale so that they are essentially ignored)
basis <- local_basis(
  manifold = plane(),
  as.matrix(expand.grid(c(0.25, 0.5, 0.75), c(0.25, 0.5, 0.75))),
  scale = rep(0.0001, 9), 
  regular = 1
)

## Construct object (various formula for testing different scenarios)
# M <- SRE(f = Z ~ 1 + xcoordinate, list(zdf), basis = basis, BAUs = BAUs, K_type = "precision")
M <- SRE(f = Z ~  1 + xcoordinate +  (1 | fct), list(zdf), basis = basis, BAUs = BAUs, K_type = "precision")
# M <- SRE(f = Z ~ 1 + xcoordinate, list(zdf), basis = basis, BAUs = BAUs, K_type = "precision")
# M <- SRE(f = Z ~ 0 + xcoordinate +  fct, list(zdf), basis = basis, BAUs = BAUs, K_type = "precision")
# M <- SRE(f = Z ~ 0 + xcoordinate +  (1 | fct), list(zdf), basis = basis, BAUs = BAUs, K_type = "precision")
# BAUs$dummy_x <- rnorm(N); M <- SRE(f = Z ~ -1 + dummy_x + (1 | fct), list(zdf), basis = basis, BAUs = BAUs, K_type = "precision")
M@Z # aggregated observations
M@X # fixed-effects design matrix
M@G # random-effects design matrix (a list)
M@S # basis-function design matrix

## These are different lengths if factor levels are missing in the data
M@G[[1]] %>% dim
M@G0[[1]] %>% dim 

## Fitting
opts_FRK$set("verbose", TRUE)
M <- SRE.fit(M, method = "TMB")

## Compare estimates to truth
M@mu_eta # close to zero as expected
alpha; M@alphahat  
gamma; M@mu_gamma  
sigma2gamma; M@sigma2gamma   
summary(xi); summary(as.numeric(M@mu_xi)) 
sqrt(sigma_epsilon^2 + sigma_xi^2); sqrt(M@sigma2fshat + M@Ve[1, 1]) 

## Prediction 
pred <- predict(M)

## Plotting
plotlist <- plot(M, pred$newdata)
ggpubr::ggarrange(plotlist = plotlist, nrow = 1, align = "hv", legend = "top")

## Coefficient estimates and uncertainty
coef_uncertainty(M)
coef_uncertainty(M, random_effects = TRUE)
coef(M)
coef(M, random_effects = TRUE) 
# NB above, we only return estimates of the observed random effects: the 
# unobserved random effects follow a mean-zero Gaussian distribution, so their 
# "estimate" is zero

## Model validation
DHARMa_object <- createDHARMa(simulate(M), zdf$Z)
plot(DHARMa_object)


# # Sanity check: Plot the G matrix versus the x coordinate 
# G <- M@G[[1]]
# X <- M@X
# head(G)
# head(X)
# plot(X[, "xcoordinate"], G[, 1], xlab = "x")
# plot(X[, "xcoordinate"], G[, 2], xlab = "x")
# plot(X[, "xcoordinate"], G[, 3], xlab = "x")
# plot(X[, "xcoordinate"], G[, 4], xlab = "x")
# plot(X[, "xcoordinate"], G[, 5], xlab = "x")