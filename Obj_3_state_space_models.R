

# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #
# 3. State-Space Models
#
# code: Matt Siegle, Brendan Connors
# last updated: 28 Jan 2019
# ----------------------------------------------------------------- #
# ----------------------------------------------------------------- #

# set working directory ####
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/GIT_R.file_organization/input files")

# Load functions and libraries for analysis ####
# packages
library(grDevices)
library(plotrix)
library(R2jags) 
library(modeest) 
library(gplots)
require(dplyr)


# Functions for analysis ####

# Posterior summary function
post.summ = function(post.samp, var) {
  post.samp = as.matrix(post.samp)
  
  # if parameter is indexed
  if(substr(var, nchar(var), nchar(var)) == "[") {
    post = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    summ = apply(post, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
    return(summ)
  }
  
  # if parameter is not indexed
  if(substr(var, nchar(var), nchar(var)) != "[") {
    post = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    summ = c(mean = mean(post), sd = sd(post), quantile(post, c(0.5, 0.025, 0.975)))
    return(summ)
  }
}


# Multi-stock simulation function

# ny <- the number of years
# Ro <- the sub-stock recruiment at time zero
# phi <- the expected correlation through time
# mat <- stock-specific maturation schedules
# alpha <- sub-stock productivity (not in log space)
# beta <- sub-stock density depedence 
# sigma.R <- recruitment variation
# U <- finite annual exploitation rate
# pm.yr <- year of simulation that pms start to be calculated over
# Rec <- estimated recruitments from last years of empirical data 
# Spw <- estimated spawers from last years of empirical data
# lst.resid <- estimated recruitment deviation from last year of empirical data

process = function(ny,Ro,phi,mat,U,alpha,beta,sigma.R,Rec,Spw,lst.resid){
  ns = length(Ro) #number of sub-stocks
  m.alpha <- alpha
  m.beta <- beta
  epi = rnorm(ny, sd= sigma.R)
  
  #Build time series of Spawners (S), abundance of returning spawners pre-harvest
  # (N), and the component of the residual that is correlated throught time (v)
  R = t(matrix(0,ns,ny))
  S = R * (1-0)
  v = R; v[,]=0
  #R[1:7,]=t(replicate(7,Ro,simplify=T))*exp(epi[1:7,])
  R[1:3,]=Rec
  N = array(0,dim=c(ny,4,ns))
  Ntot = R; Ntot[,]=0
  H = Ntot; S = Ntot
  S[4:7,] = Spw
  predR = Ntot
  
  # populate first few years with realized states
  R[4,] = exp(log(alpha[]*S[4,]*exp(-beta[]*S[4,])) + phi* lst.resid) * exp(epi[4])
  v[4,] = log(R[4,])-log(alpha[]*S[4,]*exp(-beta[]*S[4,]))
  
  for(i in 5:7){
    R[i,] = exp(log(alpha[]*S[i,]*exp(-beta[]*S[i,])) + phi* v[i-1,]) * exp(epi[i])
    v[i,] = log(R[i,])-log(alpha[]*S[i,]*exp(-beta[]*S[i,]))		
  }
  
  N[4:7,1,]=R[4:7-(3),] * mat[1]
  N[5:7,2,]=R[5:7-(4),] * mat[2]
  N[6:7,3,]=R[6:7-(5),] * mat[3]
  N[7,4,]=R[7-(6),] * mat[4]
  
  # Loop through years of simulation	  
  for(i in (7+1):ny){ 
    N[i,1,1] = R[i-(4),1] * mat[1]
    N[i,2,1] = R[i-(5),1] * mat[2]
    N[i,3,1] = R[i-(6),1] * mat[3]
    N[i,4,1] = R[i-(7),1] * mat[4]
    
    Ntot[i,1] = sum(N[i,,1])
    
    # apply harvest 
    H[i,1] =  U*Ntot[i,1]
    S_exp = Ntot[i,1]-H[i,1] ; S_exp[S_exp<0] = 0
    S[i,1] = S_exp
    
    # predict recruitment
    R[i,1] = alpha[]*S[i,1]*exp(-beta[]*S[i,1]+phi*v[i-1,1]+epi[i])
    predR[i,] = alpha[]*S[i,1]*exp(-beta[]*S[i,1])
    v[i,1] = log(R[i,1])-log(predR[i,1])
    v[v[,1]=='NaN'] <- 0
  }
  
  #Output
  S[S[,]=='NaN'] <- 0
  Ntot[Ntot[,]=='NaN'] <- 0
  p_rg <-ifelse(median(S[(ny-10):ny,])>15000,1,0)
  p_lrp <-ifelse(median(S[(ny-10):ny,])>4000,1,0)
  
  list(S=S[,],N=Ntot[,],survival=as.numeric(v),P=c(p_rg,p_lrp))
}



# Function to sample from posteriors for forward simulations
process.iteration = function(samp) {
  # 1.) extract names
  nms = names(samp)
  
  # 2.) extract elements according to the names and put them into the appropriate data structure
  
  # parameters
  alpha = unname(samp[substr(nms, 1, 5) == "alpha"])
  beta = unname(samp[substr(nms, 1, 5) == "beta"])
  last_resid = unname(samp[substr(nms, 1, 13) == "log.resid.40."])
  phi = unname(samp["phi"])
  sigma_R = unname(samp["sigma.R"])
  mat.sch = c(as.numeric(samp["pi.1."]), as.numeric(samp["pi.2."]), as.numeric(samp["pi.3."]), as.numeric(samp["pi.4."]))
  
  # states
  S = c(as.numeric(samp["S.40."]), as.numeric(samp["S.41."]), as.numeric(samp["S.42."]), as.numeric(samp["S.43."]))
  R = c(as.numeric(samp["R.44."]), as.numeric(samp["R.45."]), as.numeric(samp["R.46."]))
  
  # 3.) create output list
  output = list(
    alpha = as.numeric(alpha),
    beta = as.numeric(beta),
    phi = as.numeric(phi),
    last_resid = as.numeric(last_resid),
    sigma_R = as.numeric(sigma_R),
    mat.sch = mat.sch,
    S = S,
    R = R
  )
  
  # 4.) return output
  return(output)
  
}


# Bayes SR Model Fit #### 

# Load data
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/GIT_R.file_organization/input files")
age.df = read.csv("age_data.csv")
esc.df = read.csv("esc_data.csv")
harv.df = read.csv("harv_data.csv")

age <- age.df %>%
  filter(year != "2017")
esc <- esc.df %>%
  filter(year != "2017")
harv <- harv.df %>%
  filter(year != "2017")


# create CDN aggregate data frames
age.ag <- age %>%
  distinct(year, age_4, age_5, age_6, age_7, gooddata)
age.ag <- as.data.frame(age.ag)

esc.ag <- esc %>%
  group_by(year) %>%
  mutate(spawn = sum(Spawn)) %>%
  ungroup() %>%
  select(year, truecount, spawn) %>%
  distinct(year, spawn, truecount)
esc.ag <- as.data.frame(esc.ag)

harv.ag <- harv %>%
  group_by(year) %>%
  mutate(harvest = sum(Harvest)) %>%
  ungroup() %>%
  select(year, truecount, harvest) %>%
  distinct(year, harvest, truecount)
harv.ag <- as.data.frame(harv.ag)


# create individual substock data frames for data: age, esc, harv
age.YC <- age %>%
  filter(region == "Yukon.Carmacks") # careful about substock naming convention matt Dec 2018 note
age.YlC <- age %>%
  filter(region == "Lower Mainstem")
age.Ym <- age %>%
  filter(region == "Middle Mainstem")
age.YP <- age %>%
  filter(region == "Pelly")
age.YS <- age %>%
  filter(region == "Stewart")
age.YT <- age %>%
  filter(region == "Teslin River")
age.Yu <- age %>%
  filter(region == "Upper Lakes and Mainstem")
age.YWD <- age %>%
  filter(region == "White-Donjek")

esc.YC <- esc %>%
  filter(region == "Yukon.Carmacks")
esc.YlC <- esc %>%
  filter(region == "Lower Mainstem")
esc.Ym <- esc %>%
  filter(region == "Middle Mainstem")
esc.YP <- esc %>%
  filter(region == "Pelly")
esc.YS <- esc %>%
  filter(region == "Stewart")
esc.YT <- esc %>%
  filter(region == "Teslin River")
esc.Yu <- esc %>%
  filter(region == "Upper Lakes and Mainstem")
esc.YWD <- esc %>%
  filter(region == "White-Donjek")

harv.YC <- harv %>%
  filter(region == "Yukon.Carmacks")
harv.YlC <- harv %>%
  filter(region == "Lower Mainstem")
harv.Ym <- harv %>%
  filter(region == "Middle Mainstem")
harv.YP <- harv %>%
  filter(region == "Pelly")
harv.YS <- harv %>%
  filter(region == "Stewart")
harv.YT <- harv %>%
  filter(region == "Teslin River")
harv.Yu <- harv %>%
  filter(region == "Upper Lakes and Mainstem")
harv.YWD <- harv %>%
  filter(region == "White-Donjek")


# Set beta prior for stock ####
# for CDN aggregate
SMAX <- mean(esc.ag$spawn) # mean escapement over time series

# get mean SMAX value by substock region
df_SMAX <- esc %>%
  group_by(region) %>%
  summarize(SMAX = mean(spwn_stock))   # Spawn

# get individual SMAX values for each sub-stock, need to change substock accordingly
SMAX <- as.numeric(df_SMAX[1,2])    # for YC
SMAX <- as.numeric(df_SMAX[2,2])    # for YlC
SMAX <- as.numeric(df_SMAX[3,2])    # for Ym
SMAX <- as.numeric(df_SMAX[4,2])    # for YP
SMAX <- as.numeric(df_SMAX[5,2])    # for YS
SMAX <- as.numeric(df_SMAX[6,2])    # for YT
SMAX <- as.numeric(df_SMAX[7,2])    # for Yu
SMAX <- as.numeric(df_SMAX[8,2])    # for YWD

# turn SMAX into beta prior, need to keep track of substock accordingly
bpmu = 1/SMAX

# set CV on prior, need to keep track of substock accordingly (currently set at a 90% CV)
bptau = 1/(0.9^2)


# Format data for Bayes Model ####

# for CDN aggregate
Y = nrow(age.ag)          # number of calendar years observed
a.min = 4                 # minimum age class in data set
a.max = 7                 # maximum age class in data set
A = a.max - a.min + 1     # number of age classes
nRyrs = Y + A - 1         # number of recruitment years (see model code for details)
years = age.ag[,"year"]

# escapement: assume a 30% observation CV if directly observed, 50% otherwise
S.cv = ifelse(esc.ag[,2] == 1, 0.3, 0.5)
S.obs = esc.ag[,3]

# harvest: assume a 15% observation CV if directly observed, 30% otherwise
C.cv = ifelse(harv.ag[,2] == 1, 0.15, 0.30)
C.obs = harv.ag[,3]

# age composition: assume a ESS of 100 if directly observed, 25 otherwise
ESS = ifelse(age.ag$gooddata == 1, 100, 25) 
X <- age.ag %>%
  select(2:5)
X <- as.matrix(X)
X = t(apply(X, 1, function(x) x/sum(x)))
X = round(apply(X, 2, function(x) x * ESS))
colnames(X) = NULL
n = rowSums(X)        # n is slightly different than ESS because of rounding errors
x=X


# for each individual substock
# need to change age.df by substock accordingly

Y = nrow(age.YC)          # number of calendar years observed
a.min = 4                 # minimum age class in data set
a.max = 7                 # maximum age class in data set
A = a.max - a.min + 1     # number of age classes
nRyrs = Y + A - 1         # number of recruitment years (see model code for details)
years = age.YC[,"year"]


# escapement: assume a 30% observation CV if directly observed, 50% otherwise
# need to change age.df by substock accordingly
S.cv = ifelse(esc.YC[,5] == 1, 0.3, 0.5)
S.obs = esc.YC[,3]



# harvest: assume a 15% observation CV if directly observed, 30% otherwise
# need to change age.df by substock accordingly
C.cv = ifelse(harv.YC[,5] == 1, 0.15, 0.30)
C.obs = harv.YC[,3]



# age composition: assume a ESS of 100 if directly observed, 25 otherwise
# for each substock, need to change dataframe accordingly
ESS = ifelse(age.YC$gooddata == 1, 100, 25) 
X <- age.YC %>%
  select(3:6)
X <- as.matrix(X)
X = t(apply(X, 1, function(x) x/sum(x)))
X = round(apply(X, 2, function(x) x * ESS))
colnames(X) = NULL
n = rowSums(X)        # n is slightly different than ESS because of rounding errors
x=X


#  Bayes model ####
modelFilename = "dep_mod.txt"
cat("
    model {
    # priors for SR portion
    lnalpha ~ dunif(0, 3) 
    #beta ~ dunif(0,10)
    beta ~ dnorm(bpmu,bptau)
    beta.prior ~ dnorm(bpmu,bptau) #for plotting and checking
    tau.R ~ dgamma(0.01,0.01)  # white noise process error      
    phi ~ dunif(-0.99, 0.99)   # autocorrelation coefficient                                              
    log.resid.0 ~ dnorm(0, tau.red)  # starting residual for AR1 process
    
    # Ricker SR with AR1 process on log recruitment residuals for years with brood year spawners
    for (y in (A+a.min):nRyrs) {
    log.R[y] ~ dnorm(log.R.mean.2[y], tau.R)  # true state R is lognormally distributed around the prediction given by SR with AR1
    R[y] <- exp(log.R[y])
    log.R.mean.1[y] <- lnalpha + log(S[y-a.max]) - beta * S[y-a.max]
    log.resid.a[y] <- log.R[y] - log.R.mean.1[y]
    }             
    
    log.R.mean.2[A+a.min] <- log.R.mean.1[A+a.min] + phi * log.resid.0
    
    for (y in (A+a.min+1):nRyrs) {
    log.R.mean.2[y] <- log.R.mean.1[y] + phi * log.resid.a[y-1]
    }
    
    #derived quantities
    tau.red <- tau.R * (1 - phi * phi)
    sigma.red <- 1 / sqrt(tau.red)
    sigma.R <- 1 / sqrt(tau.R)
    alpha <- exp(lnalpha)
    log.resid <- log.resid.a[(A+a.min):nRyrs]
    
    # First `a.max` years of recruits, for which there is no spawner link
    mean.log.R0 ~ dnorm(0, 1E-4) 
    mean.R0 <- exp(mean.log.R0)
    tau.R0 ~ dgamma(0.1,0.1)
    sigma.R0 <- 1/sqrt(tau.R0)
    for (y in 1:a.max) {
    log.R[y] ~ dnorm(mean.log.R0, tau.R0)   
    R[y] <- exp(log.R[y])
    }
    
    # biological reference points: derived quantities
    lnalpha.c <- lnalpha + (sigma.R * sigma.R)/2/(1-phi * phi)
    S.max <- 1/beta
    S.eq <- lnalpha.c * S.max
    S.msy <- S.eq * (0.5 - 0.07 * lnalpha.c)
    U.msy <- lnalpha.c * (0.5 - 0.07 * lnalpha.c)
    
    # Maturity schedule: here we use a common maturation schedule to draw the brood year specific schedules;
    prob[1] ~ dbeta(1,1)
    prob[2] ~ dbeta(1,1)
    prob[3] ~ dbeta(1,1)
    pi[1]<- prob[1]
    pi[2] <- prob[2] * (1 - pi[1])
    pi[3] <- prob[3] * (1 - pi[1] - pi[2])
    pi[4] <- 1 - pi[1] - pi[2] - pi[3]
    
    D.scale ~ dunif(.045,1)
    D.sum <- 1 / (D.scale * D.scale)
    for (a in 1:A) {
    gamma[a] <- D.sum * pi[a]
    for (y in 1:(Y+A-1)) {                                                    
    g[y,a] ~ dgamma(gamma[a],1.0)
    p[y,a] <- g[y,a]/sum(g[y,])
    }
    }
    
    # Calculate the numbers at age matrix as brood year recruits at age*proportion that matured that year
    for (t in 1:Y) {
    for(a in 1:A){
    N.ta[t,a] <- R[t+A-a] * p[t+A-a,a]
    }
    }
    
    ## OBSERVATION SUBMODEL ##
    # multinomial scale sampling
    for (t in 1:Y) {
    for (a in 1:A) {
    q[t,a] <- N.ta[t,a]/N[t]
    }
    x[t,1:A] ~ dmulti(q[t,1:A], n[t])
    }
    
    for (t in 1:Y) {
    # get observation tau's from assumed CV's
    log.sigma.C[t] <- sqrt(log((C.cv[t]^2) + 1))
    log.tau.C[t] <- 1/log.sigma.C[t]^2
    log.sigma.S[t] <- sqrt(log((S.cv[t]^2) + 1))
    log.tau.S[t] <- 1/log.sigma.S[t]^2
    
    # catch model
    U[t] ~ dunif(0.01, 0.99)
    N[t] <- sum(N.ta[t,1:A])
    S[t] <- N[t] * (1 - U[t])
    
    C[t] <- N[t] * U[t]
    log.C[t] <- log(C[t])
    C.obs[t] ~ dlnorm(log.C[t], log.tau.C[t])
    
    # escapement model
    log.S[t] <- log(S[t])
    S.obs[t] ~ dlnorm(log.S[t], log.tau.S[t])
    }
    
    }
    
    ", fill=TRUE, file=modelFilename)

#file.show(modelFilename)


# Jags inputs ####

jags.data = list('Y','a.min','a.max','A','nRyrs','S.cv','S.obs','C.cv','C.obs',
                 'x','n','bpmu','bptau')

jags.parms = c("R", "N", "S", "U", "alpha", "beta", "lnalpha", "phi", "C", "log.resid",
               "log.resid.0","sigma.R", "lnalpha.c", "mean.log.R0", "pi", "q", 
               "mean.R0", "sigma.R0","S.msy", "S.max", "S.eq", "U.msy", "gamma", 
               "D.sum", "p","log.S","beta.prior")


# Fit Model ####

ptm = proc.time()

# change jagsfit.xx to reflect aggregate stock or substock
jagsfit.YC <- jags.parallel(data=jags.data,  parameters.to.save=jags.parms,n.thin=15,
                            n.iter=300000, model.file=modelFilename,n.burnin = 50000,n.chains=6)

endtime = proc.time()-ptm
endtime[3]/60



# Output files ####

# output files for aggregate stock
post.ag = as.mcmc(jagsfit.ag)
mypost.ag = as.matrix(post.ag, chain=F)
write.csv(mypost.ag,"outputs/Yuk_ag_posteriors.Dec_5_2018.csv")


# output files for each sub stock
#YC
post.YC = as.mcmc(jagsfit.YC)
mypost.YC = as.matrix(post.YC, chain=F)
write.csv(mypost.YC,"outputs/Yuk_Carmacks_posteriors.Aug_15_2018.csv")

#YlC
post.YlC = as.mcmc(jagsfit.YlC)
mypost.YlC = as.matrix(post.YlC, chain=F)
write.csv(mypost.YlC,"outputs/Yuk_low_CDN_posteriors.Aug_15_2018.csv")

#Ym
post.Ym = as.mcmc(jagsfit.Ym)
mypost.Ym = as.matrix(post.Ym, chain=F)
write.csv(mypost.Ym,"outputs/Yuk_mainstem_posteriors.Aug_15_2018.csv")

#YP
post.YP = as.mcmc(jagsfit.YP)
mypost.YP = as.matrix(post.YP, chain=F)
write.csv(mypost.YP,"outputs/Yuk_Pelly_posteriors.Aug_15_2018.csv")

#YS
post.YS = as.mcmc(jagsfit.YS)
mypost.YS = as.matrix(post.YS, chain=F)
write.csv(mypost.YS,"outputs/Yuk_Stewart_posteriors.Aug_15_2018.csv")

#YT
post.YT = as.mcmc(jagsfit.YT)
mypost.YT = as.matrix(post.YT, chain=F)
write.csv(mypost.YT,"outputs/Yuk_Teslin_posteriors.Aug_15_2018.csv")

#Yu
post.Yu = as.mcmc(jagsfit.Yu)
mypost.Yu = as.matrix(post.Yu, chain=F)
write.csv(mypost.Yu,"outputs/Yuk_upper_posteriors.Aug_20_2018.csv")

#YWD
post.YWD = as.mcmc(jagsfit.YWD)
mypost.YWD = as.matrix(post.YWD, chain=F)
write.csv(mypost.YWD,"outputs/Yuk_White-Donjek_posteriors.Aug_14_2018.csv")



# read in large, output files so the parameters can be re-generated ####
library(tidyverse)
setwd("~/Dropbox (ESSA Technologies)/EN2438_posterior_outputs")

ag <- read_csv(("Yuk_aggregate_posteriors.Aug_22_2018.csv"))
post.ag <- as.matrix(ag[,2:562]) # remove the 1st column that is added in upon import

YC <- read_csv(("Yuk_Carmacks_posteriors.Aug_15_2018.csv"))
post.YC <- as.matrix(YC[,2:562]) # remove the 1st column that is added in upon import

YlC <- read_csv(("Yuk_low_CDN_posteriors.Aug_15_2018.csv"))
post.YlC <- as.matrix(YlC[,2:562]) # remove the 1st column that is added in upon import

Ym <- read_csv(("Yuk_mainstem_posteriors.Aug_15_2018.csv"))
post.Ym <- as.matrix(Ym[,2:562]) # remove the 1st column that is added in upon import

YP <- read_csv(("Yuk_Pelly_posteriors.Aug_15_2018.csv"))
post.YP <- as.matrix(YP[,2:562]) # remove the 1st column that is added in upon import

YS <- read_csv(("Yuk_Stewart_posteriors.Aug_15_2018.csv"))
post.YS <- as.matrix(YS[,2:562]) # remove the 1st column that is added in upon import

YT <- read_csv(("Yuk_Teslin_posteriors.Aug_15_2018.csv"))
post.YT <- as.matrix(YT[,2:562]) # remove the 1st column that is added in upon import

Yu <- read_csv(("Yuk_upper_posteriors.Aug_20_2018.csv"))
post.Yu <- as.matrix(Yu[,2:562]) # remove the 1st column that is added in upon import

YWD <- read_csv(("Yuk_White-Donjek_posteriors.Aug_14_2018.csv"))
post.YWD <- as.matrix(YWD[,2:562]) # remove the 1st column that is added in upon import



# Model diagnostics and parameter summary ####

# potential scale reduction factor
#need to change dataframe for each substock
gelman.diag(post.ag, multivariate = F)

#need to change dataframe for each substock
R = post.summ(post.YWD, "R[")   # recruitment over time
S = post.summ(post.YWD, "S[")   # spawner abundance over time
N = post.summ(post.YWD, "N[")   # total run size (spawners + catch) over time
U = post.summ(post.YWD, "U[")   # harvest rate over time
resid = post.summ(post.YWD, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP = post.summ(post.YWD, "pi[")   #age composition

alpha = post.summ(post.YWD, "alpha")
beta = post.summ(post.YWD, "beta")
phi = post.summ(post.YWD, "phi")
sigma = post.summ(post.YWD, "sigma.R")
S.msy = post.summ(post.YWD, "S.msy")
S.eq = post.summ(post.YWD, "S.eq")
S.max = post.summ(post.YWD, "S.max")
U.msy = post.summ(post.YWD, "U.msy")

round(rbind(alpha, beta, sigma, phi, S.msy, S.eq, S.max, U.msy, one_over_beta), 2)

mean(mypost[,'beta'])


# save parameters that are substock specific ####
# aggregate CDN stock
R.ag = post.summ(post.ag, "R[")   # recruitment over time
S.ag = post.summ(post.ag, "S[")   # spawner abundance over time
N.ag = post.summ(post.ag, "N[")   # total run size (spawners + catch) over time
U.ag = post.summ(post.ag, "U[")   # harvest rate over time
resid.ag = post.summ(post.ag, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.ag = post.summ(post.ag, "pi[")   #age composition

alpha.ag = post.summ(post.ag, "alpha")
beta.ag = post.summ(post.ag, "beta")
phi.ag = post.summ(post.ag, "phi")
sigma.ag = post.summ(post.ag, "sigma.R")
S.msy.ag = post.summ(post.ag, "S.msy")
S.eq.ag = post.summ(post.ag, "S.eq")
S.max.ag = post.summ(post.ag, "S.max")
U.msy.ag = post.summ(post.ag, "U.msy")

# Yukon Carmacks
R.YC = post.summ(post.YC, "R[")   # recruitment over time
S.YC = post.summ(post.YC, "S[")   # spawner abundance over time
N.YC = post.summ(post.YC, "N[")   # total run size (spawners + catch) over time
U.YC = post.summ(post.YC, "U[")   # harvest rate over time
resid.YC = post.summ(post.YC, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.YC = post.summ(post.YC, "pi[")   #age composition

alpha.YC = post.summ(post.YC, "alpha")
beta.YC = post.summ(post.YC, "beta")
phi.YC = post.summ(post.YC, "phi")
sigma.YC = post.summ(post.YC, "sigma.R")
S.msy.YC = post.summ(post.YC, "S.msy")
S.eq.YC = post.summ(post.YC, "S.eq")
S.max.YC = post.summ(post.YC, "S.max")
U.msy.YC = post.summ(post.YC, "U.msy")

# Yukon low_CDN
R.YlC = post.summ(post.YlC, "R[")   # recruitment over time
S.YlC = post.summ(post.YlC, "S[")   # spawner abundance over time
N.YlC = post.summ(post.YlC, "N[")   # total run size (spawners + catch) over time
U.YlC = post.summ(post.YlC, "U[")   # harvest rate over time
resid.YlC = post.summ(post.YlC, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.YlC = post.summ(post.YlC, "pi[")   #age composition

alpha.YlC = post.summ(post.YlC, "alpha")
beta.YlC = post.summ(post.YlC, "beta")
phi.YlC = post.summ(post.YlC, "phi")
sigma.YlC = post.summ(post.YlC, "sigma.R")
S.msy.YlC = post.summ(post.YlC, "S.msy")
S.eq.YlC = post.summ(post.YlC, "S.eq")
S.max.YlC = post.summ(post.YlC, "S.max")
U.msy.YlC = post.summ(post.YlC, "U.msy")

# Yukon mainstem
R.Ym = post.summ(post.Ym, "R[")   # recruitment over time
S.Ym = post.summ(post.Ym, "S[")   # spawner abundance over time
N.Ym = post.summ(post.Ym, "N[")   # total run size (spawners + catch) over time
U.Ym = post.summ(post.Ym, "U[")   # harvest rate over time
resid.Ym = post.summ(post.Ym, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.Ym = post.summ(post.Ym, "pi[")   #age composition

alpha.Ym = post.summ(post.Ym, "alpha")
beta.Ym = post.summ(post.Ym, "beta")
phi.Ym = post.summ(post.Ym, "phi")
sigma.Ym = post.summ(post.Ym, "sigma.R")
S.msy.Ym = post.summ(post.Ym, "S.msy")
S.eq.Ym = post.summ(post.Ym, "S.eq")
S.max.Ym = post.summ(post.Ym, "S.max")
U.msy.Ym = post.summ(post.Ym, "U.msy")

# Yukon Pelly
R.YP = post.summ(post.YP, "R[")   # recruitment over time
S.YP = post.summ(post.YP, "S[")   # spawner abundance over time
N.YP = post.summ(post.YP, "N[")   # total run size (spawners + catch) over time
U.YP = post.summ(post.YP, "U[")   # harvest rate over time
resid.YP = post.summ(post.YP, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.YP = post.summ(post.YP, "pi[")   #age composition

alpha.YP = post.summ(post.YP, "alpha")
beta.YP = post.summ(post.YP, "beta")
phi.YP = post.summ(post.YP, "phi")
sigma.YP = post.summ(post.YP, "sigma.R")
S.msy.YP = post.summ(post.YP, "S.msy")
S.eq.YP = post.summ(post.YP, "S.eq")
S.max.YP = post.summ(post.YP, "S.max")
U.msy.YP = post.summ(post.YP, "U.msy")

# Yukon Stewart
R.YS = post.summ(post.YS, "R[")   # recruitment over time
S.YS = post.summ(post.YS, "S[")   # spawner abundance over time
N.YS = post.summ(post.YS, "N[")   # total run size (spawners + catch) over time
U.YS = post.summ(post.YS, "U[")   # harvest rate over time
resid.YS = post.summ(post.YS, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.YS = post.summ(post.YS, "pi[")   #age composition

alpha.YS = post.summ(post.YS, "alpha")
beta.YS = post.summ(post.YS, "beta")
phi.YS = post.summ(post.YS, "phi")
sigma.YS = post.summ(post.YS, "sigma.R")
S.msy.YS = post.summ(post.YS, "S.msy")
S.eq.YS = post.summ(post.YS, "S.eq")
S.max.YS = post.summ(post.YS, "S.max")
U.msy.YS = post.summ(post.YS, "U.msy")

# Yukon Teslin
R.YT = post.summ(post.YT, "R[")   # recruitment over time
S.YT = post.summ(post.YT, "S[")   # spawner abundance over time
N.YT = post.summ(post.YT, "N[")   # total run size (spawners + catch) over time
U.YT = post.summ(post.YT, "U[")   # harvest rate over time
resid.YT = post.summ(post.YT, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.YT = post.summ(post.YT, "pi[")   #age composition

alpha.YT = post.summ(post.YT, "alpha")
beta.YT = post.summ(post.YT, "beta")
phi.YT = post.summ(post.YT, "phi")
sigma.YT = post.summ(post.YT, "sigma.R")
S.msy.YT = post.summ(post.YT, "S.msy")
S.eq.YT = post.summ(post.YT, "S.eq")
S.max.YT = post.summ(post.YT, "S.max")
U.msy.YT = post.summ(post.YT, "U.msy")

# Yukon White-Donjek
R.YWD = post.summ(post.YWD, "R[")   # recruitment over time
S.YWD = post.summ(post.YWD, "S[")   # spawner abundance over time
N.YWD = post.summ(post.YWD, "N[")   # total run size (spawners + catch) over time
U.YWD = post.summ(post.YWD, "U[")   # harvest rate over time
resid.YWD = post.summ(post.YWD, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.YWD = post.summ(post.YWD, "pi[")   #age composition

alpha.YWD = post.summ(post.YWD, "alpha")
beta.YWD = post.summ(post.YWD, "beta")
phi.YWD = post.summ(post.YWD, "phi")
sigma.YWD = post.summ(post.YWD, "sigma.R")
S.msy.YWD = post.summ(post.YWD, "S.msy")
S.eq.YWD = post.summ(post.YWD, "S.eq")
S.max.YWD = post.summ(post.YWD, "S.max")
U.msy.YWD = post.summ(post.YWD, "U.msy")

# Yukon upper
R.Yu = post.summ(post.Yu, "R[")   # recruitment over time
S.Yu  = post.summ(post.Yu, "S[")   # spawner abundance over time
N.Yu  = post.summ(post.Yu, "N[")   # total run size (spawners + catch) over time
U.Yu  = post.summ(post.Yu, "U[")   # harvest rate over time
resid.Yu  = post.summ(post.Yu, "log.resid[")   # survival anomalies (or recruitment residuals) over time
PP.Yu  = post.summ(post.Yu, "pi[")   #age composition

alpha.Yu = post.summ(post.Yu, "alpha")
beta.Yu = post.summ(post.Yu, "beta")
phi.Yu = post.summ(post.Yu, "phi")
sigma.Yu = post.summ(post.Yu, "sigma.R")
S.msy.Yu = post.summ(post.Yu, "S.msy")
S.eq.Yu = post.summ(post.Yu, "S.eq")
S.max.Yu = post.summ(post.Yu, "S.max")
U.msy.Yu = post.summ(post.Yu, "U.msy")



# State-space Figures ####
# Figure 5: CDN aggregate stock escapement and harvest plot
head(esc.df)
head(harv.df)


# load exploitation rate
er = read.csv("er.csv")
ex_rate <- as.vector(er$er)

# create spawners dataframe
spw_df <- esc.df %>%
  group_by(year) %>%
  summarise(Spawn = sum(Spawn)) %>%
  mutate(er = ex_rate)

# create harvest dataframe
har_df <- harv.df %>%
  group_by(year) %>%
  summarise(Harvest = sum(Harvest)) %>%
  mutate(er = ex_rate)

# merge spawners and harvest data frame
w <- merge(spw_df, har_df) %>%
  gather(Stock, value, Spawn:Harvest) %>%
  filter(year < 2016)

# FIGURE, CDN aggregate harvest and spawner stock with exploitation rate
ggplot(w, aes(x=year,y=value, fill = Stock)) +
  geom_bar(stat = "identity", width=0.6, position = position_stack(reverse = FALSE)) +
  geom_line(aes(x=year, y= er*200), stat="identity", color="red") +
  xlab("Year") +
  ylab("Count (000s)") +
  scale_fill_manual(values = c("grey76", "grey28")) +
  scale_y_continuous(limits=c(0,200), breaks=c(0,50,100,150,200),
                     sec.axis = sec_axis(~ . /2 , name = "Aggregate harvest rate (%)", 
                                         breaks = seq(0,100,20), labels=seq(0,100,20))) +
  scale_x_continuous(breaks = c(1982,1984,1986,1988,1990,1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=10)) +
  theme(axis.text.y = element_text(size=10)) +
  theme(legend.justification = "center")
ggsave("Fig_5_aggregate.harvest.stock.pdf", width=7, height=5)


# Figure 6: CDN substock escapement and harvest plot 
spw_df <- esc.df %>%
  group_by(region, year) %>%
  summarise(Spawn = sum(Spawn))

har_df <- harv.df %>%
  group_by(region, year) %>%
  summarise(Harvest = sum(Harvest))

ww <- merge(spw_df, har_df) %>%
  gather(Stock, value, Spawn:Harvest) %>%
  filter(year < 2016)

# figure
ww$region2 <- factor(ww$region, levels = c("Lower Mainstem","Stewart","Pelly","White-Donjek","Middle Mainstem",
                                           "Carmacks", "Upper Lakes and Mainstem","Teslin River"))

ggplot(ww, aes(x=year,y=value, fill = Stock)) +
  geom_bar(stat = "identity", width=0.6, position = position_stack(reverse = FALSE)) +
  xlab("Year") +
  ylab("Count (000s)") +
  scale_fill_manual(values = c("grey76", "grey28")) +
  scale_x_continuous(breaks = c(1982,1986,1990,1994,1998,2002,2006,2010,2014)) +
  facet_wrap(~region2,nrow=3, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=8)) +
  theme(legend.justification = "center")
ggsave("FIg_6_substock.spawn.harvest.Jan_03.pdf", width=11, height=9)




# Figure: aggregate stock ricker ####
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/4_State_Space/matt working")
data = read.csv("input_alpha_Usmy.csv")

# get range of abundances for stock; values are recorded in sp_ag vector below
# don't need to re-run
esc.range <- esc.ag %>%
  mutate(min = min(spawn)) %>%
  mutate(max = max(spawn)) %>%
  distinct(region, min, max)

# create vector of abundances for each dataframe (esc.range data is in thousands, hence the /1000)
sp_ag <- as.vector(seq(0,150,length.out=100))

# repeat substock abundance vector 100 times for each value of alpha and beta
spw_ag <- rep(c(sp_ag),(100))

# subsample alpha and beta
plot_df <- data %>%
  filter(region == "Aggregate") %>%
  select(alpha, beta) %>%
  sample_n(100, replace=TRUE)

plot_df2 <- plot_df[rep(seq_len(nrow(plot_df)),each=100),]

# add column of abundances to dataframe
plot_df2$abund <- NA
plot_df2$abund <- spw_ag

# add on 'y' value that has the predicted ricker estimates per abundance value
plot_df2$y <- NA
plot_df2$y <- plot_df2$alpha*plot_df2$abund*exp(-plot_df2$beta*plot_df2$abund)

scen <- rep(1:100, each = 100)

plot_df2$scenario <- NA
plot_df2$scenario <- rep(c(scen),c(1))

# generate points for spawner/recruit pairs
reg <- rep(c("aggregate"),c(35))

S_med <- as.vector(S.ag[3,1:35])
S_low <- as.vector(S.ag[4,1:35])
S_upp <- as.vector(S.ag[5,1:35])

S_df <- data.frame(reg, S_med)
S_df2 <- data.frame(S_df,S_low)
S_df3 <- data.frame(S_df2,S_upp)

R_med <- as.vector(R.ag[3,1:35])
R_low <- as.vector(R.ag[4,1:35])
R_upp <- as.vector(R.ag[5,1:35])

R_df <- data.frame(reg, R_med)
R_df2 <- data.frame(R_df,R_low)
R_df3 <- data.frame(R_df2,R_upp)

df_agg <- cbind(R_df3,S_df3)


# plot figure
ggplot(data=df_agg, aes(x = S_med, y = R_med)) +
  geom_line(data=plot_df2, aes(x = abund, y = y, group=scenario),colour="gray63") +
  geom_point(size=3) +
  geom_pointrange(aes(ymin = R_low, ymax = R_upp)) +
  geom_errorbarh(aes(xmin = S_low, xmax = S_upp)) +
  xlab("Spawners (000s)") +
  ylab("Recruits (000s)") +
  scale_y_continuous(limits=c(0,600), breaks=c(0,100,200,300,400,500,600)) +
  scale_x_continuous(limits=c(0,150), breaks=c(0,25,50,75,100,125,150)) +
  theme_bw() +
  theme(axis.title = element_text(size=12)) +
  theme(axis.text = element_text(size=10))
ggsave(file="Fig_7_rec_per_spawn.pdf", width=7, height = 5)





# Figure: individual substock ricker figure ####
setwd("~/Dropbox (ESSA Technologies)/EN2438 - Yukon Chinook/Data and analyses/4_State_Space/matt working")
data = read.csv("input_alpha_Usmy.csv")

# dataframe manipulation
data2 <- data %>%
  filter(region != "Aggregate")

# get range of abundances for each substock; shown below of vectors
esc.range <- esc.df %>%
  filter(year < 2017) %>%
  group_by(region) %>%
  mutate(min = min(Spawn)) %>%
  mutate(max = max(Spawn)) %>%
  select(region, min, max) %>%
  distinct(region, min, max) %>%
  group_by(region)

# create vector of abundances for each dataframe (esc.range data is in thousands, hence the /1000)
sp_YC <- as.vector(seq(0,30000/1000,length.out=100))
sp_YlC <- as.vector(seq(0,20000/1000,length.out=100))
sp_Ym <- as.vector(seq(0,20000/1000,length.out=100))
sp_YP <- as.vector(seq(0,30000/1000,length.out=100))
sp_YS <- as.vector(seq(0,15000/1000,length.out=100))
sp_YT <- as.vector(seq(0,50000/1000,length.out=100))
sp_Yu <- as.vector(seq(0,10000/1000,length.out=100))
sp_YWD <- as.vector(seq(0,20000/1000,length.out=100))


# repeat substock abundance vector 100 times for each value of alpha and beta
spw_YC <- rep(c(sp_YC),(100))
spw_YlC <- rep(c(sp_YlC),(100))
spw_Ym <- rep(c(sp_Ym),(100))
spw_YP <- rep(c(sp_YP),(100))
spw_YS <- rep(c(sp_YS),(100))
spw_YT <- rep(c(sp_YT),(100)) 
spw_Yu <- rep(c(sp_Yu),(100))
spw_YWD <- rep(c(sp_YWD),(100))

# create full vector of all substock abundances
abund2 <- as.vector(c(spw_YC, spw_YlC, spw_Ym, spw_YP, spw_YS, spw_YT, spw_Yu, spw_YWD))

# subsample alpha and beta
plot_df <- data2 %>%
  group_by(region) %>%
  select(alpha, beta) %>%
  group_by(region) %>%
  sample_n(100, replace=TRUE)

plot_df2 <- plot_df[rep(seq_len(nrow(plot_df)),each=100),]

plot_df2$abund <- NA
plot_df2$abund <- abund2


plot_df2$y <- NA
plot_df2$y <- plot_df2$alpha*plot_df2$abund*exp(-plot_df2$beta*plot_df2$abund)

scen <- rep(1:100, each = 100)

plot_df2$scenario <- NA
plot_df2$scenario <- rep(c(scen),c(8))

# generate dataframe for points for spawner/recruit pairs
region <- rep(c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek"),each=35)

# create dataframe for all the substock S values
S_med <- as.data.frame(cbind(S.YC[3,1:35],S.YlC[3,1:35],S.Ym[3,1:35],S.YP[3,1:35],S.YS[3,1:35],S.YT[3,1:35],S.Yu[3,1:35],S.YWD[3,1:35]))
colnames(S_med) <- c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                     "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek")
S_med2 <- gather(S_med, key = "region", value = "S_med")

S_low <- as.data.frame(cbind(S.YC[4,1:35],S.YlC[4,1:35],S.Ym[4,1:35],S.YP[4,1:35],S.YS[4,1:35],S.YT[4,1:35],S.Yu[4,1:35],S.YWD[4,1:35]))
colnames(S_low) <- c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                     "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek")
S_low2 <- gather(S_low, key = "region", value = "S_low")

S_upp <- as.data.frame(cbind(S.YC[5,1:35],S.YlC[5,1:35],S.Ym[5,1:35],S.YP[5,1:35],S.YS[5,1:35],S.YT[5,1:35],S.Yu[5,1:35],S.YWD[5,1:35]))
colnames(S_upp) <- c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                     "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek")
S_upp2 <- gather(S_upp, key = "region", value = "S_upp")

S_df <- cbind(S_low2, S_upp2, by="region")
S_df2 <- cbind(S_df, S_med2, by = "region")


# create dataframe for all the substock R values
R_med <- as.data.frame(cbind(R.YC[3,1:35],R.YlC[3,1:35],R.Ym[3,1:35],R.YP[3,1:35],R.YS[3,1:35],R.YT[3,1:35],R.Yu[3,1:35],R.YWD[3,1:35]))
colnames(R_med) <- c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                     "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek")
R_med2 <- gather(R_med, key = "region2", value = "R_med")

R_low <- as.data.frame(cbind(R.YC[4,1:35],R.YlC[4,1:35],R.Ym[4,1:35],R.YP[4,1:35],R.YS[4,1:35],R.YT[4,1:35],R.Yu[4,1:35],R.YWD[4,1:35]))
colnames(R_low) <- c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                     "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek")
R_low2 <- gather(R_low, key = "region2", value = "R_low")

R_upp <- as.data.frame(cbind(R.YC[5,1:35],R.YlC[5,1:35],R.Ym[5,1:35],R.YP[5,1:35],R.YS[5,1:35],R.YT[5,1:35],R.Yu[5,1:35],R.YWD[5,1:35]))
colnames(R_upp) <- c("Yuk.Carmacks","Yuk.low.CDN","Yuk.mainstem","Yuk.Pelly",
                     "Yuk.Stewart","Yuk.Teslin","Yuk.upper","Yuk.White.Donjek")
R_upp2 <- gather(R_upp, key = "region2", value = "R_upp")

R_df <- cbind(R_low2, R_upp2, by="region2")
R_df2 <- cbind(R_df, R_med2, by = "region2")


# merge R and S substock dataframes
R_S_df <- cbind(R_df2, S_df2)
R_S_df2 <- subset(R_S_df, select=c("region", "R_low","R_upp","R_med","S_low","S_upp","S_med"))

# plot figure
ggplot(data=R_S_df2, aes(x = S_med, y = R_med, group=region)) +
  geom_line(data=plot_df2, aes(x = abund, y = y, group=scenario),colour="gray63") +
  geom_point(size=2) +
  geom_pointrange(aes(ymin = R_low, ymax = R_upp)) +
  geom_errorbarh(aes(xmin = S_low, xmax = S_upp)) +
  xlab("Spawners (000s)") +
  ylab("Recruits (000s)") +
  facet_wrap(~region, ncol=2, scales = "free") + #
  theme_bw() +
  theme(axis.title = element_text(size=12)) +
  theme(axis.text = element_text(size=10)) +
  theme(strip.background = element_rect(fill="red")) # need to work out how to match color, maybe calling aes?
ggsave(file="Fig_8_rec_per_spawn.pdf", width=7, height = 5)






# Figure: productivity for CDN aggregate run ####
label <- rep(c("mean", "sd", "50%", "2.5%", "97.5%"), 31)
brood.year <- rep(rep(1982:2012,each=5))

resid2_ag <- as.data.frame(resid.ag[,1:31])

resid.df_ag <- gather(resid2_ag, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_agg"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.ag[3],31))

resid.df_ag$lower <- resid.df_ag$`2.5%`
resid.df_ag$upper <- resid.df_ag$`97.5%`
resid.df_ag$median <- resid.df_ag$`50%`

df_full <- resid.df_ag %>%
  mutate(y = median) %>%    # + alpha_adj
  mutate(low = lower) %>%   # + alpha_adj
  mutate(up = upper)        # + alpha_adj

ggplot(df_full, aes(x=brood.year, y = y)) +
  geom_point(size=3) +
  geom_line() +
  geom_ribbon(aes(ymin=low, ymax=up),fill="grey25", alpha=0.2) +
  geom_hline(yintercept = 0, lty = "dotted") +     #geom_hline(data=df_full,aes(yintercept = alpha_adj), lty="dotted") +
  xlab("Brood year") +
  ylab("Productivity index") +
  #scale_x_continuous(breaks=c(1982, 1986, 1990, 1994, 1998, 2002, 2006, 2010)) +
  scale_y_continuous(limits=c(-3,2), breaks=c(-3,-2,-1,0,1,2,3)) +
  theme_classic() +
  theme(axis.title = element_text(size=12)) +
  theme(axis.text = element_text(size=10))
ggsave(file="Fig_9_product_aggr_Oct152018.pdf", width=7,height=5)



# Figure: productivity for all 8 substocks ####
# look in Obj_4_State_Space.R for updated code with color matching
label <- rep(c("mean", "sd", "50%", "2.5%", "97.5%"), 31)
brood.year <- rep(rep(1982:2012,each=5))


# keep relevant data from residuals output
resid2_Yu <- as.data.frame(resid.Yu[,1:31])
resid2_YC <- as.data.frame(resid.YC[,1:31])
resid2_YlC <- as.data.frame(resid.YlC[,1:31])
resid2_Ym <- as.data.frame(resid.Ym[,1:31])
resid2_YP <- as.data.frame(resid.YP[,1:31])
resid2_YS <- as.data.frame(resid.YS[,1:31])
resid2_YT <- as.data.frame(resid.YT[,1:31])
resid2_YWD <- as.data.frame(resid.YWD[,1:31])


# continue to create dataframes of residuals for each substock
resid.df_Yu <- gather(resid2_Yu, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_upper"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.Yu[3],31))

resid.df_YC <- gather(resid2_YC, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_Carmacks"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids )%>%
  mutate(alpha_adj = rep(alpha.YC[3],31))

resid.df_YlC <- gather(resid2_YlC, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_low_CDN"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.YlC[3],31))

resid.df_Ym <- gather(resid2_Ym, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_mainstem"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.Ym[3],31))

resid.df_YP <- gather(resid2_YP, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_Pelly"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.YP[3],31))

resid.df_YS <- gather(resid2_YS, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_Stewart"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.YS[3],31))

resid.df_YT <- gather(resid2_YT, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_Teslin"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.YT[3],31))

resid.df_YWD <- gather(resid2_YWD, key="year",value="resids") %>%
  mutate(stat_type = label) %>%
  mutate(brood.year = brood.year) %>%
  mutate(region = rep(c("Yukon_White-Donjek"),c(155))) %>%
  select(brood.year, region, stat_type, resids) %>%
  spread(stat_type, resids) %>%
  mutate(alpha_adj = rep(alpha.YWD[3],31))

# add values for 95% lower and upper CI and median productivity residuals to each substock dataframe
resid.df_YC$lower <- resid.df_YC$`2.5%`
resid.df_YC$upper <- resid.df_YC$`97.5%`
resid.df_YC$median <- resid.df_YC$`50%`

resid.df_YlC$lower <- resid.df_YlC$`2.5%`
resid.df_YlC$upper <- resid.df_YlC$`97.5%`
resid.df_YlC$median <- resid.df_YlC$`50%`

resid.df_Ym$lower <- resid.df_Ym$`2.5%`
resid.df_Ym$upper <- resid.df_Ym$`97.5%`
resid.df_Ym$median <- resid.df_Ym$`50%`

resid.df_YP$lower <- resid.df_YP$`2.5%`
resid.df_YP$upper <- resid.df_YP$`97.5%`
resid.df_YP$median <- resid.df_YP$`50%`

resid.df_YS$lower <- resid.df_YS$`2.5%`
resid.df_YS$upper <- resid.df_YS$`97.5%`
resid.df_YS$median <- resid.df_YS$`50%`

resid.df_YT$lower <- resid.df_YT$`2.5%`
resid.df_YT$upper <- resid.df_YT$`97.5%`
resid.df_YT$median <- resid.df_YT$`50%`

resid.df_Yu$lower <- resid.df_Yu$`2.5%`
resid.df_Yu$upper <- resid.df_Yu$`97.5%`
resid.df_Yu$median <- resid.df_Yu$`50%`

resid.df_YWD$lower <- resid.df_YWD$`2.5%`
resid.df_YWD$upper <- resid.df_YWD$`97.5%`
resid.df_YWD$median <- resid.df_YWD$`50%`

# bind all the individual substock dataframes together
full.resid <- rbind(resid.df_YC, resid.df_YlC, resid.df_Ym,
                    resid.df_YP, resid.df_YS, resid.df_YT,
                    resid.df_Yu, resid.df_YWD) 

# add columns to 'full' dataframe and change substock names
df_full <- full.resid %>%
  mutate(y = median) %>%    # + alpha_adj
  mutate(low = lower) %>%   # + alpha_adj
  mutate(up = upper) %>%    # + alpha_adj
  mutate(reg = case_when(region == "Yukon_Carmacks" ~ "Carmacks",
                         region == "Yukon_low_CDN" ~ "Lower Mainstem",
                         region == "Yukon_Pelly" ~ "Pelly",
                         region == "Yukon_Stewart" ~ "Stewart",
                         region == "Yukon_White-Donjek" ~ "White-Donjek",
                         region == "Yukon_Teslin" ~ "Teslin River",
                         region == "Yukon_mainstem" ~ "Middle Mainstem",
                         region == "Yukon_upper" ~ "Upper Lakes and Mainstem"))

# list of colors that matches map colors
"#ffc431", orange, lower mainstem
"#ffc7ae", pink, pelly
"#a9fee0", aquamarine, stewart
"#dd8071", red, carmacks
"#f2dc3b", yellow, middle mainstem
"#cadeae", moss green, upper lakes
"#288baf", dark blue, white-donjek
"#d2adea", purple, teslin


# plot productivity residuals
# create color scheme
MyColour <- c("#ffc431", "#a9fee0", "#ffc7ae", "#288baf", "#f2dc3b", "#dd8071", "#cadeae", "#d2adea")

# order substock factors
df_full$reg2 <- factor(df_full$reg, levels = c("Lower Mainstem","Stewart","Pelly","White-Donjek","Middle Mainstem",
                                               "Carmacks", "Upper Lakes and Mainstem","Teslin River"))

# plot
ggplot(df_full, aes(x=brood.year, y = y, colour = reg2), show.legend = F) +
  geom_point(size=2,show.legend = F) +
  scale_color_manual(values = MyColour) +
  geom_line(show.legend = F) +
  geom_ribbon(aes(ymin=low, ymax=up,fill=reg2), show.legend = F, alpha=0.2) +
  scale_fill_manual(values = MyColour) +
  geom_hline(yintercept = 0, lty = "dotted") +     #geom_hline(data=df_full,aes(yintercept = alpha_adj), lty="dotted") +
  xlab("Brood year") +
  ylab("Productivity index") +
  scale_x_continuous(breaks=c(1982, 1986, 1990, 1994, 1998, 2002, 2006, 2010)) +
  scale_y_continuous(limits=c(-3.5,2.5), breaks=c(-3,-2,-1,0,1,2)) +
  facet_wrap(~reg2,nrow=4) +  # ,scales = "free_y"
  theme(legend.position = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size=10))
ggsave(file="Fig_9_prod_all_substocks_JAN_03.pdf", width=11,height=9)





# Figure: alpha posteriors plot for all substocks ####
data = read.csv("outputs/input_alpha_Usmy.csv")
head(data)

ggplot(data, aes(x=region, y=alpha)) +
  geom_boxplot(outlier.shape = NA, size=1.25) +
  xlab("Substock") +
  ylab("Alpha posteriors") +
  geom_vline(xintercept = 1.5, lty="dotted") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))+
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.title = element_text(size=12))
ggsave(file="alpha_substocks_aug23.pdf", width=7,height=5)


# Figure: U.smy plot for all substocks ####
data = read.csv("outputs/input_alpha_Usmy.csv")
head(data)

ggplot(data, aes(x=region, y=U.msy)) +
  geom_boxplot(outlier.shape = NA, size=1.25) +
  #coord_flip() +
  xlab("Substock") +
  ylab("U.msy posteriors") +
  geom_vline(xintercept = 1.5, lty="dotted") +
  scale_y_continuous(limits=c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10))+
  theme(axis.text.y = element_text(size=10)) +
  theme(axis.title = element_text(size=12))
ggsave(file="Usmy_substocks_Aug23.pdf", width=7,height=5)






