remove(list=ls())

library(nimble)
library(coda)

# install.packages("nimble", repos = "http://r-nimble.org", type = "source")

# importing the data

df = readRDS("CKD-PseudoData.rds")
names(df)
# "SC_constants" "SC_data"
# df$"SC_constants": it contains all constants to be used in the nimble code
# df$SC_data: it contains the observed counts

# Nimble requires 3 main ingredients:
# the nimble code specifying the model (Proposed-Model-CKD.txt);
# a specification of the constants used in the code (df$SC_constants);
# a specification of the data (df$SC_data).


# Specifying the spatial structure:

# n_adj: number of neighbours
# adj: neighboring structure
# l_adj: dimension of spatial structure
# weights: weights to be specified for each neighboring specification
# n_regions: number of area involved in the spatial structure
n_adj=df$SC_constants$n_adj
adj =df$SC_constants$adj
l_adj=length(adj)
weights=rep(1,l_adj)
n_regions=length(n_adj)

# Number of clusters involved in the stick-breaking prior
K=20

# List of constants:
# m: number of areas
# E: expected counts
# X1, X2, X3 and X4: explanatory variables in Equation (3.3b) of the main paper
# Z1: variable in Equation (3.3e)

SC_constants <- list(m = df$SC_constants$n,
E=df$SC_constants$E,
x1=df$SC_constants$x1,
x2=df$SC_constants$x2,
x3=df$SC_constants$x3,
x4=df$SC_constants$x4,
z1=df$SC_constants$z1,
K=df$SC_constants$K,
  adj=adj,
  n_adj=n_adj,
  l_adj=l_adj,
  n_regions=n_regions,
  a = c(0.00, rep(0, K-1)), #informative prior only to parameter gamma_1
  b = c(0.05, rep(1, K-1)))
  
# y: response variable in Equation (3.3a)
y=df$SC_data$y

SC_data=list(y=y)

# loading the model specification:
source("Proposed-Model-CKD.txt")
print(ProposedSpatial_code)
SC_code=ProposedSpatial_code

# Specifying initial values for MCMC (2 chains: SC_inits1 and SC_inits2)
# If you want to use more than 2 chains, you should
# specify 3 lists of starting values

SC_inits1=list(nu.s=0.25, nu = 1, gamma = c(0.005,rep(0.015,K-1)),
              beta=c(0,0,0,0))

SC_inits2=list(nu.s=0.35, nu = 1, gamma = 1-c(0.015,rep(0.025,K-1)),
beta=c(0,0,0,0))
               
SC_inits=list(chain1=SC_inits1, chain2=SC_inits2)

# loading the model components: code, constants, data and inits

SC_model <- nimbleModel(SC_code, SC_constants, SC_data, SC_inits)


# Note: some error messages are due to the fact that some variables
# are not initialized as stressed "[Note] Any error reports that follow may simply reflect missing values in model variables."

SC_model$initializeInfo()
SC_compiled_model <- compileNimble(SC_model, resetFunctions = FALSE)

# Set up samplers: selecting parameters to be monitored

SC_mcmc_conf <- configureMCMC(SC_model,
                              monitors=c('beta','gamma','eta','g','psi',
                                         'epsilon','theta','s','u'), useConjugacy = TRUE,
enableWAIC = TRUE)

SC_mcmc <- buildMCMC(SC_mcmc_conf)
SC_compiled_mcmc <- compileNimble(SC_mcmc, project = SC_model, resetFunctions = FALSE)

# Running the model
time <- proc.time()

# iterations: number of MCMC iterations
# burnin: number of MCMC iterations to be removed
# thin: thinning rate of MCMC samples

iterations=110000
burn=10000
thin=100

# number of iterations to be increased (pay attention to the execution time!)
# choose your seed:
seed=123
SC_samples=runMCMC(SC_compiled_mcmc, inits=SC_inits,
                       nchains = 2, nburnin=burn, niter = iterations,
                       samplesAsCodaMCMC = TRUE, thin=thin,
                       summary = FALSE, WAIC = TRUE,
                       setSeed=c(seed*1,seed*2))

samps1=data.frame(SC_samples$samples$chain1)
samps2=data.frame(SC_samples$samples$chain2)
print(SC_samples$WAIC)

# samples of parameters to be monitored
# are stored in samps1 and samps2

# Look at the posterior distribution of beta parameters

beta.post1=samps1[,grep(colnames(samps1),pattern="beta")]
colnames(beta.post1)=c("DemDens","Age","Unempl.","House ownership")
boxplot(data.frame(beta.post1),range=0,las=2,cex.axis=0.8)
abline(h=0,col=2,lty=2,lwd=2)

beta.post2=samps2[,grep(colnames(samps2),pattern="bet")]
colnames(beta.post2)=c("DemDens","Age","Unempl.","House ownership")
boxplot(data.frame(beta.post2),range=0,names=colnames(beta.post2),cex.axis=0.8,las=2)
abline(h=0,col="grey",lwd=2,lty=2)

# convergence check
par(mfrow=c(3,2))
nomi=paste("beta",c(1,2,3,4,5))
plot(beta.post1[,1],col=2,type="l",xlab="Iteration",main=expression(beta[1]),ylab="")
points(beta.post2[,1],col=3,type="l")
plot(beta.post1[,2],col=2,type="l",xlab="Iteration",main=expression(beta[2]),ylab="")
points(beta.post2[,2],col=3,type="l")
plot(beta.post1[,3],col=2,type="l",xlab="Iteration",main=expression(beta[3]),ylab="")
points(beta.post2[,3],col=3,type="l")
plot(beta.post1[,4],col=2,type="l",xlab="Iteration",main=expression(beta[4]),ylab="")
points(beta.post2[,4],col=3,type="l")


#potential scale reduction factor

b1=mcmc(beta.post1)
b2=mcmc(beta.post2)
combinedchains = mcmc.list(b1, b2)
psrf.beta=gelman.diag(combinedchains,transform=T)$psrf


##### Gamma: reporting probabilities

gamma.post1=samps1[,grep(colnames(samps2),pattern="gamma")]
gamma.post2=samps2[,grep(colnames(samps2),pattern="gamma")]
boxplot(gamma.post1,range=0)
boxplot(gamma.post2,range=0)

g1=mcmc(gamma.post1)
g2=mcmc(gamma.post2)
combinedchains = mcmc.list(g1, g2)
psrf.gamma=gelman.diag(combinedchains,transform=T)$psrf

### Relative risks: theta
theta.post1=samps1[,grep(colnames(samps2),pattern="theta")]
theta.post2=samps2[,grep(colnames(samps2),pattern="theta")]
m=ncol(theta.post1)

boxplot(theta.post1, main="Relative Risk",xlab=expression(theta),axes=F,
range=0,ylim=c(0,3))
axis(2,seq(0,2,0.5))
axis(1,seq(1,m,9))
box()

boxplot(theta.post2, main="Relative Risk",xlab=expression(theta),axes=F,
range=0,ylim=c(0,3))
axis(2,seq(0,3,0.5))
axis(1,seq(1,m,9))
box()

t1=mcmc(theta.post1)
t2=mcmc(theta.post2)
combinedchains = mcmc.list(t1, t2)
psrf.theta=gelman.diag(combinedchains,transform=T)$psrf

# combind psrf for model parameters:
n.par=length((c(psrf.beta[,1],psrf.gamma[,1],psrf.theta[,1])))
plot(c(psrf.beta[,1],psrf.gamma[,1],psrf.theta[,1]),type="b",
xlab="Model parameters",ylab="PSRF",pch=19,cex=0.6,axes=F,
ylim=c(1,1.3))
axis(2,seq(1,1.3,0.05))
axis(1,seq(1,n.par,1),
tick=F,labels=F)
box()
abline(h=1.05,lwd=2,lty=2,col="grey")

# Maps:
# loading longitude, latitude and shape of the Apulia region:

library(ggplot2)
geo.info=readRDS("geoinfo.rds")

Relative_risk=apply(theta.post1,2,mean)

# Map with colors: notice that dark colors correpond to small
# relative risks.

t1 = ggplot(data=geo.info) +
  geom_sf(aes(fill = Relative_risk)) +
  scale_fill_viridis_c(option="inferno") +
ggtitle(expression(paste("Relative risks ", hat(theta))),
           subtitle = "Apulia Region, n = 258 Cities")

