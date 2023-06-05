
rm(list=ls())

library(nimble)

# Data simulation:
set.seed(133)
seed=133

lambda.true=list()
lambda.proposed=list()
lambda.lopes.true=list()
lambda.lopes.quantiles=list()
lambda.poisson=list()

waic.poisson=list()
waic.proposed=list()
waic.lopes.quantiles=list()
waic.lopes.true=list()


# Number of simulated data set
n.data=50

# fixed quantities
n=75
beta=c(0.1,-0.5)
# covariate effects:
x1=rnorm(n,0,1)
x2=rnorm(n,0,1)

K=4
theta.true=list()
theta.proposed=list()
theta.lopes.true=list()
theta.lopes.quantiles=list()
theta.poisson=list()

for(sim in 1:n.data){
E=runif(n,15, 2000)
u=rnorm(n,0, sd=0.1)
theta=rep(NA,n)
lambda=rep(NA,n)
for(i in 1:n){
theta[i]=beta[1]*x1[i]+ beta[2]*x2[i]+u[i]
lambda[i]=E[i]*exp(theta[i])
}
 theta.true[[sim]]=theta
 lambda.true[[sim]]=lambda

# Z generation:
w=round(c(0.31,0.3,0.2,0.2)*n,0)
mu=c(-3,0,3,8)
sigma=c(2,1,2,3)
z=c()
for(i in 1:length(w)){
    z=c(z,rnorm(w[i],mu[i],sd=sigma[i]))
}
#####
# True groups:
#####
z1=scale(z)[,1]+rnorm(n,sd=2)
g=cut(z1,quantile(z1,c(0,0.3,0.6,0.8,1),include.lowest=T))
levels(g)=seq(1,K,1)
g[which.min(z1)]=1
g=as.numeric(g)
g=(K+1)-g
g.true=g

# Groups based on sestiles:
K.quantiles=6
p=quantile(z,seq(1,6,1)/6)
g=cut(z,quantile(z,seq(0,6,1)/6),include.lowest=T)
levels(g)=seq(1,K.quantiles,1)
g=as.numeric(g)
g.lopes=(K.quantiles+1)-g


# TRUE hAI:

epsilon=rep(NA,n)
            ind=diag(K)
            ind[lower.tri(ind)]=1
            hAI.true=matrix(NA,nrow=n,ncol=K)
            for(i in 1:n){
            hAI.true[i,]=ind[g.true[i],]
                }

# hAI in De Oliveira et al. based on sestiles:
ind=diag(K.quantiles)
ind[lower.tri(ind)]=1
hAI=matrix(NA,nrow=n,ncol=K.quantiles)
for(i in 1:n){
hAI[i,]=ind[g.lopes[i],]
    }

# Simulation scheme 2:
# Changing epsilon you can obtain other simulation scenarios
gamma=rep(NA,K.quantiles)
# reporting probabilities:
epsilon=rep(NA,n)
epsilon[g.true==1]=0.98
epsilon[g.true==2]=0.95
epsilon[g.true==3]=0.90
epsilon[g.true==4]=0.80

# data generation
y=rep(NA,n)
sme=rep(NA,n)
for(i in 1:n){
sme[i]=E[i]*exp(theta[i])*epsilon[i]
y[i]=rpois(1,sme[i])
}

######
## Fit the models
######

##############
# Model 1: Lopes  model with true group number of components
##############

source("Lopes-Model.txt")
SC_code_Lopes=Lopes_code


n_regions=n
K=4

a = c(0.00, rep(0, K-1))
b = c(0.04, rep(1, K-1))

SC_constants <- list(
n = 75,
E=E,
x1=x1,x2=x2,
K=K, n_regions=n_regions,a = a,
                  b = b)

SC_data=list(y=y)

SC_constants$hAI=hAI.true

SC_inits1=list(gamma = c(0.1,rep(0.1,K-1)), beta=beta)

SC_inits2=list(beta=beta,gamma = c(0.1,rep(0.01,K-1)))

SC_inits=list(chain1=SC_inits1, chain2=SC_inits2)

SC_model <- nimbleModel(SC_code_Lopes, SC_constants, SC_data, SC_inits)
SC_model$initializeInfo()
SC_compiled_model <- compileNimble(SC_model, resetFunctions = TRUE)


# Set up samplers.
SC_mcmc_conf <- configureMCMC(SC_model,monitors=c('beta','mu','epsilon','theta','lambda','gamma','u'), useConjugacy = TRUE,
enableWAIC = TRUE)

SC_mcmc <- buildMCMC(SC_mcmc_conf)

SC_compiled_mcmc <- compileNimble(SC_mcmc, project = SC_model, resetFunctions = TRUE)

seed=123

it=80000
bu=5000
th=5

SC_samples_Lopes=runMCMC(SC_compiled_mcmc, inits=SC_inits,
                   nchains = 2, nburnin=bu, niter = it, samplesAsCodaMCMC = TRUE,
                   thin=th, summary = FALSE, WAIC = TRUE,
                   setSeed=c(seed*1,seed*2))

samps1=(SC_samples_Lopes$samples$chain1)

samps2=(SC_samples_Lopes$samples$chain2)

#### estimation of epsilon (reporting probability):
# some convergence check
beta.post1=samps1[,grep(colnames(samps2),pattern="beta")]
boxplot(data.frame(beta.post1))
points(beta,col=2)
beta.post2=samps2[,grep(colnames(samps2),pattern="beta")]

lambda.post1=samps1[,grep(colnames(samps2),pattern="lambda")]
lambda.post2=samps2[,grep(colnames(samps2),pattern="lambda")]

lambda.lopes.true[[sim]]=rbind(lambda.post1,lambda.post2)
theta.post1=samps1[,grep(colnames(samps2),pattern="theta")]
theta.post2=samps2[,grep(colnames(samps2),pattern="theta")]

waic.lopes.true[[sim]]=SC_samples_Lopes$WAIC$WAIC

theta.lopes.true[[sim]]=log(rbind(theta.post1,theta.post2))


##########
## Lopes model with groups defined by quantiles
##########

source("Lopes-Model.txt")
SC_code_Lopes=Lopes_code

n_regions=n


SC_constants <- list(
n = 75,
E=E,
x1=x1,x2=x2,
K=K.quantiles, n_regions=n_regions,a = c(a[1], rep(0, K.quantiles-1)), #informative prior only to parameter gamma_1
                  b = c(b[1], rep(1, K.quantiles-1))) #informative prior only to parameter gamma_1

SC_data=list(y=y)

SC_constants$hAI=hAI

SC_inits1=list(gamma = c(0.1,rep(0.1,K.quantiles-1)), beta=beta)


SC_inits2=list(beta=beta,gamma = c(0.1,rep(0.01,K.quantiles-1)))

SC_inits=list(chain1=SC_inits1, chain2=SC_inits2)

SC_model <- nimbleModel(SC_code_Lopes, SC_constants, SC_data, SC_inits)
SC_model$initializeInfo()
SC_compiled_model <- compileNimble(SC_model, resetFunctions = TRUE)


# Set up samplers.
SC_mcmc_conf <- configureMCMC(SC_model,monitors=c('beta','mu','epsilon','theta','lambda','gamma','u'), useConjugacy = TRUE,
enableWAIC = TRUE)

SC_mcmc <- buildMCMC(SC_mcmc_conf)

SC_compiled_mcmc <- compileNimble(SC_mcmc, project = SC_model, resetFunctions = TRUE)

it=15000
bu=5000
th=10
SC_samples_Lopes_quantiles=runMCMC(SC_compiled_mcmc, inits=SC_inits,
                   nchains = 2, nburnin=bu, niter = it, samplesAsCodaMCMC = TRUE,
                   thin=th, summary = FALSE, WAIC = TRUE,
                   setSeed=c(seed*1,seed*2))

samps1=(SC_samples_Lopes_quantiles$samples$chain1)

samps2=(SC_samples_Lopes_quantiles$samples$chain2)

beta.post1=samps1[,grep(colnames(samps1),pattern="bet")]
apply(beta.post1,2,mean)
boxplot(data.frame(beta.post1))
points(beta,col=2,pch=19)

lambda.post1=samps1[,grep(colnames(samps2),pattern="lambda")]
lambda.post2=samps2[,grep(colnames(samps2),pattern="lambda")]
lambda.lopes.quantiles[[sim]]=rbind(lambda.post1,lambda.post2)

waic.lopes.quantiles[[sim]]=SC_samples_Lopes_quantiles$WAIC$WAIC

theta.post1=samps1[,grep(colnames(samps2),pattern="theta")]
theta.post2=samps2[,grep(colnames(samps2),pattern="theta")]
theta.lopes.quantiles[[sim]]=log(rbind(theta.post1,theta.post2))


#######
## Proposed model
#######

source("Proposed-Model.txt")
print(Proposed_code)
SC_code_proposed=Proposed_code

n_regions=n
K.prop=10


SC_constants <- list(a0=0.96,a1=1,
n = 75,
E=E,
x1=x1,x2=x2,
K=K.prop, n_regions=n_regions) #informative prior only to parameter gamma_1

#informative prior only to parameter gamma_1

SC_constants$z1=scale(z,center=T,scale=T)[,1]


SC_data=list(y=y)

#x1 = SC_constants$x1
#x2 =  SC_constants$x2

x1 = SC_constants$x1
x2 = SC_constants$x2

y=SC_data$y

mod0=glm(y/E ~ x1+x2-1, family=poisson("log"))
summary(mod0)

SC_inits1=list(tau_u=0.05, gamma=c(sort(unique(epsilon),decreasing=T),rep(0,K.prop-4)),beta=beta,
eta=c(1,1,1,1,rep(0,K.prop-4)))
SC_inits2=list(tau_u=0.01,gamma=c(sort(unique(epsilon),decreasing=T),rep(0,K.prop-4)),beta=beta,
eta=c(1,1,1,1,rep(0,K.prop-4)))

SC_inits=list(chain1=SC_inits1, chain2=SC_inits2)

SC_model <- nimbleModel(SC_code_proposed, SC_constants, SC_data, SC_inits)

SC_model$initializeInfo()

SC_compiled_model <- compileNimble(SC_model, resetFunctions = TRUE)

# Set up samplers.
SC_mcmc_conf <- configureMCMC(SC_model,monitors=c('beta','mu','epsilon','theta','lambda','g','u','gamma','psi','eta'), useConjugacy = TRUE,
enableWAIC = TRUE)

SC_mcmc <- buildMCMC(SC_mcmc_conf)

SC_compiled_mcmc <- compileNimble(SC_mcmc, project = SC_model, resetFunctions = TRUE)

it=55000
bu=5000
th=10

SC_samples_proposed=runMCMC(SC_compiled_mcmc, inits=SC_inits,
nchains = 2, nburnin=bu, niter = it, samplesAsCodaMCMC = TRUE,
thin=th, summary = FALSE, WAIC = TRUE,setSeed=c(seed*1,seed*2))

                   
samps1=data.frame(SC_samples_proposed$samples$chain1)
samps2=data.frame(SC_samples_proposed$samples$chain2)

beta.post1=samps1[,grep(colnames(samps1),pattern="bet")]
beta.post2=samps2[,grep(colnames(samps1),pattern="bet")]
boxplot(data.frame(beta.post1,beta.post2))
points(rep(beta,2),col=2)

par(mfrow=c(1,2))
plot(beta.post1[,1],type="l")
points(beta.post2[,1],type="l",col=2)
plot(beta.post1[,2],type="l")
points(beta.post2[,2],type="l",col=2)

lambda.post1=samps1[,grep(colnames(samps2),pattern="lambda")]
lambda.post2=samps2[,grep(colnames(samps2),pattern="lambda")]
lambda.proposed[[sim]]=rbind(lambda.post1,lambda.post2)

waic.proposed[[sim]]=SC_samples_proposed$WAIC$WAIC

theta.post1=samps1[,grep(colnames(samps2),pattern="theta")]
theta.post2=samps2[,grep(colnames(samps2),pattern="theta")]

theta.proposed[[sim]]=log(rbind(theta.post1,theta.post2))

########
### Simple Poisson model
########
source("Poisson-Model.txt")
print(Poisson_code)
SC_code_Poisson=Poisson_code

n_regions=n

SC_constants <- list(n = n, E=E,x1=x1,x2=x2)

SC_data=list(y=y)


x1 = SC_constants$x1
x2 =  SC_constants$x2
y=SC_data$y

mod0=glm(y/E ~ x1+x2-1, family=poisson("log"))
summary(mod0)

SC_inits1=list(tau_u=0.25,beta=mod0$coefficients)

SC_inits2=list(tau_u=0.25,beta=mod0$coefficients)

SC_inits=list(chain1=SC_inits1, chain2=SC_inits2)


SC_model <- nimbleModel(SC_code_Poisson, SC_constants, SC_data, SC_inits)
SC_model$initializeInfo()
SC_compiled_model <- compileNimble(SC_model, resetFunctions = TRUE)


# Set up samplers.
SC_mcmc_conf <- configureMCMC(SC_model,monitors=c('beta','mu','epsilon','theta','lambda'), useConjugacy = TRUE,
enableWAIC = TRUE)

SC_mcmc <- buildMCMC(SC_mcmc_conf)
SC_compiled_mcmc <- compileNimble(SC_mcmc, project = SC_model, resetFunctions = TRUE)

# Set up samplers.

## Setting seed ----
#set.seed(123) # Seed for R and NIMBLE, 140921 is used for results in the article.
#seed=123
time <- proc.time()  #Iniciando tempo
it=85000
bu=15000
th=100
SC_samples_Poisson=runMCMC(SC_compiled_mcmc, inits=SC_inits,
                   nchains = 2, nburnin=bu, niter = it, samplesAsCodaMCMC = TRUE,
                   thin=th, summary = FALSE, WAIC = TRUE,
                   setSeed=c(seed*1,seed*2))

samps1=data.frame(SC_samples_Poisson$samples$chain1)
samps2=data.frame(SC_samples_Poisson$samples$chain2)

beta.post1=samps1[,grep(colnames(samps1),pattern="bet")]
apply(beta.post1,2,mean)
boxplot(data.frame(beta.post1))
points(beta,col=2,pch=19)

#### estimation of epsilon (reporting probability):

mu.post1=samps1[,grep(colnames(samps2),pattern="mu")]
mu.post2=samps2[,grep(colnames(samps2),pattern="mu")]


lambda.post1=mu.post1
lambda.post2=mu.post2

waic.poisson[[sim]]=SC_samples_Poisson$WAIC$WAIC

theta.post1=samps1[,grep(colnames(samps2),pattern="theta")]
theta.post2=samps2[,grep(colnames(samps2),pattern="theta")]
theta.poisson[[sim]]=log(rbind(theta.post1,theta.post2))

lambda.poisson[[sim]]=rbind(lambda.post1,lambda.post2)
}


waic.tot=cbind(apply(matrix(unlist(waic.lopes.true),ncol=n.data,byrow=F),1,median),
apply(matrix(unlist(waic.poisson),ncol=n.data,byrow=F),1,median),
apply(matrix(unlist(waic.lopes.quantiles),ncol=n.data,byrow=F),1,median),
apply(matrix(unlist(waic.proposed),ncol=n.data,byrow=F),1,median))
colnames(waic.tot)=c("Lopes true","Poisson","Lopes","Proposed")

simulation.results=list(waic.tot,theta.true,theta.lopes.true,theta.lopes.quantiles,theta.proposed,theta.poisson, lambda.true,lambda.lopes.true,lambda.lopes.quantiles,lambda.proposed,lambda.poisson)
names(simulation.results)=c("waic.tot","theta.true","theta.lopes.true","theta.lopes.quantiles","theta.proposed","theta.poisson", "lambda.true","lambda.lopes.true","lambda.lopes.quantiles","lambda.proposed","lambda.poisson")


saveRDS(simulation.results, file = "SimulationStudy-Scenario2.rds")

### final graphs:



lambda.true=simulation.results[[7]]
lambda.lopes.true=simulation.results[[8]]
lambda.proposed=simulation.results[[10]]
lambda.lopes.quantiles=simulation.results[[9]]
lambda.poisson=simulation.results[[11]]

rmse.lambda.true=list()
rmse.lambda.proposed=list()
rmse.lambda.q=list()
rmse.lambda.pois=list()

for(sim in 1:n.data){
rmse.lambda.true[[sim]]=sqrt((apply(lambda.lopes.true[[sim]],2,mean)-lambda.true[[sim]])^2/lambda.true[[sim]])
rmse.lambda.proposed[[sim]]=sqrt((apply(lambda.proposed[[sim]],2,mean)-lambda.true[[sim]])^2/lambda.true[[sim]])
rmse.lambda.q[[sim]]=sqrt((apply(lambda.lopes.quantiles[[sim]],2,mean)-lambda.true[[sim]])^2/lambda.true[[sim]])
rmse.lambda.pois[[sim]]=sqrt((apply(lambda.poisson[[sim]],2,mean)-lambda.true[[sim]])^2/lambda.true[[sim]])
}

rmse.lambda.true=apply(matrix(unlist(rmse.lambda.true),ncol=n.data),1,mean)[1:n]
rmse.lambda.proposed=apply(matrix(unlist(rmse.lambda.proposed),ncol=n.data),1,mean)[1:n]
rmse.lambda.q=apply(matrix(unlist(rmse.lambda.q),ncol=n.data),1,mean)[1:n]
rmse.lambda.pois=apply(matrix(unlist(rmse.lambda.pois),ncol=n.data),1,mean)[1:n]

par(mfrow=c(1,2))
plot(rmse.lambda.true,rmse.lambda.proposed,xlim=c(0,6),ylim=c(0,6),ylab="Proposed model",xlab="Reference model")
abline(0,1,col="grey",lwd=2,lty=2)
plot(rmse.lambda.q,rmse.lambda.proposed,xlim=c(0,6),ylim=c(0,6),pch=19,ylab="Proposed model",xlab="Competing models")
points(rmse.lambda.pois,rmse.lambda.proposed,pch=19,col="grey")
abline(0,1,col="grey",lwd=2,lty=2)

# on log scale
par(mfrow=c(1,2))
plot(log(rmse.lambda.true),log(rmse.lambda.proposed),ylab="Proposed model",xlab="Reference model")
abline(0,1,col="grey",lwd=2,lty=2)
plot(log(rmse.lambda.q),log(rmse.lambda.proposed),pch=19,ylab="Proposed model",xlab="Competing models")
 points(log(rmse.lambda.pois),log(rmse.lambda.proposed),pch=19,col="grey")
abline(0,1,col="grey",lwd=2,lty=2)
abline(0,1,col="grey",lwd=2,lty=2)






