rm(list=ls())
library(mvtnorm)
library(tidyverse)
library(dFCHMM)

set.seed(2021)

# Initialize parameters
N <- 50 # N - number of subjects
TIME <- 500 # Length of time series
ATTEMPTS <- 5 # Number of random initializations for the EM algorithm
num.states <- 3 # Number of states
A <- rbind(c(.8,.1,.1),
           c(.05,.9,.05),
           c(.1,.05,.85)) # True transition probability matrix
P <- 10 # Dimensionality

# Covariance matrices per state
Sigma <- list()
Sigma[[1]] <- matrix(.5,nrow=P,ncol=P)
diag(Sigma[[1]]) <- 1
Sigma[[2]] <- diag(1,P);
Sigma[[2]][abs(row(Sigma[[2]])-col(Sigma[[2]]))==1] <- .3
Sigma[[3]] <- diag(1,P)

# Simulate state per subject
st.all <- matrix(NA,nrow=N,ncol=TIME)
for(n in 1:N){
  st <- c(1)
  while(length(st) < TIME){
    state <- sample(1:3,1,prob=A[st[length(st)],])
    st <- c(st,state)
  }
  st.all[n,] <- st
}

# Simulate time series
Xt <- array(NA,dim=c(TIME,P,N))
for(n in 1:N){
  for(time in 1:TIME){
    Xt[time,,n] <- rmvnorm(1,rep(0,P),Sigma[[st.all[n,time]]]) # Simulate from zero-mean multivariate normal
  }
}

# Fit the HMM
result.hmm <- fit.hmm(Xt,num.states,ATTEMPTS=ATTEMPTS,tol=1e-5)
# Estimate state sequence
st.hat <- hmm_viterbi(result.hmm) %>% 
  relabel.sequence(new.order = c(1,3,2))

# Obtain pseudo-residuals
resid <- pseudo_resid(Xt,result.hmm)
# QQ plot for the first subject
resid[,1] %>% 
  pchisq(df=P) %>% 
  qnorm() %>%
  qqnorm()
abline(a=0,b=1)