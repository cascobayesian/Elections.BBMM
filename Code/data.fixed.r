library(Rcpp)
library(microbenchmark)
library(MultiRNG)
library(HMP)
library(VGAM)
library(gtools)
library(MCMCpack)
library(TailRank)
library(DirichletMultinomial)

source("Code/sim.data.r")
sourceCpp("Code/DMM_MCMC_Rcpp.cpp")
source("Code/birth.death.r")
source("Code/model.funcs.r")
source("Code/DMM_MCMC_R.r")

args <- commandArgs(trailingOnly =TRUE)
#file.in <- args[1]
#file.out <- args[2]
file.in <- "Data/PortlandME.Ref/PortlandME_ref_full.csv"
file.out <- "test.RData"
print(c(file.in,file.out))
raw.data <- read.table(file.in,header=TRUE,sep=",")
counts <- as.matrix(raw.data[,-c(1)])
towns <- as.character(raw.data[,1])


S <- 2;
N <- dim(raw.data)[1]
Q <- (dim(raw.data)[2]-1)/2
N <- dim(raw.data)[1]
dataSet <- array(0,dim=c(N,Q,S))
r <- 1:N

for(q in 1:Q)
{
  set <- (2*q-1):(2*q)
  dataSet[1:N,q,] <- counts[r,set]
}

K <- 10
print('initializing...')
model <- init.model.mix(dataSet,K);
print('run initializing...')
start <- proc.time()
model <- mh(model,dataSet,T=20,diff=10,verbose=TRUE)
end <- proc.time()
print(end-start)

T <- 5000
print(model$llk)

outs <- rep(list(),T)
outs[[1]] <- list(model,1);

for(j in 2:T)
{
  print(c(j,model$K,sum(table(model$r)>0),sum(table(factor(model$r,1:model$K))<1)))
  print(model$r)
  out<-birth.death(model,dataSet)
  
  model<-out[[1]]
  model <- mh(model,dataSet,T=20,diff=5,verbose=FALSE)
  outs[[j]] <- list(model,out[[2]])
}

save(file = file.out,list=ls())
