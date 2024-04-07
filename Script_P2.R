library(deSolve)
library(ggplot2)
library(phaseR)

#### zone de tests fct LV_SPNR ####
source('fct_LV_SPNR.R')

S0 <- 15
P0 <- 30
N0 <- 45
R0 <- 100
CI <- c(S=S0, P=P0, N=N0, R=R0)

a_N <- 1
a_P <- 1
a_S <- 1
alpha <- 0.5
b_N <- 1
b_P <- 1
b_S <- 1
beta <- 1
delta <- 1
K <- 1
gamma <- 1
m_N <- 0.1
m_P <- 0.1
m_S <- 0.1
phi <- 0.5
psi <- 0.5
r <- 2
para = c(a_N=a_N,a_P=a_P,a_S=a_S,alpha=alpha,b_N=b_N,b_P=b_P,
         b_S=b_S,beta=beta,delta=delta,K=K,gamma=gamma,m_N=m_N,
         m_P=m_P,m_S=m_S,phi=phi,psi=psi,r=r)

LV_soln <- ode(CI, seq(1,30), LV_SPNR, para)
