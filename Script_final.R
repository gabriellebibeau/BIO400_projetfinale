library(deSolve)
library(ggplot2)

#Premi√®re figures ----
source('fct_LV5.R')

#Conditions initiales arbitraires
P0 <- 0.15
N0 <- 0.3
R0 <- 1 - P0 - N0
CI <- c(P=P0, N=N0, R=R0)

#ParamËtres selon la figure 1
a         <- 1
a_prime   <- 1 #entre 0 et 1 dans la fig1
alpha     <- 0.5
b         <- 1
b_prime   <- 0.25 #entre 0 et 1 dans la fig1
beta      <- 1
K         <- 1
m         <- 0.5
m_prime   <- 0.5
r         <- 1

parms <- c(a=a,a_prime=a_prime,alpha=alpha,b=b,
          b_prime=b_prime,beta=beta,K=K,m=m,
          m_prime=m_prime,r=r)

LV5_soln <- ode(CI, seq(1,30), LV5, parms)
plot(LV5_soln[,'time'], LV5_soln[,'P'])
plot(LV5_soln[,'N'], LV5_soln[,'time'])
plot(LV5_soln[,'time'], LV5_soln[,'R'])

#Deuxi√®me set de figures ----
source('fct_S8.R')

#Condition initiales arbitraires
P0 <- 0.3
N0 <- 1 - P0
Condition_Initiale <- c(P=P0, N=N0)

#paramËtres selon la Figure 1 et mon imagination
alpha     <- 0.5
b         <- 1
b_prime   <- 1
beta      <- 1
e         <- 1
e_prime   <- 1
I         <- 1
m         <- 0.5
m_prime   <- 0.5

para <- c(alpha=alpha,b=b,b_prime=b_prime,beta=beta,e=e,
          e_prime=e_prime,I=I,m=m,m_prime=m_prime)

S8_soln <- ode(Condition_Initiale, seq(1,30), S8, para)
plot(S8_soln[,'time'], S8_soln[,'P'])
plot(S8_soln[,'N'], S8_soln[,'time'])
