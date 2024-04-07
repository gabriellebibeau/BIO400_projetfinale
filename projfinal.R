library(deSolve)
library(ggplot2)
library(phaseR)

#Figure 1
#Fonction
LV5 <- function(t, ConI, parms5 = c(a,a_prime,alpha,b,b_prime,
                                    beta,K,m,m_prime,r)){
  
  with(as.list(ConI, parms5), {
    # Lotka-voltera
    dP <- P*(b_prime*a_prime*R + beta*alpha*N - m_prime) # dP/dt
    dN <- N*(a*b*R - m - alpha*P) #dN/dt
    dR <- R*(r*(1- (R/K)) - a*N - a_prime*P) #dR/dt
    
    # Resultat
    res <- c(dP = dP, dN = dN, dR = dR)
    return(list(res))
  })
}


#Conditions initiales
P0 <- 0.3333
N0 <- 0.3333
R0 <- 1 - P0 - N0
CI_LV5 <- c(P=P0, N=N0, R=R0)

#Paramètres
parametre_LV5 <- c(a         <- 1,
           a_prime   <- 0,
           alpha     <- 0.5,
           b         <- 1,
           b_prime   <- 0,
           beta      <- 1,
           K         <- 1,
           m         <- 0.5,
           m_prime   <- 0.5,
           r         <- 5)

#Solution
LV5_sol <- ode(y=CI_LV5, times= seq(1,100), func= LV5, parms= parametre_LV5)

#Les équations n'arrivent pas à génerer de l'instabilité

#Résultat
figure_1 <- matrix(0.5, 11, 11)

lim <- 0.01
steps <- 100
i <- 10
j <- 10


for(i in 0:10){
  for(j in 0:10){
    
    parametre_LV5 <- c(a <- 1,
               a_prime <- i/10,
               alpha     <- 0.5,
               b         <- 1,
               b_prime <- j/10,
               beta      <- 1,
               K         <- 1,
               m         <- 0.5,
               m_prime   <- 0.5,
               r         <- 1)

    LV5_sol <- ode(y=CI_LV5, times= seq(1,steps), func= LV5, parms= parametre_LV5)
    
    if(LV5_sol[steps,'P'] < lim & LV5_sol[steps,'N'] > lim){ #IGPredator excluded (cyan)
      figure_1[i+1,j+1] <- 0.5
      
    }else if(LV5_sol[steps,'P'] > lim & LV5_sol[steps,'N'] < lim){ #IGPrey excluded (rose)
      figure_1[i+1,j+1] <- 1.5
      
    }else if(LV5_sol[steps,'P']  > lim & LV5_sol[steps,'N']  > lim){ #Stable coexistence (vert)
      figure_1[i+1,j+1] <- 2.5
      
    }else{                                                         #Unstable (brun)
      figure_1[i+1,j+1] <- 3.5
      
    }
  }
}

image (figure_1, col = c("cyan2", "salmon", "green", "brown4"), breaks = 0:4, xlab = 'b_prime' , ylab = 'a_prime')
legend(x = 'bottomleft', legend = c('Prédateur IG exclus', 'Proie IG exclus','Coexistance stable', 'Instable'), fill = 1:4)

ggplot() +
  geom_line(aes(LV5_sol[,'time'], LV5_sol[,'P']), color = 'red') +
  geom_line(aes(LV5_sol[,'time'], LV5_sol[,'N']), color = 'blue') +
  geom_line(aes(LV5_sol[,'time'], LV5_sol[,'R']), color = 'green') +
  labs(title = "Modèle Lotka-Volterra", x = "Temps écoulé (pas de temps)", y = "Nombre d'individus") +
  theme_minimal()


#Figure 2

#Fonction
S8 <- function(t, ConI, parms8 = c(alpha,b,b_prime,
                                   beta,e,e_prime,I,m,m_prime)){
  
  with(as.list(ConI, parms8), {
    # Schoener
    dP <- P*((b_prime*e_prime*I)/(e_prime*P + e*N) + beta*alpha*N - m_prime) # dP/dt
    dN <- N*((b*e*I)/(e_prime*P + e*N) - m - alpha*P) #dN/dt
    
    # Resultat
    res <- c(dP = dP, dN = dN)
    return(list(res))
  })
}

#Condition initiales arbitraires
P0 <- 0.3
N0 <- 0.3
CI_S8 <- c(P=P0, N=N0)

#paramètres selon la Figure 1 et mon imagination
alpha     <- 0.5
b         <- 1
b_prime   <- 1
beta      <- 1
e         <- 1
e_prime   <- 1
I         <- 1
m         <- 0.5
m_prime   <- 0.5

parametre_S8 <- c(alpha=alpha,b=b,b_prime=b_prime,beta=beta,e=e,
          e_prime=e_prime,I=I,m=m,m_prime=m_prime)


#Solution
S8_sol <- ode(y=CI_S8, times= seq(1,100), func= S8, parms= parametre_S8)

ggplot() +
  geom_line(aes(S8_sol[,'time'], S8_sol[,'P']), color = 'red') +
  geom_line(aes(S8_sol[,'time'], S8_sol[,'N']), color = 'blue') +
  labs(title = "Modèle Schoener", x = "Temps écoulé (pas de temps)", y = "Nombre d'individus") +
  theme_minimal()








############ Partie 2

SuperPred <- function(t, ConI, parms = c(a,a_prime,alpha,b,b_prime,
                                    beta,K,m,m_prime,r, 
                                    a_doubleprime, gamma, delta, b_doubleprime, phi, psi, m_doubleprime)){
  
  with(as.list(ConI, parms), {
    # Lotka-voltera
    dS <- S*(b_doubleprime*a_doubleprime*R + phi*delta*N + psi*gamma*P - m_doubleprime) #dS/dt
    dP <- P*(b_prime*a_prime*R + beta*alpha*N - m_prime - gamma*S) # dP/dt
    dN <- N*(a*b*R - m - alpha*P - delta*S) #dN/dt
    dR <- R*(r*(1- (R/K)) - a*N - a_prime*P - a_doubleprime*S) #dR/dt
    
    # Resultat
    res <- c(dS = dS, dP = dP, dN = dN, dR = dR)
    return(list(res))
  })
}


#Conditions initiales
S0 <- 0.2
P0 <- 0.2
N0 <- 0.2
R0 <- 1 - P0 - N0 - S0
CI_SuperPred <- c(S=S0, P=P0, N=N0, R=R0)


#Paramètres
parametre_SuperPred <- c(a   <- 1,
                   a_prime   <- 0.6,
                   alpha     <- 1,
                   b         <- 1,
                   b_prime   <- 0.6,
                   beta      <- 1,
                   K         <- 1,
                   m         <- 0.5,
                   m_prime   <- 0.5,
                   r         <- 1,
                   a_doubleprime <- 1,
                   gamma     <- 0.1,
                   delta     <- 1,
                   b_doubleprime <- 1,
                   phi       <- 1,
                   psi       <- 0.5, 
                   m_doubleprime <- 0.2)

#Solution
SuperPred_sol <- ode(y=CI_SuperPred, times= seq(1,100), func= SuperPred, parms= parametre_SuperPred)

ggplot() +
  geom_line(aes(SuperPred_sol[,'time'], SuperPred_sol[,'S']), color = 'yellow') +
  geom_line(aes(SuperPred_sol[,'time'], SuperPred_sol[,'P']), color = 'red') +
  geom_line(aes(SuperPred_sol[,'time'], SuperPred_sol[,'N']), color = 'blue') +
  geom_line(aes(SuperPred_sol[,'time'], SuperPred_sol[,'R']), color = 'green') +
  labs(title = "Modèle Lotka-Volterra avec super-prédateur", x = "Temps écoulé (pas de temps)", y = "Nombre d'individus") +
  theme_minimal()

#Résultat
figure_1 <- matrix(0.5, 11, 11)

lim <- 0.05
steps <- 100

for(i in 0:10){
  for(j in 0:10){
    
    parametre_SuperPred <- c(a   <- 1,
                             a_prime   <- 0.8,
                             alpha     <- 1,
                             b         <- 1,
                             b_prime   <- 0.7,
                             beta      <- 1,
                             K         <- 1,
                             m         <- 0.5,
                             m_prime   <- 0.5,
                             r         <- 1,
                             a_doubleprime <- i/10,
                             gamma     <- 0.4,
                             delta     <- 0.4,
                             b_doubleprime <- j/10,
                             phi       <- 0.8,
                             psi       <- 0.8, 
                             m_doubleprime <- 0.2)
    
    SuperPred_sol <- ode(y=CI_SuperPred, times= seq(1,100), func= SuperPred, parms= parametre_SuperPred)
    
    if(SuperPred_sol[steps,'S'] > 0 & SuperPred_sol[steps,'P'] > lim & LV5_sol[steps,'N'] > lim){ #1 (red)
      figure_1[i+1,j+1] <- 0.5
      
    }else if(SuperPred_sol[steps,'S'] < 0 & SuperPred_sol[steps,'P'] > lim & LV5_sol[steps,'N'] > lim){ #2 (yellow)
      figure_1[i+1,j+1] <- 1.5
      
    }else if(SuperPred_sol[steps,'S'] > 0 & SuperPred_sol[steps,'P'] < lim & LV5_sol[steps,'N'] > lim){ ##3 (blue)
      figure_1[i+1,j+1] <- 2.5
      
    }else if(SuperPred_sol[steps,'S'] > 0 & SuperPred_sol[steps,'P'] > lim & LV5_sol[steps,'N'] < lim){ #4 (green)
      figure_1[i+1,j+1] <- 3.5
      
    }else if(SuperPred_sol[steps,'S'] < 0 & SuperPred_sol[steps,'P'] < lim & LV5_sol[steps,'N'] > lim){ #5 (pink)
      figure_1[i+1,j+1] <- 4.5
      
    }else if(SuperPred_sol[steps,'S'] > 0 & SuperPred_sol[steps,'P'] < lim & LV5_sol[steps,'N'] < lim){ #6 (orange)
      figure_1[i+1,j+1] <- 5.5
      
    }else if(SuperPred_sol[steps,'S'] < 0 & SuperPred_sol[steps,'P'] > lim & LV5_sol[steps,'N'] < lim){ #7 (cyan)
      figure_1[i+1,j+1] <- 6.5
      
    }
  }
}


image (figure_1, col = c("red", "yellow", "blue", "green", "pink", "orange", "cyan" ), breaks = 0:7)


#1. Super pred: O
#   Pred: O
#   Prey: O

#2. Super pred: X
#   Pred: O
#   Prey: O

#3. Super pred: O
#   Pred: X
#   Prey: O

#4. Super pred: O
#   Pred: O
#   Prey: X

#5. Super pred: X
#   Pred: X
#   Prey: O

#6. Super pred: O
#   Pred: X
#   Prey: X

#7. Super pred: X
#   Pred: O
#   Prey: X
