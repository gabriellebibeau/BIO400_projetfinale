#Fonction lotka voltera avec IGP (5) ----
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

#Zone de test LV5 ----

#Conditions initiales arbitraires en proportion
P0 <- 15
N0 <- 30
R0 <- 100 - P0 - N0
Condition_Initiale <- c(P=P0, N=N0, R=R0)

#Paramètres selon la figure 1
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

para <- c(a=a,a_prime=a_prime,alpha=alpha,b=b,
          b_prime=b_prime,beta=beta,K=K,m=m,
          m_prime=m_prime,r=r)

#test fct
test_LV5 <- LV5(5, Condition_Initiale, para);test_LV5
#souvent nombres négatifs, mais j'ai l'impression que c'est juste un 
#problème de paramètres... À garder en tête pour la suite