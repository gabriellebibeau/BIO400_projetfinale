#Fonction Schoener avec IGP (8) ----
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

#Zone de test S8 ----

#Condition initiales arbitraires
P0 <- 15
N0 <- 30
Condition_Initiale <- c(P=P0, N=N0)

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

para <- c(alpha=alpha,b=b,b_prime=b_prime,beta=beta,e=e,
          e_prime=e_prime,I=I,m=m,m_prime=m_prime)

test_S8 <- S8(5, Condition_Initiale, para);test_S8
#même chose que l'autre fonction fct...