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
