#Fonction lotka voltera avec IGP et prédateur suprême ----
LV_SPNR <- function(t, ConI, parmsS = c(a_N,a_P,a_S,alpha,b_N,b_P,
                                        b_S,beta,delta,K,gamma,m_N,
                                        m_P,m_S,phi,psi,r)){
  
  with(as.list(ConI, parmsS), {
    # Lotka-voltera
    dS <- S*(b_S*a_S*R + phi*delta*N + psi*gamma*P - m_S) # dS/dt
    dP <- P*(b_P*a_P*R + beta*alpha*N - gamma*S - m_P) # dP/dt
    dN <- N*(a_N*b_N*R- alpha*P - delta*S - m_N) #dN/dt
    dR <- R*(r*(1- (R/K)) - a_N*N - a_P*P - a_S*S) #dR/dt
    
    # Resultat
    res <- c(dS = dS, dP = dP, dN = dN, dR = dR)
    return(list(res))
  })
}
