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
