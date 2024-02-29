library(deSolve)
library(ggplot2)

#lotka-voltera de base ----

f <- function(N, parmsf = c(r, v)){
  
    with(as.list(parmsf), {
      res_f <- r*(1-N/v) #v est N_max
  
      return(res_f)
  })
}

g <- function(N, P, parmsg = c(b,B,k)){
    with(as.list(parmsg), {
    
      res_g <- b*N/(B+k*P+N)
  
      return(res_g)
    })
}

LV_base <- function(t, ConI, parmsf = c(r, v), parmsg = c(b,B,k), h){ #À compléter!
  
    with(as.list(ConI), {
    # Lotka-voltera
    dN <- f(N, parmsf)*N - g(N,P, parmsg)*P # dN/dt
    dP <- h*g(N,P, parmsg)*P #dP/dt

    # Resultat
    res <- c(dN = dN, dP = dP)
    return(list(res))
    })
}

