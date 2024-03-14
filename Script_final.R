library(deSolve)
library(ggplot2)

#lotka-voltera de base (1) ----
#fonctions à définir : a', alpha et a (functional response);  phi (recrutement ressource de base);
  #fonction de croissance des proies
  f <- function(N, parmsf = c(r, v)){
  
      with(as.list(parmsf), {
        res_f <- r*(1-N/v) #v est N_max
  
        return(res_f)
    })
  }

  #fonction du taux de consommation des prédateurs
  g <- function(N, P, parmsg = c(b,k)){
      with(as.list(parmsg), {
    
        res_g <- b*N/(k*P+N) #nous assumons que B=0, donc qu'il n'y a pas de constante de saturation
  
        return(res_g)
      })
  }

  #Fonction lotka voltera de base
  LV_base <- function(t, ConI, parmsf = c(r, v), parmsg = c(b,k), parmsLV = c(e,mu)){
  
      with(as.list(ConI, parmsLV), {
      # Lotka-voltera
      dN <- f(N, parmsf)*N - g(N,P, parmsg)*P # dN/dt
      dP <- e*g(N,P, parmsg)*P - mu*P #dP/dt

      # Resultat
      res <- c(dN = dN, dP = dP)
      return(list(res))
      })
  }

#Zone de test LV_base ----
t <- 10 #aléatoire...

#Conditions initiales arbitraires en proportion
N0 <- 0.7 
P0 <- 0.3
Condition_Initiale <- c(N=N0, P=P0)

#Paramètres selon la figure 1
  #f()
  r <- 1
  v <- 1
  parametres_f <- c(r=r,v=v)

  #g()
  b <- 2
  k <- 1
  parametres_g <- c(b=b,k=k)

  #LV_base()
  e <- 0.5
  mu <- 0.6
  parametres_LV <- c(e=e,mu=mu)

#test fct
test_LV_base <- LV_base(t, Condition_Initiale, parametres_f, parametres_g, parametres_LV)
#donne un nombre négatif pour dN, peut être flaw, à vérif éventuellement si la suite est erroné aussi


