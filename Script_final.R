library(deSolve)
library(ggplot2)
library(phaseR)



#PLAN
#
#1. Créer un algorithme pour que le resultat de l'équation avec les paramètres varie en fonction de a prime et b prime (pour chacunes des combinaisons
# aprime bprime du quadrant (0-1))

# When Stable = IG prey out

#Stable

#R* > 0, si a'*m*beta + r*alpha*beta - a*m' > 0
#Après manipulations algébriques et en assumant que a'*m*beta + r*alpha*beta - a*m' = 0, on arrive au résultat suivant:

#a_prime <- (a*m_prime - r*alpha*beta)/(m*beta)

#a' = 0. Donc, logiquement, a' > 0

#N* > 0 si m*a'*alpha + K*a*a'*b*m' - K*(a'^2)*b'*m - K*a'*b'*r*alpha > 0

#P* > 0 si K*a*b*r*alpha*beta + K*a*a'*b'*m - K*(a^2)*b*m' - m*r*alpha*beta > 0
#Après manipulations algébriques et en assumant que K*a*b*r*alpha*beta + K*a*a'*b'*m - K*(a^2)*b*m' - m*r*alpha*beta = 0, 
#on arrive au résultat suivant:

#a_prime <- (K*(a^2)*b*m_prime + m*r*alpha*beta - K*a*b*r*alpha*beta)/(K*a*b_prime*m)

#b_prime <- (K*(a^2)*b*m_prime + m*r*alpha*beta - K*a*b*r*alpha*beta)/(K*a*a_prime*m)


#(a_prime*m*beta + r*alpha*beta - a*m_prime > 0) &
#  (m*a_prime*alpha + K*a*a_prime*b*m_prime - K*(a_prime^2)*b_prime*m - K*a_prime*b_prime*r*alpha) > 0 &
# (K*a*b*r*alpha*beta + K*a*a_prime*b_prime*m - K*(a^2)*b*m_prime - m*r*alpha*beta > 0))


# 5 equilibria:
# (i): All species at 0 densitiy. If K > 0, unstable
# (ii): Basal ressource at K, but IGpredator and IGprey absent. Unstable if K > m/ab, or K > m'/a'b
# (iii): Resource and IGprey are present, with e.d. of m/ab and (r/a)(1-m/abK). Unstable if a'b'(m/ab)+alphabeta(r/a)(1-m/abK)-m' > 0
# (iv): Ressorce and IGpredator are present, e.d. m'/a'b' and (r/a')(1-m'/a'b'K). 
# (v): all species present, with densities x


param <- c(a         <- 1,
           a_prime   <- 1, #Variable
           alpha     <- 0.5,
           b         <- 1,
           b_prime   <- 1, #Variable
           beta      <- 1,
           K         <- 1,
           m         <- 0.5,
           m_prime   <- 0.5,
           r         <- 1)


#Three-species equilibrium
  

D <- K*a*a_prime*(b*beta - b_prime) + r*alpha*beta

R_star <- function (param){
  D <- K*a*a_prime*(b*beta - b_prime) + r*alpha*beta
  R <- K*(r*alpha*beta + a_prime*m*beta - a*m_prime)/D
  return(R)
}

N_star <- function (param){
  D <- K*a*a_prime*(b*beta - b_prime) + r*alpha*beta
  N <- (K*a*a_prime*b*m_prime + m_prime*r*alpha - K*(a_prime^2)*b_prime*m - K*a_prime*b_prime*r*alpha)/(alpha*D)
  return(N)
}

P_star <- function (param){
  D <- K*a*a_prime*(b*beta - b_prime) + r*alpha*beta
  P <- (K*a*a_prime*b_prime*m + K*a*b*r*alpha*beta - K*(a^2)*b*m_prime - m*r*alpha*beta)/(alpha*D)
  return(P)
}

Three_species <- function(param, D, R_star, N_star, P_star){

  Jacobian_matrix <- matrix(data = c(0, -
                                     alpha*N_star(param), 
                                   -a_prime*R_star(param), 
                                   beta*alpha*P_star(param), 
                                   0, 
                                   -a*R_star(param), 
                                   b_prime*a_prime*P_star(param), 
                                   b*a*N_star(param),
                                   -r*R_star(param)/K), 3,3)

  ev <- eigen(Jacobian_matrix)
  (values <- ev$values)
  real_part <- c(1,2,3)

  x <- 0
  answer <- F

  for(i in 1:3){
   real_part[i] <- Re(sqrt(as.complex(values[i])))
  
   if(real_part[i] < 0 ){
      x <- x+1
   }
  }

  if(x ==3 ){
    answer <- T
  }

  return(answer)
}

figure_1 <- matrix(0.5, 101, 101)

for(i in 0:100){
  for(j in 0:100){

    a_prime <- i/100
    b_prime <- j/100
    D <- K*a*a_prime*(b*beta - b_prime) + r*alpha*beta
    
    if((a_prime*m*beta + r*alpha*beta - a*m_prime > 0) &
      (m*a_prime*alpha + K*a*a_prime*b*m_prime - K*(a_prime^2)*b_prime*m - K*a_prime*b_prime*r*alpha) > 0 &
      (K*a*b*r*alpha*beta + K*a*a_prime*b_prime*m - K*(a^2)*b*m_prime - m*r*alpha*beta > 0)& 
      R_star(param) > 0 &
      N_star(param) > 0 &
      P_star(param) > 0){ #Three species coexistence 1st circumstance (yellow)
      figure_1[i+1,j+1] <- 1.5
      
    }else if (Three_species(param, D, R_star, N_star, P_star) == T){ #Three species coexistence 2nd circumstance (yellow)
      figure_1[i+1,j+1] <- 1.5
      
    }else if((a*b*(m_prime/(a_prime*b_prime))-alpha*(r/a_prime)*(1-(m_prime/(a_prime*b_prime*K)))-m) < 0){ #IGprey excluded (blue)
      figure_1[i+1,j+1] <- 2.5
      
    }else if((a_prime*b_prime*(m/(a*b))+alpha*beta*(r/a)*(1-(m/(a*b*K)))-m_prime) < 0){ #IGpredator excluded (green)
      figure_1[i+1,j+1] <- 3.5
      
    }else if((a*b*(m_prime/(a_prime*b_prime))-alpha*(r/a_prime)*(1-(m_prime/(a_prime*b_prime*K)))-m) > 0 |
            (a_prime*b_prime*(m/(a*b))+alpha*beta*(r/a)*(1-(m/(a*b*K)))-m_prime) > 0){ #Instability (pink)
      figure_1[i+1,j+1] <- 4.5
  
    }
  }  
} 

#rolf <- c(x1,x2,x3,x4)
#print(rolf)

#Commentaire (DEBUG)
# Vérifie les équations. 
# Si les conditions pour l'instabilité sont en premier, toujours instable

# Donc, il semble que l'ordre de classement de l'état aille un impact. Plus contraindre le modèle.

image (figure_1, col = c("red", "yellow", "blue", "green", "pink", "orange"), breaks = 0:6)
 #contour(figure_1, add = TRUE, lwd = 1)

# 5 equilibria:
# (i): All species at 0 densitiy. If K > 0, unstable
# (ii): Basal ressource at K, but IGpredator and IGprey absent. Unstable if K > m/ab, or K > m'/a'b
# (iii): Resource and IGprey are present, with e.d. of m/ab and (r/a)(1-m/abK). Unstable if a'b'(m/ab)+alphabeta(r/a)(1-m/abK)-m' > 0
# (iv): Ressorce and IGpredator are present, e.d. m'/a'b' and (r/a')(1-m'/a'b'K). 
# (v): all species present, with densities x


((a_prime*m*beta + r*alpha*beta - a*m_prime > 0) &
  (m*a_prime*alpha + K*a*a_prime*b*m_prime - K*(a_prime^2)*b_prime*m - K*a_prime*b_prime*r*alpha) > 0 &
  (K*a*b*r*alpha*beta + K*a*a_prime*b_prime*m - K*(a^2)*b*m_prime - m*r*alpha*beta > 0) & 
  R_star(param) > 0 &
  N_star(param) > 0 &
  P_star(param) > 0)

else if(a*b*(m_prime/(a_prime*b_prime))-alpha*(r/a_prime)*(1-(m_prime/(a_prime*b_prime*K)))-m > 0 |
        a_prime*b_prime*(m/(a*b))+alpha*beta*(r/a)*(1-(m/(a*b*K)))-m_prime > 0){ #Instability (green)
  figure_1[i,j] <- 3.5
  print("3.5")
  
}



  #2. Créer un algorithme qui interprète les résultats et, selon les équations d'équilibre, applique un état au résultat (type d'équilibre ou absence
  # d'équilibre)
  #
  #3. Tracer une ligne aux frontières des états pour le représenter graphiquement
  

#Deuxième set de figures ----
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
N0 <- 1 - P0
Condition_Initiale <- c(P=P0, N=N0)

#paramètres selon la Figure 2
alpha     <- 0.5
b         <- 1
b_prime   <- 1
beta      <- 1
e         <- 1
e_prime   <- 1
I         <- 1
m         <- 0.1
m_prime   <- 1

para <- c(alpha=alpha,b=b,b_prime=b_prime,beta=beta,e=e,
          e_prime=e_prime,I=I,m=m,m_prime=m_prime)

S8_soln <- ode(Condition_Initiale, seq(1,30), S8, para)

P <- c(S8_soln[,'P'])
N <- c(S8_soln[,'N'])


x <- 1:30
y1 <- P
y2 <- N

df <- data.frame(x = x, y = c(y1,y2), group = rep(c("P", "N"), each = length(x)))

ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line() +
  labs(title = "Données", x = "X-axis", y = "Y-axis") +
  theme_minimal()

#add Nullclines

?nullclines

