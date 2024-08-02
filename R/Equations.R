### Equations
# Volume by L and D
Eq1V <- function(L,D){
  (pi/4)*L*(D^2)
}

# Mass (assumption: constant)
Eq2m <- function(V,gd){
  V*gd
}

# Shape (assumption: constant)
Eq3S <- function(L,D){
  L/D
}

# Volume by S and L
Eq4V <- function(S,L) {
  (pi/(4*(S^2)))*(L^3)
}

Eq5 <- function(){
  
}

L <- 0.2927
D <- 0.0183
gd <- 1.121
rD <- 0.03 #Assumed rate between 0.01 and 0.05
Eq6 <- function(r, t){
  S <- Eq3S(L,D)
  V <- Eq4V(S, L)
  m <- Eq2m(V, gd)
  m * exp(-r*t)
}


Eq7 <- function(){}

Eq8 <- function(){}

L <- 0.2927
D <- 0.0183
gd <- 1.121
fd <- 1.025
vis <- 0.0189
coef <- 0.0790
g <- 981
p <- -1.644

eq9 <- function(L, D, rhoF, rhoW, visc = vis, 
                grav = g, pow = p, const = coef){
  const / visc * {rhoF - rhoW} * grav * L^2 * {L / D}^pow
}

eq10 <- function(L, D, rhoF, rhoW, visc = vis, 
                 grav = g, pow = p, const = coef){
  S <- Eq3S(L,D)
  V <- Eq4V(S, L)
  m <- Eq2m(V, rhoF)
  const / visc * {rhoF - rhoW} * grav * S^pow * 
    {{4 * S^2} / {pi * rhoF}}^{2/3} * m^{2/3}
}

eq10.2 <- function(L, D, m, rhoF, rhoW, visc = vis, 
                   grav = g, pow = p, const = coef){
  S <- Eq3S(L,D)
  const / visc * {rhoF - rhoW} * grav * S^pow * 
    {{4 * S^2} / {pi * rhoF}}^{2/3} * m^{2/3}
}




