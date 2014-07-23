rm(list=ls())
source('C:/Users/YINLAB/Documents/Jay/Reactions.R')
library(deSolve)

prod <- Reaction('Production')
setReactant(prod, name='A', stoich=1, order=1)
setReactant(prod, name='B', stoich=1, order=1)
setProduct(prod, name='B', stoich=1)
setProduct(prod, name='C', stoich=1)
makeRates(prod)
print(prod)

deg <- Reaction('Degradation')
setReactant(deg, name='C', stoich=1, order=1)
setProduct(deg, name='D', stoich=1)
print(deg)

VSV <- ODEModel('VSV')
addReaction(VSV, prod)
addReaction(VSV, deg)
addDE(VSV, dX='dB', eq='-1*k_Degradation*F')
addDE(VSV, dX='dE', eq='k_Production')
addDE(VSV, dX='dF', eq='0')
VSVFunc <- makeODEFunc(VSV)

params <- list(
    k_Degradation=1e-6, 
    k_Production=1e-6
)

# Be sure that the initial conditions are in the same order as in the VSVFunc
Xini <- list(

     F=10,
     A=10, 
     B=10, 
     C=0, 
     D=0,
     E=0
   
)
print(Xini)
Xini <- getOrderedICVector(VSV,Xini)
print(Xini)

times <- seq(0, 200, by = 1)
out   <- ode(Xini, times, VSVFunc, params)
summary(out)
plot(out)