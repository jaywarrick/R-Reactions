rm(list=ls())
source('/Users/jaywarrick/Documents/Yin_Lab/Papers/VSV Modeling/Reactions.R')
library(deSolve)


#######################################
######## Translation        ###########
#######################################
transcript <- 'Nt'
protein <- 'N'
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='Rib', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
setProduct(rxn, name='Rib', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_Nt <- rxn

transcript <- 'Pt'
protein <- 'P'
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='Rib', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
setProduct(rxn, name='Rib', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_Pt <- rxn

transcript <- 'Mt'
protein <- 'M'
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='Rib', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
setProduct(rxn, name='Rib', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_Mt <- rxn

transcript <- 'Gt'
protein <- 'G'
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='Rib', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
setProduct(rxn, name='Rib', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_Gt <- rxn

transcript <- 'Lt'
protein <- 'L'
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='Rib', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
setProduct(rxn, name='Rib', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_Lt <- rxn





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

Xini <- c(
    A=10, 
    B=10, 
    C=0, 
    D=0,
    E=0,
    F=10
)

times <- seq(0, 200, by = 1)
out   <- ode(Xini, times, VSVFunc, params)
summary(out)
dev.off()
plot(out)