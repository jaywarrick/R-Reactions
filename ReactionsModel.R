source('~/Documents/GitHub/R-Reactions/Reactions.R')
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



# 
# 
# deg <- Reaction('Degradation')
# setReactant(deg, name='C', stoich=1, order=1)
# setProduct(deg, name='D', stoich=1)
# print(deg)
# 
# VSV <- ODEModel('VSV')
# addReaction(VSV, prod)
# addReaction(VSV, deg)
# addDE(VSV, dX='dB', eq='-1*k_Degradation*F')
# addDE(VSV, dX='dE', eq='k_Production')
# addDE(VSV, dX='dF', eq='0')
# VSVFunc <- makeODEFunc(VSV)

RxnRates <- function(t, X, params)
{
    with(as.list(c(X, params)), {
        
        
        dA = -1*kon*A*B + koff*AB
        dB = -1*kon*A*B + koff*AB
        dAB = -koff*AB
        
        return(list(c( dA,dB,dAB )))
    })
}

params <- list(
    kon=1e5, 
    koff=1e-4
)

Xini <- c(
    A=C_Ltot, 
    B=C_in, 
    AB=0
)

times <- seq(0, 200, by = 1)
out   <- ode(Xini, times, RxnRates, params)
summary(out)
# dev.off()
plot(out)
