source('C:/Users/YINLAB/Documents/Jay/Reactions.R')


#######################################
######## Translation        ###########
#######################################

############### Reaction ##############
transcript <- 'NT'
protein <- 'N'
rxn <- Reaction(paste('Translate_', transcript, sep=''))
# Reactants
setReactant(rxn, name='RIB', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
# Products
setProduct(rxn, name='RIB', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_NT <- rxn

############### Reaction ##############
transcript <- 'POLT'
protein <- 'POL_Tot'
# Reactants
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='RIB', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
# Products
setProduct(rxn, name='RIB', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_POLT <- rxn

############### Reaction ##############
transcript <- 'HOSTT'
protein <- 'HOST'
# Reactants
rxn <- Reaction(paste('Translate_', transcript, sep=''))
setReactant(rxn, name='RIB', stoich=1, order=1)
setReactant(rxn, name=transcript, stoich=1, order=1)
# Products
setProduct(rxn, name='RIB', stoich=1)
setProduct(rxn, name=transcript, stoich=1)
setProduct(rxn, name=protein, stoich=1)
Translate_HOSTT <- rxn




#######################################
######## Transcription      ###########
#######################################

############### Reaction ##############
rxn <- Reaction('Transcription')
# Reactants
setReactant(rxn, name='POL', stoich=1, order=1)
setReactant(rxn, name='GEN', stoich=1, order=1)
# Products
setProduct(rxn, name='POL', stoich=1)
setProduct(rxn, name='GEN', stoich=1)
setProduct(rxn, name='NT', stoich=1)
setProduct(rxn, name='POLT', stoich=(0.75^4)) # 4 intergenetic attentuations of 75% each
Transcription <- rxn




#######################################
######## AGEN Production    ###########
#######################################

############### Reaction ##############
rxn <- Reaction('GEN2AGEN')
# Reactants
setReactant(rxn, name='POL', stoich=1, order=1)
setReactant(rxn, name='GEN', stoich=1, order=1)
setReactant(rxn, name='N', stoich=1258, order=1)
# Products
setProduct(rxn, name='POL', stoich=1)
setProduct(rxn, name='GEN', stoich=1)
setProduct(rxn, name='AGEN', stoich=1)
GEN2AGEN <- rxn




#######################################
######## GEN Production     ###########
#######################################

############### Reaction ##############
rxn <- Reaction('AGEN2GEN')
# Reactants
setReactant(rxn, name='POL', stoich=1, order=1)
setReactant(rxn, name='AGEN', stoich=1, order=1)
setReactant(rxn, name='N', stoich=1258, order=1)
# Products
setProduct(rxn, name='POL', stoich=1)
setProduct(rxn, name='AGEN', stoich=1)
setProduct(rxn, name='GEN', stoich=1)
AGEN2GEN <- rxn




#######################################
######## Degradation        ###########
#######################################

############### Reaction ##############
rxn <- Reaction('Degrade_NT')
# Reactants
setReactant(rxn, name='NT', stoich=1, order=1)
# Products
setProduct(rxn, name='NT', stoich=0)
Degrade_NT <- rxn

############### Reaction ##############
rxn <- Reaction('Degrade_POLT')
# Reactants
setReactant(rxn, name='POLT', stoich=1, order=1)
# Products
setProduct(rxn, name='POLT', stoich=0)
Degrade_POLT <- rxn

############### Reaction ##############
rxn <- Reaction('Degrade_HOSTT')
# Reactants
setReactant(rxn, name='HOSTT', stoich=1, order=1)
# Products
setProduct(rxn, name='HOSTT', stoich=0)
Degrade_HOSTT <- rxn

