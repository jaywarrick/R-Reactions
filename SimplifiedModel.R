rm(list=ls())
source('C:/Users/YINLAB/Documents/Jay/SimplifiedModelReactions.R')
library(deSolve)


#######################################
######## Parameters         ###########
#######################################
k_Transcription=5e-4
k_Translation=5e-20
k_Degradation = 1/(60*60*2.5)
params <- list(
     
     # Host Constants
     HOSTT <- 2.5e5,
     RIB_Initial <- 5e6,
     RIB_Spacing <- 238.5,
     RIB_Rate <- 18,
     
     # Viral Constants
     POL_Spacing <- 170,
     POL_Rate <- 3.7,
     MOI <- 1,
     lNT <- 1333,
     lPT <- 822,
     lMT <- 838,
     lGT <- 1672,
     lPOLT <- 6380,
     lGEN_Transcription <- lNT + 0.75*lPT + 0.75*0.75*lMT + 0.75*0.75*0.75*lGT + 0.75*0.75*0.75*0.75*lPOLT, # average length of a genome for a polymerase performing transcription
     lGEN_Replication <- lNT + lPT + lMT + lGT + lPOLT, # average length of a genome/agenome for a polymerase performing replication
     
     #Translation
     k_Translation = k_Translation,
     k_Translate_HOSTT=k_Translation, 
     k_Translate_NT=k_Translation, 
     k_Translate_POLT=k_Translation,
     
     #Transcription
     k_Transcription=k_Transcription,
     
     # AGEN Production
     k_GEN2AGEN = 1e-20*k_Transcription,
     
     # GEN Production
     k_AGEN2GEN = 1e-20*k_Transcription,
     
     # RIB Equilibrium
     k_RIB_EQ = 1e-2,
     
     # POL Equilibrium
     k_POL_EQ = 1e-2,
     
     # Degradation (2.5 hr half life)
     k_Degradation = k_Degradation,
     k_Degrade_NT = k_Degradation,
     k_Degrade_POLT = k_Degradation,
     k_Degrade_HOSTT = k_Degradation
)


#######################################
######## Initial Conditions ###########
#######################################
Xini <- list(
     
     # Proteins
     RIB=RIB_Initial,
     N=0,
     POL=50*MOI,
     HOST=0,
     POL_Tot=50*MOI,
     RIB_Tot=RIB_Initial,
     
     # Transcripts
     NT=0, 
     POLT=0,
     HOSTT=HOSTT,
     
     # Templates
     GEN=MOI,
     AGEN=0
)


#######################################
######## Simplified Model   ###########
#######################################

SimplifiedModel <- ODEModel('SimplifiedModel')

# preCode is evaluated at each timestep and has all Parameters and state variables available to it
preCode <- '

# Calculations for RIB_EQ
k2 = k_Translation
PolTot = RIB_Tot
lHOSTT_sum = HOSTT*lHOSTT
lNT_sum = NT*lNT
lPT_sum = 0.75*NT*lPT
lMT_sum = 0.75*0.75*NT*lMT
lGT_sum = 0.75*0.75*0.75*NT*lGT
lPOLT_sum = POLT*lPOLT
ltot2 = lHOSTT_sum + lNT_sum + lPT_sum + lMT_sum + lGT_sum + lPOLT_sum
HOSTT_frac = lHOSTT_sum/ltot2
NT_frac = lNT_sum/ltot2
POLT_frac = lPOLT_sum/ltot2
PolSpacing = RIB_Spacing
PolRate = RIB_Rate
tsp2 = PolSpacing/PolRate
RIB_ON <- PolTot - (sqrt((k2*PolRate*PolTot*tsp2)^2+(2*k2*PolRate^2-2*k2^2*ltot2*PolRate)*PolTot*tsp2+PolRate^2+2*k2*ltot2*PolRate+(k2*ltot2)^2)+k2*PolRate*PolTot*tsp2-PolRate-k2*ltot2)/(2*k2*PolRate*tsp2)
RIB_UNITS_ON <- RIB_ON/ltot2

# Calculations for POL_EQ
k2 = k_Transcription
PolTot = POL_Tot
lTranscription_sum = lGEN_Transcription
lGEN2AGEN_sum = GEN*lGEN_Replication*lGEN*k_GEN2AGEN/k_Transcription
lAGEN2GEN_sum = AGEN*lGEN_Replication*k_AGEN2GEN/k_Transcription
ltot2 = ltranscription + lGEN2AGEN + lAGEN2GEN
Transcription_frac = lTranscription_sum/ltot2
AGEN2GEN_frac = lAGEN2GEN_sum/ltot2
GEN2AGEN_frac = lGEN2AGEN_sum/ltot2
PolSpacing = POL_Spacing
PolRate = POL_Rate
tsp2 = PolSpacing/PolRate
POL_ON <- PolTot - (sqrt((k2*PolRate*PolTot*tsp2)^2+(2*k2*PolRate^2-2*k2^2*ltot2*PolRate)*PolTot*tsp2+PolRate^2+2*k2*ltot2*PolRate+(k2*ltot2)^2)+k2*PolRate*PolTot*tsp2-PolRate-k2*ltot2)/(2*k2*PolRate*tsp2)
POL_UNITS_ON <- POL_ON/ltot2

'
setPreCalcCode(SimplifiedModel, preCode)

# Translation
addReaction(SimplifiedModel, Translate_HOSTT)
addReaction(SimplifiedModel, Translate_NT)
addReaction(SimplifiedModel, Translate_POLT)

# Transcription
addReaction(SimplifiedModel, Transcription)

# AGEN Production
addReaction(SimplifiedModel, GEN2AGEN)

# GEN Production
addReaction(SimplifiedModel, AGEN2GEN)

# Degradation
#addReaction(SimplifiedModel, Degrade_NT)
#addReaction(SimplifiedModel, Degrade_POLT)
#addReaction(SimplifiedModel, Degrade_HOSTT)

# Equilibrium Kinetics
# Needed temporarily to build up simulation from translation... addDE(SimplifiedModel, dX='dPOL_Tot', eq='0')addDE(SimplifiedModel, dX='dAGEN', eq='0')
# addDE(SimplifiedModel, dX='dN', eq='RIB_Bound*NT_Fraction')
addDE(SimplifiedModel, dX='dN', eq='RIB_ON*NT_frac')
addDE(SimplifiedModel, dX='dRIB', eq='0')
#addDE(SimplifiedModel, dX='dRIB', eq='k_RIB_EQ*(RIB_EQ-RIB)')
addDE(SimplifiedModel, dX='dPOL', eq='k_POL_EQ*(POL_EQ-POL)')




#######################################
######## Run The Simulation ###########
#######################################

# Get the simulation function for the model
VSVFunc <- makeODEFunc(SimplifiedModel)
# Get the initial conditions vector for the simulation function
Xini_OrderedVector <- getOrderedICVector(SimplifiedModel,Xini)
times <- seq(0, 60*60*6, by = 10)
 
dev.off()
out[,'time'] <- out[,'time']/(60*60)
plot(out)
plot(out[,'time'],out[,'HOSTT']/(out[,'HOSTT'] + out[,'NT']), type='l')