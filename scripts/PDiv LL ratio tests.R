# PDiv LL Ratio Tests ####
# Established as separate script on 1/27/2021 to test full factorial of interaction removals

# Packages ####
library(tidyverse)
library(plotrix)
library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(ggpubr)
library(rcompanion)
library(multcomp)
library(MuMIn)
library(glmmTMB)
library(lmtest)
library(corrplot)


# Data Importation #### 
pdmods <- read.csv("../data/Pdiv_cleaned.csv") %>% 
  mutate(proptop = if_else(Adtot == 0, nymtop,(((adtop*Adtot)+(nymtop*Nymtot))/Poptot)),
         cage = paste(Treat,Block,Rep), offprop = (1-centprop))

# Base Models ####
# Aphid Population
pop.int <- glmmTMB(Poptot ~ (HC + C7 + PT + Day)^2 + (1|cage),  
                            family = "nbinom1", 
                            data = pdmods)
pop.base<- glmmTMB(Poptot ~ HC*Day + C7*Day + PT*Day + (1|cage),  
                                  family = "nbinom1", 
                                  data = pdmods)

# Off-plant dispersal
off <- glmer(offprop ~ (HC + C7 + PT + Day)^2 + (1|cage), 
                       family = "binomial", 
                       weights = Poptot, data = pdmods)

# Aphid feeding location
top <- glmer(proptop ~ (HC + C7 + PT + Day)^2 + (1|cage), 
                       family = "binomial", 
                       weights = Poptot, data = pdmods)

# Null hypothesis: There is no interaction between predator treatments

# Population Tests

pop.int <- glmmTMB(Poptot ~ (HC + C7 + PT + Day)^2 + (1|cage),  
                   family = "nbinom1",
                   data = pdmods)

pop.base<- glmmTMB(Poptot ~ HC*Day + C7*Day + PT*Day + (1|cage),  
                   family = "nbinom1", 
                   data = pdmods)

pop.hcc7 <- glmmTMB(Poptot ~ HC*Day + C7*Day + PT*Day + HC*PT + + C7*PT + (1|cage),  
                  family = "nbinom1",                           
                  data = pdmods)
                                                              
pop.hcpt <- glmmTMB(Poptot ~ HC*Day + C7*Day + PT*Day + HC*C7 + C7*PT + (1|cage),  
                  family = "nbinom1", 
                  data = pdmods)
                                                                      
pop.c7pt <- glmmTMB(Poptot ~ HC*Day + C7*Day + PT*Day + HC*C7 + HC*PT + (1|cage),  
                  family = "nbinom1", 
                  data = pdmods)

anova(pop.int,pop.base, type = 'LLR') # p = .015  *
anova(pop.int,pop.hcc7, type = 'LLR') # p = .07   ' 
anova(pop.int,pop.hcpt, type = 'LLR') # p = .64
anova(pop.int,pop.c7pt, type = 'LLR') # p = .0049 **

# Dispersal Tests

off.int <-  glmer(offprop ~ (HC + C7 + PT + Day)^2 + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods)

off.base <- glmer(offprop ~ HC*Day + C7*Day + PT*Day + (1|cage), 
                 family = "binomial", 
                 weights = Poptot, data = pdmods)

off.hcc7 <- glmer(offprop ~ HC*Day + C7*Day + PT*Day + HC*PT + + C7*PT + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods) 

off.hcpt <- glmer(offprop ~ HC*Day + C7*Day + PT*Day + HC*C7 + C7*PT + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods) 

off.c7pt <- glmer(offprop ~ HC*Day + C7*Day + PT*Day + HC*C7 + HC*PT + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods)

anova(off.int,off.base, type = 'LLR') # p = .38
anova(off.int,off.hcc7, type = 'LLR') # p = .41
anova(off.int,off.hcpt, type = 'LLR') # p = .16 eh
anova(off.int,off.c7pt, type = 'LLR') # p = .97


# Feeding Location Tests 

top.int <- glmer(proptop ~ (HC + C7 + PT + Day)^2 + (1|cage), 
                 family = "binomial", 
                 weights = Poptot, data = pdmods)

top.base <- glmer(proptop ~ HC*Day + C7*Day + PT*Day + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods)

top.hcc7 <- glmer(proptop ~ HC*Day + C7*Day + PT*Day + HC*PT + + C7*PT + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods) 

top.hcpt <- glmer(proptop ~ HC*Day + C7*Day + PT*Day + HC*C7 + C7*PT + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods) 

top.c7pt <- glmer(proptop ~ HC*Day + C7*Day + PT*Day + HC*C7 + HC*PT + (1|cage), 
                  family = "binomial", 
                  weights = Poptot, data = pdmods)

anova(top.int,top.base, type = 'LLR') # p = .82
anova(top.int,top.hcc7, type = 'LLR') # p = .47
anova(top.int,top.hcpt, type = 'LLR') # p = .57
anova(top.int,top.c7pt, type = 'LLR') # p = .81








