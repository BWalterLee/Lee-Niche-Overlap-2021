# P Div GLMMs 

# Initially Assembled on 1/4/21
# 

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
library(aod)
library(broom)
install.packages("broom.mixed")
pdmods <- read.csv("../data/Pdiv_cleaned.csv") %>% 
  mutate(proptop = if_else(Adtot == 0, nymtop,(((adtop*Adtot)+(nymtop*Nymtot))/Poptot)),
         cage = paste(Treat,Block,Rep), offprop = (1-centprop))

Day6dat <- pdmods %>% 
  filter(Day == 6)%>%
  mutate(cage = paste(Treat,Block,Rep),
         plants = 9,offprop = (1-centprop))

# Quick autocorr test
aph_corcheck <- cor(Day6dat %>% dplyr::select(HC,C7,PT,Poptot,proptop,offprop,propinf))
corrplot.mixed(aph_corcheck)

# Final Models for Reporting ####

# Aphid Population
poptot.glmer.nb1 <- glmmTMB(Poptot ~ (HC + C7 + PT + Day)^2 + (1|cage),  
                            family = "nbinom1", 
                            data = pdmods)
Anova(poptot.glmer.nb1)
summary(poptot.glmer.nb1)

# Off-plant dispersal
offprop.glmer <- glmer(offprop ~ (HC + C7 + PT + Day)^2 + (1|cage), 
                        family = "binomial", 
                        weights = Poptot, data = pdmods)
summary(offprop.glmer)
Anova(offprop.glmer)


# Aphid feeding location
proptop.glmer <- glmer(proptop ~ (HC + C7 + PT + Day)^2 + (1|cage), 
                       family = "binomial", 
                       weights = Poptot, data = pdmods)
Anova(proptop.glmer)
summary(proptop.glmer)

# Predator Effects on Infection
inf_bin_glmer.b <- glmmTMB(propinf ~ (HC + C7 + PT)^2 + (1|cage), 
                            family = "binomial",
                            weights = plants,
                            data = Day6dat)

summary(inf_bin_glmer.b)
Anova(inf_bin_glmer.b)

# Aphid response Effects on Infection
inf_bin_glmer.bb.aph <- glmmTMB(propinf ~ Poptot + offprop + proptop + (1|cage), 
                                 family = "binomial",
                                 weights = plants,
                                 data = Day6dat)

summary(inf_bin_glmer.bb.aph)
Anova(inf_bin_glmer.bb.aph)

inf_bin_glmer.bb.aph.i <- glmmTMB(propinf ~ (Poptot + offprop + proptop)^2 + (1|cage),  # To report why no int included
                                family = "binomial",
                                weights = plants,
                                data = Day6dat)
summary(inf_bin_glmer.bb.aph.i)
BIC(inf_bin_glmer.bb.aph,inf_bin_glmer.bb.aph.i)


pop.t<- broom.mixed::tidy(poptot.glmer.nb1)
off.t<- broom.mixed::tidy(offprop.glmer)
top.t<- broom.mixed::tidy(proptop.glmer)
inf.t<- broom.mixed::tidy(inf_bin_glmer.bb.aph)
inf.pred.t<- broom.mixed::tidy(inf_bin_glmer.b)
write.csv(pop.t, "../tables/pop_t.csv")
write.csv(off.t, "../tables/off_t.csv")
write.csv(top.t, "../tables/top_t.csv")
write.csv(inf.t, "../tables/inf_t.csv")
write.csv(inf.pred.t, "../tables/inf_pred_t.csv")




