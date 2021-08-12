# New Plotting Scripts for Aphid and Infection Responses

# Some updated post lab meeting 2 19 2021

library(tidyverse)
library(ggpubr)
library(plotrix)
library(glmmTMB)

pdiv_clean <- read.csv("../data/Pdiv_cleaned.csv") %>% 
  mutate(proptop = if_else(Adtot == 0, nymtop,(((adtop*Adtot)+(nymtop*Nymtot))/Poptot)), offprop = (1-centprop),cage = paste(Treat,Block,Rep),
         plants = 9)
str(pdiv_clean)
pdiv_aphdata <- pdiv_clean %>% 
  group_by(Treat,Day) %>% 
  filter(!is.na(Adtot), !is.na(Nymtot),!is.na(adtop)) %>% 
  dplyr::summarize(mean_pop = mean(Adtot+Nymtot),popSE = std.error(Adtot+Nymtot),
                   mean_proptop = mean(proptop), proptopSE = std.error(proptop),
                   mean_centprop = mean(centprop), centpropSE = std.error(centprop))



pdiv_clean$Treat = factor(pdiv_clean$Treat, levels = c("Control", "HCHC", "C7C7", "PTPT", "HCC7", "HCPT", "C7PT"))
inf_boxplot <- ggplot(data =pdiv_clean %>% filter(!is.na(totinf))) + geom_boxplot(mapping = aes(x = Treat, y = totinf), size = .75) + 
  theme_classic() + labs(y = "PEMV Infected Hosts") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.position = "none", axis.text.y = element_text(size = 12)) + 
  scale_y_continuous(breaks = 1:9) 
inf_boxplot 


# Figure 2 ####
feeding_pos <- pdiv_clean %>% 
  dplyr::mutate(Top = (Poptot*proptop),
                Bottom = (Poptot - Top)) %>% 
  dplyr::select(Treat,HC,C7,PT,Block,Day,Top,Bottom) %>% 
  pivot_longer(cols = c(Top,Bottom), names_to = "position", values_to = "aphids")
feeding_pos$Treat = factor(feeding_pos$Treat, levels = c("Control", "HCHC", "C7C7", "PTPT", "HCC7", "HCPT", "C7PT"))
feeding_pos$position = factor(feeding_pos$position, levels = c("Top","Bottom"))
head(feeding_pos)
#Viewq(feeding_pos)

feeding_sum <- feeding_pos %>% 
  mutate(Position = if_else(position == "Top", "Upper","Lower")) %>% 
  dplyr::group_by(Treat,Day, Position) %>% 
  dplyr::summarise(mean_aph = mean(aphids), se_aph = std.error(aphids))
#Viewq(feeding_sum)
feeding_sum$Treat = factor(feeding_sum$Treat, levels = c("Control","HCHC", "C7C7", "PTPT", "na","HCC7", "HCPT", "C7PT"))
feeding_sum$Position = factor(feeding_sum$Position, levels = c("Upper","Lower"))

write.csv( feeding_sum,"../data/feeding_sum.csv")
# I cheated using excel 

control_embedded <- read.csv("../data/control_embedded.csv", header = T, sep = ",")
control_embedded <- control_embedded %>% 
  mutate(Position = if_else(position == "Top", "Upper","Lower"))
control_embedded$Treat = factor(control_embedded$Treat, levels = c( "HCHC", "C7C7", "PTPT", "HCC7", "HCPT", "C7PT"))
control_embedded$position = factor(control_embedded$position, levels = c("Upper","Lower"))

aphid_abun_summary <- ggplot(data = feeding_sum, aes(x = Day, y = mean_aph, fill = Position)) +
  geom_bar(stat = "identity", position = 'dodge') + 
  geom_errorbar(aes(ymin = mean_aph - se_aph, ymax = mean_aph + se_aph,color = Position), 
                position = position_dodge(1.8), stat = "identity", width = .4, size = 1) +
  facet_wrap(~Treat, nrow = 2) + theme_classic() + 
  #geom_point(aes(x = Day, y = mean_con, color = Position), size = 2,position = position_dodge(1.8)) + 
  #geom_line(aes(x = Day, y = mean_con, color = Position),position = position_dodge(1.8)) +
  # geom_errorbar(aes(ymin = mean_con - se_con, ymax = mean_con + se_con, color = Position),
         #       position = position_dodge(1.8), stat = "identity", width = .3)+ 
  scale_fill_grey(start = .2, end = .6) + scale_color_grey(start = .2, end = .6, guide = "none")+
  theme(strip.text = element_text(size = 12, face= "bold", hjust = -.0004), legend.position = "right",
        legend.text = element_text(size = 14),legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 14),  axis.title.x = element_text(size = 14), text = element_text(size = 12),
        axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14)) +
  labs(fill = "Location on Plant", y = "Mean Aphid Abundance") + theme(legend.position  = c(1,.175), legend.justification = c(1.25,0))
   
ggsave("../Figures/Figure_2_base.png", aphid_abun_summary, width = 11, height = 7.5) 

# Fig 3 ####

# Incorporating predicted model values
library("effects")
pdmods <- read.csv("../data/Pdiv_cleaned.csv") %>% 
  mutate(proptop = if_else(Adtot == 0, nymtop,(((adtop*Adtot)+(nymtop*Nymtot))/Poptot)),
         cage = paste(Treat,Block,Rep), offprop = (1-centprop))

Day6dat <- pdmods %>% 
  filter(Day == 6)%>%
  mutate(cage = paste(Treat,Block,Rep),
         plants = 9,offprop = (1-centprop))

# Prevalence Moddel (run in GLMM script first) 
inf_bin_glmer.bb.aph <- glmmTMB(propinf ~ Poptot + offprop + proptop + (1|cage), 
                                family = "binomial",
                                weights = plants,
                                data = Day6dat)

#The Interaction
ac.abundance.raw <- effect('Poptot', inf_bin_glmer.bb.aph, se=TRUE, xlevels=100)
#Data Frame
ac.abundance <-as.data.frame(ac.abundance.raw)
#Viewq(ac.abundance)


#Create plot
Plot.abun <-ggplot(data=ac.abundance, aes(x=Poptot, y=fit))+
  geom_line(size=2) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=.2)+
  #ggtitle("Interaction Plot for Weevil Treatment by Nymph Production")+
  theme_bw() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top") +
  xlab("Aphid Abundance") +
  ylab("PEMV Prevalence") +
  geom_point(data=Day6dat, aes(x=Poptot, y=propinf), position=position_jitter(w=0, h=0.05))
#geom_smooth(method='lm',formula=y~x, aes(color=Weevil),level=0.95)
Plot.abun
# plot(ac.interaction.raw, multiline = TRUE)

# Feeding Location
#The Interaction
ac.feeding.raw <- effect('proptop', inf_bin_glmer.bb.aph, se=TRUE, xlevels=100)
#Data Frame
ac.feeding <-as.data.frame(ac.feeding.raw)
#Viewq(ac.feeding)

#Create plot
Plot.feeding <- ggplot(data=ac.feeding, aes(x=proptop*100, y=fit)) +
  geom_line(size=2) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=.2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color = "white"),
        legend.position="top",
        axis.title.y = element_blank()) +
  xlab("% Aphids Feeding High") +
  geom_point(data=Day6dat, aes(x=proptop*100, y=propinf), position=position_jitter(w=0, h=0.05))
Plot.feeding

# dispersal Location
#The Interaction
ac.dispersal.raw <- effect('offprop', inf_bin_glmer.bb.aph, se=TRUE, xlevels=100)
#Data Frame
ac.dispersal <-as.data.frame(ac.dispersal.raw)



#Create plot
Plot.dispersal <- ggplot(data=ac.dispersal, aes(x=offprop*100, y=fit)) +
  geom_line(size=2) +
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se),alpha=.2) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color = "white"),
        legend.position="top",
        axis.title.y = element_blank()) +
  xlab("% Aphids Dispersed") +
  geom_point(data=Day6dat, aes(x=offprop*100, y=propinf), position=position_jitter(w=0, h=0.05))
Plot.dispersal

# Inf Boxplot adjusted
pdiv_clean$Treat = factor(pdiv_clean$Treat, levels = c("Control", "HCHC", "C7C7", "PTPT", "HCC7", "HCPT", "C7PT"))
inf_boxplot <- ggplot(data =pdiv_clean %>% filter(!is.na(propinf))) + geom_boxplot(mapping = aes(x = Treat, y = propinf), size = .75) + 
  theme_classic() + labs(y = "PEMV Prevalence") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 14), axis.text.x = element_text(size = 14),
        legend.position = "none", axis.text.y = element_text(size = 12)) 
inf_boxplot 


new_aph_responses_fig <- ggarrange(inf_boxplot, ggarrange(Plot.abun,Plot.feeding,Plot.dispersal,ncol = 3, common.legend = T, legend = "bottom"), ncol = 1)
ggsave("../Figures/Figure_3_base.png", new_aph_responses_fig, width = 11, height = 7.5) 
