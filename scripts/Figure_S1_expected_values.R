# P Div Script for Expected Values to plot for aphid responses


library(tidyverse)
library(plotrix)
library(lme4)
library(car)
library(emmeans)
library(ggpubr)

basedata <- read.csv("../data/Pdiv_cleaned.csv", header = T, sep = ",")

basedata <- basedata %>% 
  mutate(proptop = if_else(Adtot == 0, nymtop,(((adtop*Adtot)+(nymtop*Nymtot))/Poptot)),
         offprop = (1-centprop))
head(basedata)

# 1/25/2021 New Test for values for expected lines
# Species Split time Data now
HCC7data.full <- basedata %>% 
  filter(Treat == "HCHC" | Treat == "C7C7" | Treat == "HCC7") %>% 
  dplyr::select(-addist,-nymdist,-centtot,centprop,-propcol)

HCPTdata.full <- basedata %>% 
  filter(Treat == "HCHC" | Treat == "PTPT" | Treat == "HCPT") %>% 
  dplyr::select(-addist,-nymdist,-centtot,centprop,-propcol)

C7PTdata.full <- basedata %>% 
  filter(Treat == "C7C7" | Treat == "PTPT" | Treat == "C7PT") %>% 
  dplyr::select(-addist,-nymdist,-centtot,centprop,-propcol)
# New Ben mini-function for expected data
exp_function <- function(data,pred12,pred1,pred2,value){
  
  # Calculate single species group averages and halve
  pred_halfsum_2 <- (mean((data %>% filter(Treat == pred1 & Day == 2))[,value])+
                       mean((data %>% filter(Treat == pred2& Day == 2))[,value]))/2
  
  pred_halfsum_4 <- (mean((data %>% filter(Treat == pred1 & Day == 4))[,value])+
                       mean((data %>% filter(Treat == pred2& Day == 4))[,value]))/2
  
  pred_halfsum_6 <- (mean((data %>% filter(Treat == pred1 & Day == 6))[,value])+
                       mean((data %>% filter(Treat == pred2& Day == 6))[,value]))/2
  
  full_frame <- data.frame(c(2,4,6), c(pred_halfsum_2,pred_halfsum_4,pred_halfsum_6))
  names(full_frame) <- c("Day", "response")
  return(full_frame)
 
}

HCC7popday <- exp_function(data = HCC7data.full,
                        pred12 = "HCC7", pred1 = "HCHC", pred2 = "C7C7", value = "Poptot")
HCPTpopday <- exp_function(data = HCPTdata.full, 
                        pred12 = "HCPT", pred1 = "HCHC", pred2 = "PTPT", value = "Poptot")
C7PTpopday <- exp_function(data = C7PTdata.full, 
                        pred12 = "C7PT", pred1 = "C7C7", pred2 = "PTPT", value = "Poptot")
HCC7popday
HCPTpopday
C7PTpopday

HCC7dispday <- exp_function(data = HCC7data.full,
                    pred12 = "HCC7", pred1 = "HCHC", pred2 = "C7C7", value = "offprop")
HCPTdispday <- exp_function(data = HCPTdata.full, 
                      pred12 = "HCPT", pred1 = "HCHC", pred2 = "PTPT", value = "offprop")
C7PTdispday <- exp_function(data = C7PTdata.full, 
                           pred12 = "C7PT", pred1 = "C7C7", pred2 = "PTPT", value = "offprop")
HCC7dispday
HCPTdispday
C7PTdispday

HCC7topday <- exp_function(data = HCC7data.full,
                            pred12 = "HCC7", pred1 = "HCHC", pred2 = "C7C7", value = "proptop")
HCPTtopday <- exp_function(data = HCPTdata.full, 
                            pred12 = "HCPT", pred1 = "HCHC", pred2 = "PTPT", value = "proptop")
C7PTtopday <- exp_function(data = C7PTdata.full, 
                            pred12 = "C7PT", pred1 = "C7C7", pred2 = "PTPT", value = "proptop")
HCC7topday
HCPTtopday
C7PTtopday

# Working on multi-panel expected vs observed results

pdiv_aphdata <- basedata %>% 
  group_by(Treat,Day) %>% 
  filter(!is.na(Adtot), !is.na(Nymtot),!is.na(adtop)) %>% 
  dplyr::summarize(mean_pop = mean(Adtot+Nymtot),popSE = std.error(Adtot+Nymtot),
                   mean_top = mean(proptop), topSE = std.error(proptop),
                   mean_disp = mean(offprop), dispSE = std.error(offprop))

pd = position_dodge(.1)
pdiv_aphdata$Treat = factor(pdiv_aphdata$Treat, levels = c("Control", "HCHC", "C7C7", "PTPT", "HCC7", "HCPT", "C7PT"))

# Pop
# If visualizing points is ever wanted (I think it's too busy)
# geom_point(data = HCC7data.full %>% filter(Treat == "HCC7"), aes(x = Day, y = Poptot)) +

HCC7pop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCC7"), aes(x=Day, y=mean_pop)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) + 
  geom_text(label = "*", aes(x = 4, y = 200), size = 12)+ # indicating significance according to LRT
  geom_errorbar(aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response), data = HCC7popday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,300) + scale_x_continuous(breaks = c(2,4,6))+labs(y = "Aphid Population") 
HCC7pop_plot

HCPTpop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCPT"), aes(x=Day, y=mean_pop)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response), data = HCPTpopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,300) + scale_x_continuous(breaks = c(2,4,6))+theme(axis.text.y = element_blank(),axis.title.y = element_blank())
HCPTpop_plot

C7PTpop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "C7PT"), aes(x=Day, y=mean_pop)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_text(label = "*", aes(x = 4, y = 200), size = 12)+ # indicating significance according to LRT
  geom_errorbar(aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response), data = C7PTpopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,300) + scale_x_continuous(breaks = c(2,4,6))+theme(axis.text.y = element_blank(),axis.title.y = element_blank())
C7PTpop_plot


# Disp
HCC7disp_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCC7"), aes(x=Day, y=mean_disp*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCC7dispday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,100) + scale_x_continuous(breaks = c(2,4,6)) +labs(y = "% Aphid Dispersed")
HCC7disp_plot

HCPTdisp_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCPT"), aes(x=Day, y=mean_disp*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCPTdispday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())
HCPTdisp_plot

C7PTdisp_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "C7PT"), aes(x=Day, y=mean_disp*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = C7PTdispday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())
C7PTdisp_plot


# top
HCC7top_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCC7"), aes(x=Day, y=mean_top*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCC7topday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(50,100) + scale_x_continuous(breaks = c(2,4,6)) +labs(y = "% Aphids Feeding High")
HCC7top_plot

HCPTtop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCPT"), aes(x=Day, y=mean_top*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCPTtopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(50,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())
HCPTtop_plot

C7PTtop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "C7PT"), aes(x=Day, y=mean_top*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = C7PTtopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(50,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())
C7PTtop_plot

obs_exp_figure <- ggarrange(
HCC7pop_plot,HCPTpop_plot,C7PTpop_plot,
HCC7disp_plot,HCPTdisp_plot,C7PTdisp_plot,
HCC7top_plot,HCPTtop_plot,C7PTtop_plot,
ncol = 3, nrow = 3, legend = "none", common.legend = T, align = "v",
labels = c("HCC7","HCPT", "C7PT","HCC7","HCPT", "C7PT","HCC7","HCPT", "C7PT"),
hjust = -1.5, vjust = 11, font.label = list(size = 12))

ggsave("../Figures/Figure_S1.png", obs_exp_figure, width = 9, height = 6) 



# With Control Line Added 
# RUN THESE TWO SEPARATELY!!!!

# Pop
# If visualizing points is ever wanted (I think it's too busy)
# geom_point(data = HCC7data.full %>% filter(Treat == "HCC7"), aes(x = Day, y = Poptot)) +

HCC7pop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCC7"), aes(x=Day, y=mean_pop)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) + 
  geom_errorbar(aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response), data = HCC7popday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,400) + scale_x_continuous(breaks = c(2,4,6))+labs(y = "Aphid Population") +
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_pop ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_pop),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd,color = "grey60")


HCC7pop_plot

HCPTpop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCPT"), aes(x=Day, y=mean_pop)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response), data = HCPTpopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,400) + scale_x_continuous(breaks = c(2,4,6))+theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_pop ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_pop),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd,color = "grey60")

HCPTpop_plot

C7PTpop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "C7PT"), aes(x=Day, y=mean_pop)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response), data = C7PTpopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,400) + scale_x_continuous(breaks = c(2,4,6))+theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_pop ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_pop),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_pop-popSE, ymax=mean_pop+popSE), width=.15, position = pd,color = "grey60")

C7PTpop_plot


# Disp
HCC7disp_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCC7"), aes(x=Day, y=mean_disp*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCC7dispday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,100) + scale_x_continuous(breaks = c(2,4,6)) +labs(y = "% Aphid Dispersed") +
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_disp*100 ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_disp*100),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd,color = "grey60")
HCC7disp_plot

HCPTdisp_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCPT"), aes(x=Day, y=mean_disp*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCPTdispday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_disp*100 ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_disp*100),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd,color = "grey60")
HCPTdisp_plot

C7PTdisp_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "C7PT"), aes(x=Day, y=mean_disp*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = C7PTdispday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(0,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_disp*100 ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_disp*100),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_disp*100-dispSE*100, ymax=mean_disp*100+dispSE*100), width=.15, position = pd,color = "grey60")
C7PTdisp_plot


# top
HCC7top_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCC7"), aes(x=Day, y=mean_top*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCC7topday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(50,100) + scale_x_continuous(breaks = c(2,4,6)) +labs(y = "% Aphids Feeding High")+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_top*100 ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_top*100),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd,color = "grey60")
HCC7top_plot

HCPTtop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "HCPT"), aes(x=Day, y=mean_top*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = HCPTtopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(50,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_top*100 ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_top*100),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd,color = "grey60")
HCPTtop_plot

C7PTtop_plot <- ggplot(pdiv_aphdata %>%  filter(Treat == "C7PT"), aes(x=Day, y=mean_top*100)) +            
  geom_line(size = 1, position=pd) +
  geom_point(size = 2, position=pd, aes(shape = Treat)) +  
  geom_errorbar(aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd) +
  geom_smooth(mapping = aes(x = Day, y = response*100), data = C7PTtopday, formula = y ~x, size = 1, color = "grey40", linetype = "dashed", method = "glm",se = F) +
  theme_classic() + ylim(50,100) + scale_x_continuous(breaks = c(2,4,6)) +theme(axis.text.y = element_blank(),axis.title.y = element_blank())+
  geom_line(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(x=Day, y = mean_top*100 ),color = "grey60")+
  geom_point(data = pdiv_aphdata %>% filter(Treat == "Control"), size = 2, position=pd, aes(x=Day, y = mean_top*100),color = "grey60") + 
  geom_errorbar(data = pdiv_aphdata %>% filter(Treat == "Control"), aes(ymin=mean_top*100-topSE*100, ymax=mean_top*100+topSE*100), width=.15, position = pd,color = "grey60")
C7PTtop_plot

obs_exp_figure_controls <- ggarrange(
  HCC7pop_plot,HCPTpop_plot,C7PTpop_plot,
  HCC7disp_plot,HCPTdisp_plot,C7PTdisp_plot,
  HCC7top_plot,HCPTtop_plot,C7PTtop_plot,
  ncol = 3, nrow = 3, legend = "none", common.legend = T, align = "v",
  labels = c("HC:C7","HC:PT", "C7:PT","HC:C7","HC:PT", "C7:PT","HC:C7","HC:PT", "C7:PT"),
  hjust = -1.5, vjust = 11.5, font.label = list(size = 12))

ggsave("../obs_exp_figure_controls.png", obs_exp_figure_controls, width = 9, height = 6) 





