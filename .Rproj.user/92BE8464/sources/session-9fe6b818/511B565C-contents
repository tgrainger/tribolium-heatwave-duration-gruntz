library(ggplot2)
library(tidyverse)
library(broom)
library(purrr)
library(DHARMa)
library(cowplot)
library(ggpubr)
library(car)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(patchwork)

getwd()

#my ggplot theme
theme_tess <- function () { 
  theme_cowplot()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))+
    theme(axis.text.x=element_text(size=22))+
    theme(axis.text.y=element_text(size=22))+
    theme(axis.title.x=element_text(size=22))+
    theme(axis.title.y=element_text(size=22))+
    theme(plot.title = element_text(hjust = 0.5,size=20))}

theme_tess_small <- function () { 
  theme_cowplot()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
    theme(axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))+
    theme(axis.text.x=element_text(size=10))+
    theme(axis.text.y=element_text(size=10))+
    theme(axis.title.x=element_text(size=10))+
    theme(axis.title.y=element_text(size=10))+
    theme(plot.title = element_text(hjust = 0.5,size=10))+
    theme(axis.title.y=element_text(size=10))}

# load and manipulate population size data
data <- read_csv("./population counts.csv") %>%
  mutate(ID = paste(treatment, replicate, sep = "_")) %>%
  group_by(ID) %>%
  mutate(live_lag1 = lag(live, n = 1, order_by = count),
         ftreatment = factor(treatment, levels = c("control", "5 day", "10 day", "15 day")),
         treatment = factor(treatment, levels = c("control", "5 day", "10 day", "15 day")))%>%
  mutate(deathrate=dead/live,
         birthrate=births/live_lag1)
         
#get fecundity data in
fec <-read.csv("./fecundity.csv",stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","") )

#order fecundity treatments
fec$treatment <- factor(fec$treatment, 
                        levels = c("control", "5 day", "10 day", "15 day"))

fec<-fec[complete.cases(fec),] #remove rows with NA for eggs

#### FIGURE 1 ####

## prep for figure

#get predicted values to add lines for gompertz model

gompertz_data <- data %>%
  arrange(treatment, ID, count) %>%
  group_by(treatment, ID) %>%
  mutate(live_lag1 = lag(live)) %>%
  ungroup() %>%
  filter(!is.na(live_lag1),
         live > 0,
         live_lag1 > 0) #remove ones with zeros

#fit gompertz model (linearize Ives et al. 2003 one)
gompertz_lm <- lm(log(live / live_lag1) ~ log(live_lag1) * treatment,
  data = gompertz_data)

# get coefficients (intercept (a) and slope (b))
coefs <- coef(gompertz_lm)

# a and b estimmates for each treatment
gompertz_fits <- gompertz_data %>%
  distinct(treatment) %>%
  mutate(
    a = coefs["(Intercept)"] +
      if_else(paste0("treatment", treatment) %in% names(coefs),
        coefs[paste0("treatment", treatment)],0),
    b = coefs["log(live_lag1)"] +
      if_else(paste0("log(live_lag1):treatment", treatment) %in% names(coefs),
        coefs[paste0("log(live_lag1):treatment", treatment)],0))

#get starting N value for each treatment (N at count 1)
start_N <- data %>%
  filter(count == 1) %>%
  group_by(treatment) %>%
  summarise(N0 = mean(live, na.rm = TRUE), .groups = "drop")

#offsets to line up start and end of lines with points on plot
offsets <- c("15 day"  = -0.225,
            "10 day"  = -0.075,
            "5 day"   =  0.075,
            "control" =  0.225)

# get Gompertz predictions
gompertz_pred <- gompertz_fits %>%
  left_join(start_N, by = "treatment") %>%
  mutate(
    pred = pmap(
      list(a, b, N0),
      function(a, b, N0) {
        t <- seq(1, 5, length.out = 200)
        logN0 <- log(N0)
        logN <- a * (1 - (b + 1)^(t - 1)) / (1 - (b + 1)) +((b + 1)^(t - 1)) * logN0 #closed form linearized gompertz solution
        N <- exp(logN)
        tibble(count = t,fit = N)})) %>%
  select(treatment, pred) %>%
  unnest(pred) %>%
  mutate(x = count + offsets[treatment])

### PLOT

##Population dynamics (panel a)

live<-ggplot(data, aes(x = count, y = live, color=treatment)) +
  geom_point(position = position_dodge(width = 0.6),
             size = 3, alpha=0.8,shape=16) + #mean points
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heat wave\nduration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  expand_limits(y = 0)+
  xlab("") +
  ylab("Population size")+
  scale_x_discrete(limits = c("1", "2", "3","4","5"),
                   labels=c("0", "1", "3","5","7"))+
  scale_y_continuous(limits = c(0, 530), breaks = seq(0, 500, by = 100))+
  theme_tess()+
  theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16))+
  geom_line(data = gompertz_pred,
  aes(x = x, y = fit, color = treatment, group = treatment),
  linewidth = 0.9, alpha=0.7,
  inherit.aes = FALSE)
  #theme(plot.margin = margin(15, 5, 5, 5))

legend <- get_legend(live)
live  <- live   + theme(legend.position = "none")

#windows();live

ggsave(file="Figures/live beetles.pdf", live, 
       width = 30, height = 21, units = "cm")

## Birth rate (panel b)

#calculate means and standard errors
summary_data_births<-data%>%
  group_by(treatment,count) %>%
  summarise(mean = mean(birthrate,na.rm=TRUE),
            n=sum(!is.na(birthrate)),
            sd = sd(birthrate,na.rm=TRUE),
            se = sd / sqrt(n))

sig_df <- data.frame(
  count = c("1","1","1","1",
            "2","2","2","2",
            "3","3","3","3",
            "4","4","4","4",
            "5","5","5","5"),
  treatment = c("control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day"),
  label = c("A","A","A","A",
            "B","B","A","A", #order was weird - had to fiddle
            "B","AB","A","AB",
            "A","A","A","A", 
            "A","A","A","A"))

sig_df<-sig_df%>%
  mutate(count=as.integer(count),
          y_fixed = ifelse(row_number() <= 4, 11,
            ifelse(row_number()<= 8,8.7,0.6)))

births_main<-ggplot(summary_data_births, aes(x = count, y = mean, color=treatment)) +
  geom_point(size = 4.5,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = data, aes(x=count,y = birthrate, color = treatment), #all points
             position = position_dodge(width = 0.6), size = 2.5, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heat wave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(override.aes = list(color = "white")))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0,
                linewidth = 1.2) +
  xlab("") +
  ylab("Birth rate")+
  scale_x_discrete(limits = c("1", "2", "3","4","5"),
                   labels=c("0", "1", "3","5","7"))+
  scale_y_continuous(limits=c(0,9.5))+
  geom_text(data = sig_df,
            aes(x = count, y = y_fixed, label = label, group = treatment),
            position = position_dodge(width = 0.6),
            show.legend = FALSE,
            size = 4.5,
            vjust = 0,
            color="#d4791e") +
  theme_tess()+
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16))

births_main<-births_main+theme(legend.position = "none")

#fecundity inset

summary_fec<-fec%>%
  group_by(treatment,monthssinceheatwave) %>%
  summarise(mean = mean(eggs),
            n=n(),
            sd = sd(eggs),
            se = sd / sqrt(n))

summary_fec_pops<-fec%>%
  group_by(treatment,monthssinceheatwave, replicate) %>%
  summarise(mean = mean(eggs),
            n=n(),
            sd = sd(eggs),
            se = sd / sqrt(n))

three<-fec%>%
  filter(monthssinceheatwave==3)
threesum<-summary_fec%>%
  filter(monthssinceheatwave==3)
threesumpops<-summary_fec_pops%>%
  filter(monthssinceheatwave==3)

sig_df <- data.frame(
  treatment = c("control","5 day", "10 day", "15 day"),
  label = c("A","A","AB","B"))%>%
  mutate(y_fixed = c(9.6, 9.6, 9.6, 9.6))

fecundityone<-ggplot(threesum, aes(x = treatment, y = mean, color=treatment)) +
  geom_point(size = 2.5,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = threesumpops, aes(x=treatment,y = mean, color = treatment), #all points
             position = position_dodge(width = 0.6), size = 1.5, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heatwave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0,
                linewidth = 1.2) +
  scale_x_discrete(labels=c("15 days", "10 days", "5 days", "Control"),
                   breaks = c("15 day","10 day","5 day","control"))+
  geom_text(data = sig_df,
            aes(x = treatment, y = y_fixed, label = label, group = treatment),
            position = position_dodge(width = 0.6),
            show.legend = FALSE,
            size = 3.5,
            vjust = 0,
            color="#d4791e") +
  xlab("Heatwave duration") +
  ylab("Fecundity 0 months post heatwave (eggs/48h)")+
  theme_tess_small()+
  theme(legend.position = "none")

#windows();fecundityone

births<-births_main + inset_element(fecundityone,
                left = 0.5,
                bottom = 0.4,
                right = 0.98,
                top = 0.98)
  
#windows();births

ggsave(file="Figures/births.pdf", births, 
       width = 30, height = 21, units = "cm")

## Death rate (panel c)

#calculate means and standard errors
summary_data_dead<-data%>%
  group_by(treatment,count) %>%
  summarise(mean = mean(deathrate,na.rm=TRUE),
            n=sum(!is.na(deathrate)),
            sd = sd(deathrate,na.rm=TRUE),
            se = sd / sqrt(n))

sig_df <- data.frame(
  count = c("1","1","1","1",
            "2","2","2","2",
            "3","3","3","3",
            "4","4","4","4",
            "5","5","5","5"),
  treatment = c("control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day",
                "control","5 day", "10 day", "15 day"),
  label = c("A","A","A","A",
            "A","A","A","A",
            "A","A","A","A",
            "B","AB","A","AB", #the ordering is messed up here for some reason
            "C", "B", "A", "AB")) #the ordering is messed up here for some reason

sig_df$count<-as.integer(sig_df$count)

#where to put significance letters
sig_df <- sig_df %>%
  mutate(y_fixed = ifelse(row_number() <= 4, 0.04,
                          ifelse(row_number() <= 8, 0.04,
                            ifelse(row_number()<=12,0.04,
                                 ifelse(row_number()<=16,0.15,
                                        ifelse(row_number()<=20,0.5))))))

#plot
dead<-ggplot(summary_data_dead, aes(x = count, y = mean, color=treatment)) +
  geom_point(size = 4,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = data, aes(x=count,y = deathrate, color = treatment), #all points
             position = position_dodge(width = 0.6), size = 3, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heatwave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(override.aes = list(color = "white")))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0,
                linewidth = 1.2) +
  xlab("Months since heatwave") +
  ylab("Mortality rate")+
  scale_x_discrete(limits = c("1", "2", "3","4","5"),
                   labels=c("0", "1", "3","5","7"))+
  # scale_y_continuous(limits=c(0,140),
  #                    breaks = seq(0, 140, by = 30))+
  geom_text(data = sig_df,
            aes(x = count, y = y_fixed, label = label, group = treatment),
            position = position_dodge(width = 0.7),
            show.legend = FALSE,
            size = 4.5,
            vjust = 0,
            color="#d4791e") +
  theme_tess()+
  theme(legend.position = "none")

#windows();dead

ggsave(file="Figures/dead beetles.pdf", 
       dead, width = 30, height = 21, units = "cm")

##PUT ALL PLOTS TOGETHER

livedead<-plot_grid(live, legend, births, NULL,dead,NULL,
                    labels = c("A", "","B","", "C",""),
                    label_size = 20,
                    label_fontface = "bold",
                    nrow=3,ncol=2,
                    rel_heights=c(1,1,1), 
                    align="v",
                    rel_widths=c(1,0.15))

ggsave(file="Figures/Figure 1.pdf",
       livedead, width = 28, height = 48, units = "cm")

#### POPULATION SIZE ANALYSIS ####

# gompertz model
gompertz_lm <- lm(log(live/live_lag1) ~ log(live_lag1)*treatment, data = data)

# logistic model
logistic_lm <- lm(log(live/live_lag1) ~ live_lag1*treatment, data = data)

AIC(gompertz_lm,logistic_lm) #gompertz better

#Get Ks for each replicate

# run model for each cage
models <- data %>%
  group_by(ID, treatment, replicate) %>% 
  nest() %>%
  mutate(model = map(data, ~ lm(log(live) ~ log(live_lag1), data = .x)),
         tidy_model = map(model, tidy, conf.int = T)) %>%
  unnest(tidy_model) %>%
  rename(coefficient = term) 

# get parameter estimates for each cage
params <- models %>%
  select(treatment, replicate, ID, coefficient, estimate) %>%
  pivot_wider(names_from = coefficient, values_from = estimate) %>%
  rename(a = `(Intercept)`, b = `log(live_lag1)`) %>% #intercept is the intrinsic growth rate and slope is the strength of density-dependence
  mutate(K = a/(1-b)) # calculate carry capacity (from equation 3 in Ives et al. 2003), which is already on the log scale

# make a new dataframe
write_csv(params, file = "params.csv")

#Analyze Ks

# load and manipulate data
params <- read_csv("params.csv") %>%
  mutate(ID = paste(treatment, replicate, sep = "_")) %>%
  group_by(ID) %>%
  mutate(ftreatment = factor(treatment, levels = c("control", "5 day", "10 day", "15 day")))%>%
  mutate(unloggedK=exp(K),
         unloggeda=exp(a))#add columns with the unlogged K and a estimates

#boxplot of effect of treatment on K
# params %>%
#   ggplot(aes(x = ftreatment, y = unloggedK)) +
#   geom_boxplot() 

# linear model for K
lm_K <- lm(unloggedK ~ ftreatment, data = params)
Anova(lm_K, type=2) #strong treatment effect

emmeans(lm_K, pairwise ~ ftreatment, adjust = "tukey")
#everything is different from one another except control vs 5 day and 10 vs 15 day (same result as before) (A A B B)

#boxplot of effect of treatment on growth rate
# params %>%
#   ggplot(aes(x = ftreatment, y = unloggeda)) +
#   geom_boxplot() 

# linear model for a
lm_a <- lm(a ~ ftreatment, data = params)
Anova(lm_a, type=2) #treatment effect

emmeans(lm_a, pairwise ~ ftreatment, adjust = "tukey")
#only control vs. 15 day are different

### DEATH RATE ANALYSIS ####

#month zero - cannot analyze, no deaths

#month 1 - very few dead
monthone <- data%>%
  filter(count=="2")

lm<-lm(log(deathrate+0.001)~treatment,data=monthone)
Anova(lm, type=2) 
#no treatment effect

#check assumptions
#plot(lm) #looks ok

#month 3

monththree <- data%>%
  filter(count=="3")

lm<-lm(log(deathrate+0.001)~treatment,data=monththree)
Anova(lm, type=2) 
#no treatment effect

#check assumptions
#plot(lm) #looks OK

#month 5

monthfive <- data%>%
  filter(count=="4")

lm<-lm(log(deathrate+0.001)~treatment,data=monthfive)
Anova(lm, type=2) 
#significant effect

#check assumptions
#plot(lm) #looks OK

anova_model<-aov(log(deathrate+0.001)~treatment,data=monthfive) 

tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
#only 15 day vs. control different

#month 7

monthseven <- data%>%
  filter(count=="5")

lm<-lm(log(deathrate+0.001)~treatment,data=monthseven)
Anova(lm, type=2) 
#treatment effect

#check assumptions
#plot(lm) #looks ok

anova_model<-aov(log(deathrate+0.001)~treatment,data=monthseven)

tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
#5 day and control the same, 5 and 10 day the same, 10 and 15 day the same (A A AB B)



#### BIRTH RATE ANALYSIS ####

#month zero - not applicable (no data)

#month 1 
monthone <- data%>%
  filter(count=="2")

lm<-lm(log(birthrate+0.001)~treatment,data=monthone)
Anova(lm, type=2) 
#strong treatment effect

#check assumptions
#plot(lm) #looks fine (18 is an outlier but it's real)

anova_model<-aov(log(birthrate+0.001)~treatment,data=monthone) 

tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
#everything different except 5 day same as control and 10 day same as 15 day 

#month 3
monththree <- data%>%
  filter(count=="3")

lm<-lm(log(birthrate+0.001)~treatment,data=monththree)
Anova(lm, type=2) 
#treatment effect

#check assumptions
#plot(lm) #looks OK (some outliers caused by 0 births)

#try removing 3 zero value "outliers"
monththree_trimmed<-data%>%
  filter(count=="3",births!="0")

lm<-lm(log(birthrate+0.001)~treatment,data=monththree_trimmed)
Anova(lm, type=2) #still significant - removing outliers doesnt change result

#plot(lm) #looks better

#concluded that the zeros are real (not experimental error) and kept them in

anova_model<-aov(log(birthrate+0.001)~treatment,data=monththree) 

tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
#15 day different from control (note: more births in 15 day than in control)

#month 5

monthfive <- data%>%
  filter(count=="4")

lm<-lm(log(birthrate+0.001)~treatment,data=monthfive)
Anova(lm, type=2) 
#no treatment effect

#check assumptions
#plot(lm) #looks OK

#month 7

monthseven <- data%>%
  filter(count=="5")

lm<-lm(log(birthrate+0.001)~treatment,data=monthseven)
Anova(lm, type=2) 
#no treatment effect

#check assumptions
#plot(lm) #looks OK


#### FECUNDITY ####

##Plot

#get means

summary_fec<-fec%>%
  group_by(treatment,monthssinceheatwave) %>%
  summarise(mean = mean(eggs),
            n=n(),
            sd = sd(eggs),
            se = sd / sqrt(n))

summary_fec_pops<-fec%>%
  group_by(treatment,monthssinceheatwave, replicate) %>%
  summarise(mean = mean(eggs),
            n=n(),
            sd = sd(eggs),
            se = sd / sqrt(n))

zero<-fec%>%
  filter(monthssinceheatwave==0)
zerosum<-summary_fec%>%
  filter(monthssinceheatwave==0)
zerosumpops<-summary_fec_pops%>%
  filter(monthssinceheatwave==0)

one<-fec%>%
  filter(monthssinceheatwave==1)
onesum<-summary_fec%>%
  filter(monthssinceheatwave==1)
onesumpops<-summary_fec_pops%>%
  filter(monthssinceheatwave==1)

three<-fec%>%
  filter(monthssinceheatwave==3)
threesum<-summary_fec%>%
  filter(monthssinceheatwave==3)
threesumpops<-summary_fec_pops%>%
  filter(monthssinceheatwave==3)

feczero<-ggplot(zerosum, aes(x = treatment, y = mean, color=treatment)) +
  geom_point(size = 4,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = zerosumpops, aes(x=treatment,y = mean, color = treatment), #all points
             position = position_dodge(width = 0.6), size = 3, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heatwave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0) +
  #scale_y_continuous(limits=c(0,5))+
  scale_x_discrete(labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"))+
  xlab("Heatwave duration") +
  ylab("Fecundity (eggs laid in 48h)")+
  labs(title = "Zero months post heatwave")+
  theme_tess()+
  theme(legend.position = "none")

#windows();feczero

fecone<-ggplot(onesum, aes(x = treatment, y = mean, color=treatment)) +
  geom_point(size = 4,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = onesumpops, aes(x=treatment,y = mean, color = treatment), #all points
    position = position_dodge(width = 0.6), size = 3, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heatwave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0) +
  #scale_y_continuous(limits=c(0,5))+
  scale_x_discrete(labels=c("15 days", "10 days", "5 days", "Control"),
                   breaks = c("15 day","10 day","5 day","control"))+
  xlab("Heatwave duration") +
  ylab("Fecundity (eggs laid in 48h)")+
  labs(title = "One month post heatwave")+
  theme_tess()+
  theme(legend.position = "none")

#windows();fecone

fecthree<-ggplot(threesum, aes(x = treatment, y = mean, color=treatment)) +
  geom_point(size = 4,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = threesumpops, aes(x=treatment,y = mean, color = treatment), #all points
    position = position_dodge(width = 0.6), size = 3, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heatwave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0) +
  #scale_y_continuous(limits=c(0,10))+
  scale_x_discrete(labels=c("15 days", "10 days", "5 days", "Control"),
                   breaks = c("15 day","10 day","5 day","control"))+
  xlab("Heatwave duration") +
  ylab("Fecundity (eggs laid in 48h)")+
  labs(title = "Three months post heatwave")+
  theme_tess()+
  theme(legend.position = "none")

#windows();fecthree

##PUT PLOTS TOGETHER for supp mat figure

fecundity<-plot_grid(feczero,NULL, fecthree,
                     align="v",ncol=3,rel_widths=c(1,0.1,1))

ggsave(file="Figures/fecundity_zero_three.pdf",
       fecundity, width = 40, height = 20, units = "cm")

##Analysis

#month zero

monthzero <- fec%>%
  filter(monthssinceheatwave=="0")

hist(monthzero$eggs) #not normal (more poisson)

m1<-glmer(eggs ~ treatment + (1 | replicate),
          data = monthzero,
          family = poisson)
Anova(m1, type=2) #treatment effect

#check assumptions
sim <- simulateResiduals(m1)
plot(sim) #fine

emm <- emmeans(m1, ~ treatment)
pairs(emm, adjust = "tukey")
#control and 5 day are both different from 15 day


#month 1

monthone <- fec%>%
  filter(monthssinceheatwave=="1")

hist(monthone$eggs) #not normal

#poisson
m1<-glmer(eggs ~ treatment + (1 | replicate),
          data = monthone,
          family = poisson)

#check assumptions
sim <- simulateResiduals(m1)
plot(sim) #not good- overdispersed and zero inflated

#try zero inflated negative binomial distribution
m2 <- glmmTMB(eggs ~ treatment + (1 | replicate),
              ziformula = ~1,
            data = monthone,
            family = nbinom2)
Anova(m2, type=2) #no effect of treatment

#check assumptions
sim <- simulateResiduals(m2)
plot(sim) #all fine

#month 3

monththree <- fec%>%
  filter(monthssinceheatwave=="3")

m1<-glmer(eggs ~ treatment + (1 | replicate),
          data = monththree,
          family = poisson)

#check assumptions
sim <- simulateResiduals(m1)
plot(sim) #not good- overdispersed and zero inflated!

#try negative binomial
m2 <- glmmTMB(eggs ~ treatment + (1 | replicate),
              ziformula = ~1,
              data = monththree,
              family = nbinom2)
Anova(m2, type=2) #no effect of treatment

#check assumptions
sim <- simulateResiduals(m2)
plot(sim) #all fine


#### BODY SIZE ####

#get data in
body <-read.csv("body size.csv",stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","") )

#orders the treatments in the order that we want them
body$treatment <- factor(body$treatment, 
                        levels = c("control", "5 day", "10 day", "15 day"))

body<-body[complete.cases(body),] #remove rows with NA for weight

#make a column for weight in milligrams (for nicer axis labels)
body$weightmicro<-body$weight*1000

body$sex <- factor(body$sex,
  levels = c("m", "f"),
  labels = c("Female", "Male"))

summary_body_pops<-body%>%
  group_by(sex,treatment,monthssinceheatwave, replicate) %>%
  summarise(meanweightmicro = mean(weightmicro),
            n=n(),
            sd = sd(weightmicro),
            se = sd / sqrt(n))

summary_body<-summary_body_pops%>%
  group_by(sex,treatment,monthssinceheatwave) %>%
  summarise(mean = mean(meanweightmicro),
            n=n(),
            sd = sd(meanweightmicro),
            se = sd / sqrt(n))

threesumpops<-summary_body_pops%>%
  filter(monthssinceheatwave==3)
threesum<-summary_body%>%
  filter(monthssinceheatwave==3)

#plot

bodythree<-ggplot(threesum, aes(x = treatment, y = mean, color=treatment)) +
  geom_point(size = 4,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = threesumpops, aes(x=treatment,y = meanweightmicro, color = treatment), #all points
             position = position_dodge(width = 0.6), size = 3, alpha=0.5,shape=16) +
  facet_wrap(~ sex, strip.position = "bottom") +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heat wave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0) +
  scale_x_discrete(labels = c(
      "control" = "Control",
      "5 day"  = "5 day",
      "10 day" = "10 day",
      "15 day" = "15 day"))+
  #scale_y_continuous(limits=c(0,5))+
  ylab("Body size (mg)")+
  labs(title = "Three months post heatwave")+
  theme_tess()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(strip.placement = "outside",
    strip.background = element_blank(),
    axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(
    panel.spacing = unit(2, "lines"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20)
  )

#windows();bodythree

ggsave(file="Figures/body size.pdf",
       bodythree, width = 40, height = 24, units = "cm")

#Analysis

m3<-subset(threesumpops, sex=="Male")
f3<-subset(threesumpops,sex=="Female")

#Month 3 males
lm<-lm(meanweightmicro~treatment,data=m3)
Anova(lm, type=2) 
#no treatment effect

#windows();plot(lm) #fine

#Month 3 females
lm<-lm(meanweightmicro~treatment,data=f3)
Anova(lm, type=2) 
#no treatment effect

#plot(lm) #fine


#### SEX RATIOS ####

#get data in
sex <-read.csv("sex ratios.csv",stringsAsFactors = FALSE,
                strip.white = TRUE, na.strings = c("NA","") )

sex<-sex%>%
  mutate(ratioftom=females/(females+males))

#orders the treatments in the order that we want them
sex$treatment <- factor(sex$treatment, 
                         levels = c("control", "5 day", "10 day", "15 day"))

summary_sex<-sex%>%
  group_by(treatment) %>%
  summarise(mean = mean(ratioftom),
            n=n(),
            sd = sd(ratioftom),
            se = sd / sqrt(n))

#Plot
ratiosplot<-ggplot(summary_sex, aes(x = treatment, y = mean, color=treatment)) +
  geom_point(size = 4,position = position_dodge(width = 0.6)) + #mean points
  geom_point(data = sex, aes(x=treatment,y = ratioftom, color = treatment), #all points
             position = position_jitterdodge(jitter.width = 0.3,
                                             dodge.width  = 0),
             size = 3, alpha=0.5,shape=16) +
  scale_color_manual(values = c("#894400","#c46200", "#ffb56b", "#ffd1a1"), 
                     name = "Heatwave duration", 
                     labels=c("15 days", "10 days", "5 days", "Control"),
                     breaks = c("15 day","10 day","5 day","control"),
                     guide = guide_legend(reverse=TRUE))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                position = position_dodge(width = 0.6),
                width = 0) +
  expand_limits(y = 0)+
  xlab("Heatwave duration") +
  ylab("Sex ratio (proportion female)")+
  scale_x_discrete (labels=c("15 days", "10 days", "5 days", "Control"),
                    breaks = c("15 day","10 day","5 day","control"))+
  theme_tess()

#windows();ratiosplot

ggsave(file="Figures/sex ratio.pdf",
       ratiosplot, width = 20, height = 17, units = "cm")

#analysis

lm<-lm(log(ratioftom) ~ treatment, data=sex)
Anova(lm, type=2) 

#check assumptions
#plot(lm) #looks ok - one outlier but it's real (ratio of 0.1)

