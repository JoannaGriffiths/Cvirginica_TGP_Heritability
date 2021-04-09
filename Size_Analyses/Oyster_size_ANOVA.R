

setwd("~/LSU/Research/Oyster Heritability")

library(lme4)
library(car)
library(emmeans)

################################adult sizes reared at LUM and GI
setwd("~/Documents/Oyster Heritability Data")


growth <- read.table("parent_oyster_size.txt", header = T)
attach(growth)
shapiro.test(Size_mm)

growth.model1 = lm(Size_mm ~ Site, data = growth)
summary(aov(growth.model1))
'''
             Df Sum Sq Mean Sq F value Pr(>F)    
Site          1  51589   51589   313.6 <2e-16 ***
Residuals   278  45734     165                   
'''

#graph LUM vs GI acclimation site on adult size
LUM=subset(growth$Size_mm, growth$Site=="LUM")
mean(LUM) #88.22
sd(LUM) #13.87
sd(LUM)/sqrt(length(LUM)) #1.24

GI=subset(growth$Size_mm, growth$Site=="GI")
mean(GI) #115.51
sd(GI) #11.91
sd(GI)/sqrt(length(GI)) #0.96

dat <- data.frame(
  Salinity = factor(c("LUMCON", "Grand Isle")),
  Length = c(88.22, 115.51),
  stdev = c(13.87, 11.91),
  ster = c(1.24,0.96)
)


quartz()
ggplot(data=dat, aes(x=Salinity, y=Length, fill= Salinity)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
  geom_errorbar(data=dat, aes(ymin=Length-ster, ymax=Length+ster), width=.2, position=position_dodge(.7)) +
  scale_fill_manual(values=c("#FF9900", "#3399FF")) +
  labs(y=expression("Length (mm)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
  coord_cartesian(ylim = 75:120)

#############################larvae data below
  
growth <- read.table("oyster_size.txt", header = T)
attach(growth)
shapiro.test(Length)
lglength=log(Length)
shapiro.test(lglength)
  
hist(lglength, main="Histogram of Growth")
  
growth.model1.6 = lmer(lglength ~ Salinity*Sire*Dam + (1|DamID) + (1|SireID), data = growth)
summary(growth.model1.6)
Anova(growth.model1.6,test.statistic = "F") 
'''
                        F Df  Df.res    Pr(>F)    
Salinity          39.4959  1 2607.63 3.841e-10 ***
Sire               0.1015  1   12.11   0.75541    
Dam                0.8766  1   17.17   0.36212    
Salinity:Sire      3.6816  1 2611.38   0.05512 .  
Salinity:Dam      43.7088  1 2613.01 4.603e-11 ***
Sire:Dam          35.0812  1 1271.40 4.065e-09 ***
Salinity:Sire:Dam  2.1108  1 2609.67   0.14638    
'''

lsmeans(growth.model1.6, list(pairwise ~ Salinity*Dam), adjust = "tukey") 
'''
$`pairwise differences of Salinity, Dam`
 contrast      estimate          SE      df t.ratio p.value
H,H - L,H  0.004634881 0.007203095 2602.76   0.643  0.9178
H,H - H,L -0.074288720 0.049484550   17.59  -1.501  0.4577
H,H - L,L  0.012795790 0.049630170   17.79   0.258  0.9938
L,H - H,L -0.078923600 0.049563354   17.67  -1.592  0.4078
L,H - L,L  0.008160910 0.049708866   17.87   0.164  0.9984
H,L - L,L  0.087084510 0.009823901 2615.76   8.865  <.0001
'''

lsmeans(growth.model1.6, list(pairwise ~ Salinity*Sire), adjust = "tukey")
'''
$`pairwise differences of Salinity, Sire`
 contrast      estimate          SE      df t.ratio p.value
H,H - L,H  0.054587242 0.007617812 2600.85   7.166  <.0001
H,H - H,L -0.042828497 0.046736709   12.81  -0.916  0.7967
H,H - L,L -0.005696349 0.047004264   13.15  -0.121  0.9993
L,H - H,L -0.097415739 0.046710702   12.79  -2.086  0.2094
L,H - L,L -0.060283591 0.046934939   13.07  -1.284  0.5879
H,L - L,L  0.037132148 0.009288071 2617.51   3.998  0.0004
'''


#looking for effects of non-additive genetic variation
growth.model1.7 = lm(lglength ~ DamID*SireID, data = growth)
summary(aov(growth.model1.7))
  '''
              Df Sum Sq Mean Sq F value Pr(>F)    
DamID          21  23.30  1.1094   63.19 <2e-16 ***
SireID         11   6.39  0.5808   33.09 <2e-16 ***
DamID:SireID    6   2.37  0.3946   22.48 <2e-16 ***
Residuals    2588  45.43  0.0176                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
'''

#check that block does not exlpain residuals in our model. It does not, thus Block was not included as a random effect in model
model1.6.resid <- resid(growth.model1.6)
windows()
plot(growth$Block, model1.6.resid, ylab="Residuals", xlab="Block")

  
  #graph low vs high dam acclimation salinity and rearing salinity at low vs high
  low_sal_low_dam=subset(growth$Length, growth$Salinity=="L" & growth$Dam=="L")
  mean(low_sal_low_dam) #0.074
  sd(low_sal_low_dam) #0.01094428
  sd(low_sal_low_dam)/sqrt(length(low_sal_low_dam)) #0.0005211562
  
  low_sal_high_dam=subset(growth$Length, growth$Salinity=="L" & growth$Dam=="H")
  mean(low_sal_high_dam) #0.07467286
  sd(low_sal_high_dam) # 0.01318167
  sd(low_sal_high_dam)/sqrt(length(low_sal_high_dam)) #0.0004640167
  
  high_sal_low_dam=subset(growth$Length, growth$Salinity=="H" & growth$Dam=="L")
  mean(high_sal_low_dam) #0.0807464
  sd(high_sal_low_dam) # 0.01228544
  sd(high_sal_low_dam)/sqrt(length(high_sal_low_dam)) #0.0005210186
  
  high_sal_high_dam=subset(growth$Length, growth$Salinity=="H" & growth$Dam=="H")
  mean(high_sal_high_dam) #0.07584933
  sd(high_sal_high_dam) #0.0137566
  sd(high_sal_high_dam)/sqrt(length(high_sal_high_dam)) #0.000479526
  
  dat5 <- data.frame(
    Salinity = factor(c("LL", "LH", "HL", "HH")),
    Parental_sal = factor(c("Dam Low", "Dam High", "Dam Low", "Dam High")),
    Larval_sal = factor(c("Offspring Low", "Offspring Low", "Offspring High", "Offspring High"), levels=c("Offspring Low","Offspring High")),
    Length = c(74.0, 74.67, 80.74, 75.85),
    stdev = c(10.94, 13.18, 12.2, 13.73),
    ster = c(0.52, 0.46, 0.52, 0.479)
  )
  
  quartz()
  ggplot(data=dat5, aes(x=Parental_sal, y=Length, fill= Larval_sal)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
    geom_errorbar(data=dat5, aes(ymin=Length-ster, ymax=Length+ster), width=.2, position=position_dodge(.7)) +
    scale_fill_manual(values=c("#3399FF", "#FF9900")) +
    labs(y=expression("Length (Âµm)"), x="", fill="") +
    #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
          strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
    coord_cartesian(ylim = 65:81)
  
  
  #graph low vs high sire acclimation salinity and rearing salinity at low vs high
  low_sal_low_dam=subset(growth$Length, growth$Salinity=="L" & growth$Sire=="L")
  mean(low_sal_low_dam) #0.07705947
  sd(low_sal_low_dam) #0.0130551
  sd(low_sal_low_dam)/sqrt(length(low_sal_low_dam)) #0.0006127064
  
  low_sal_high_dam=subset(growth$Length, growth$Salinity=="L" & growth$Sire=="H")
  mean(low_sal_high_dam) #0.07293451
  sd(low_sal_high_dam) # 0.01181811
  sd(low_sal_high_dam)/sqrt(length(low_sal_high_dam)) #0.0004194092
  
  high_sal_low_dam=subset(growth$Length, growth$Salinity=="H" & growth$Sire=="L")
  mean(high_sal_low_dam) #0.07854545
  sd(high_sal_low_dam) #0.01284704
  sd(high_sal_low_dam)/sqrt(length(high_sal_low_dam)) #0.0005000704
  
  high_sal_high_dam=subset(growth$Length, growth$Salinity=="H" & growth$Sire=="H")
  mean(high_sal_high_dam) #0.07716134
  sd(high_sal_high_dam) #0.01385657
  sd(high_sal_high_dam)/sqrt(length(high_sal_high_dam)) #0.0005167629
  
  dat5 <- data.frame(
    Salinity = factor(c("LL", "LH", "HL", "HH")),
    Parental_sal = factor(c("Sire Low", "Sire High", "Sire Low", "Sire High")),
    Larval_sal = factor(c("Offspring Low", "Offspring Low", "Offspring High", "Offspring High"), levels=c("Offspring Low","Offspring High")),
    Length = c(77.05, 72.93, 78.94, 77.16),
    stdev = c(13.05, 11.81, 12.8, 13.85),
    ster = c(0.61, 0.42, 0.50, 0.516)
  )
  
  quartz()
  ggplot(data=dat5, aes(x=Parental_sal, y=Length, fill= Larval_sal)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
    geom_errorbar(data=dat5, aes(ymin=Length-ster, ymax=Length+ster), width=.2, position=position_dodge(.7)) +
    scale_fill_manual(values=c("#3399FF", "#FF9900")) +
    labs(y=expression("Length (Âµm)"), x="", fill="") +
    #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
          strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
    coord_cartesian(ylim = 65:81)
  
  #graph low vs high rearing salinity
  
  low_sal=subset(growth$Length, growth$Salinity=="L")
  mean(low_sal) #0.0744351
  sd(low_sal) #0.01243681
  sd(low_sal)/sqrt(length(low_sal)) #0.0003520478
  
  high_sal=subset(growth$Length, growth$Salinity=="H")
  mean(high_sal) #0.07782379
  sd(high_sal) #0.0133959
  sd(high_sal)/sqrt(length(high_sal)) #0.0003607376
  
  dat <- data.frame(
    Salinity = factor(c("Low", "High")),
    Length = c(74.43, 77.82),
    stdev = c(12.43, 13.39),
    ster = c(0.35,0.36)
  )
  
  
  quartz()
  ggplot(data=dat, aes(x=Salinity, y=Length, fill= Salinity)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
    geom_errorbar(data=dat, aes(ymin=Length-ster, ymax=Length+ster), width=.2, position=position_dodge(.7)) +
    scale_fill_manual(values=c("#FF9900", "#3399FF")) +
    labs(y=expression("Length (Âµm)"), x="", fill="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
          strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
    coord_cartesian(ylim = 65:80)
  
  
  
  
