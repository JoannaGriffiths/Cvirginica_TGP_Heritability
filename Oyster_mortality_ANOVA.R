
library(lme4)
library(car)
library(emmeans)

setwd("~/Documents/Oyster Heritability Data")

mort=read.table("oyster_mortality_LC_VB.txt", header = T)
mort$Cross <- as.factor(mort$Cross)
mort$Salinity <- as.factor(mort$Salinity)
mort$Sire <- as.factor(mort$Sire)
mort$Dam <- as.factor(mort$Dam)

attach(mort)
shapiro.test(X.survival)  # pvalue less than 0.05, data not normal

sqrtmort=sqrt(X.survival) #log and sqrt makes pvalue smaller/worse
shapiro.test(sqrtmort)

lgmort=log(X.survival) #log and sqrt makes pvalue smaller/worse
shapiro.test(lgmort)

mort.model1 = lmer(X.survival ~ Salinity*Sire*Dam + (1|Cross), data = mort)
summary(mort.model1)

Anova(mort.model1,test.statistic = "F") 
'''
                       F Df Df.res Pr(>F)
Salinity          0.1564  1 28.258 0.6955
Sire              0.6370  1 36.856 0.4299
Dam               0.0183  1 37.010 0.8932
Salinity:Sire     1.2155  1 28.381 0.2795
Salinity:Dam      1.0103  1 28.829 0.3232
Sire:Dam          1.2531  1 37.096 0.2702
Salinity:Sire:Dam 0.0213  1 29.298 0.8849
'''

