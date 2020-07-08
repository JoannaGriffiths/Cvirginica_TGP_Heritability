
library(nlme)
library(car)
library(emmeans)

setwd("~/LSU/Research/Oyster Heritability/phys_data")
biochem <- read.delim("biochem_data.txt", header = T)

attach(biochem)
hist(Protein, main="Histogram")
shapiro.test(Protein)
shapiro.test(HC)
shapiro.test(TG)
shapiro.test(FA)
shapiro.test(CH)
shapiro.test(PL)

lgProtein=log(Protein)
lgHC=log(HC)
lgTG=log(TG)
lgFA=log(FA)
lgCH=log(CH)
lgPL=log(PL)

shapiro.test(lgProtein) #logs look better than sqrt
shapiro.test(lgHC)
shapiro.test(lgTG)
shapiro.test(lgFA)
shapiro.test(lgCH)
shapiro.test(lgPL)


protein.model = lm(lgProtein ~ Site*Pop, data=biochem)
protein.model = lm(lgProtein ~ Site, data=biochem)
aov1 = aov(protein.model)
summary(aov1)
 
HC.model = lm(lgHC ~ Site*Pop, data=biochem)
HC.model = lm(lgHC ~ Site, data=biochem)
aov1 = aov(HC.model)
summary(aov1)

TG.model = lm(lgTG ~ Site*Pop, data=biochem)
TG.model = lm(lgTG ~ Site, data=biochem)
aov1 = aov(TG.model)
summary(aov1)

FA.model = lm(lgFA ~ Site*Pop, data=biochem)
FA.model = lm(lgFA ~ Site, data=biochem)
aov1 = aov(FA.model)
summary(aov1)

CH.model = lm(lgCH ~ Site*Pop, data=biochem)
CH.model = lm(lgCH ~ Site, data=biochem)
aov1 = aov(CH.model)
summary(aov1)

PL.model = lm(lgPL ~ Site*Pop, data=biochem)
PL.model = lm(lgPL ~ Site, data=biochem)
aov1 = aov(PL.model)
summary(aov1)

#testing for effect of whether this female successfully spawned
protein.model = lm(lgProtein ~ Spawn, data=biochem)
aov1 = aov(protein.model)
summary(aov1)

HC.model = lm(lgHC ~ Spawn, data=biochem)
aov1 = aov(HC.model)
summary(aov1)

TG.model = lm(lgTG ~ Spawn, data=biochem)
aov1 = aov(TG.model)
summary(aov1)

FA.model = lm(lgFA ~ Spawn, data=biochem)
aov1 = aov(FA.model)
summary(aov1)

CH.model = lm(lgCH ~ Spawn, data=biochem)
aov1 = aov(CH.model)
summary(aov1)

PL.model = lm(lgPL ~ Spawn, data=biochem)
aov1 = aov(PL.model)
summary(aov1)



Data <- read.delim("oyster_size_eggQual_VB_LC.txt", header = T)
attach(Data)
lglength <- log(Data$Length)
lgProtein <- log(Data$Protein)
lgHC <- log(Data$HC)
lgTG <- log(Data$TG)
lgFA <- log(Data$FA)
lgCH <- log(Data$CH)
lgPL <- log(Data$PL)

plot(Data$Protein, Data$Length)
abline(fit <- (lm(Data$Length ~ Data$Protein)))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)
aov(fit)

plot(lgProtein, lglength)
abline(fit <- (lmer(lglength ~ lgProtein + (1|DamID) + (1|SireID), data = Data)))
protein.model = lmer(lglength ~ lgProtein + (1|DamID) + (1|SireID), data = Data)
Anova(protein.model,test.statistic = "F")

Data <- read.delim("oyster_size_eggQual_VB_LC.txt", header = T)
plot(Data$HC, Data$Length)
abline(fit <- (lm(Data$Length ~ Data$HC)))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)
aov(fit)
HC.model = lmer(lglength ~ lgHC + (1|DamID) + (1|SireID), data = Data)
Anova(HC.model,test.statistic = "F")

Data <- read.delim("oyster_size_eggQual_VB_LC.txt", header = T)
plot(Data$TG, Data$Length)
abline(fit <- (lm(Data$Length ~ Data$TG)))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)
aov(fit)
TG.model = lmer(lglength ~ lgTG + (1|DamID) + (1|SireID), data = Data)
Anova(TG.model,test.statistic = "F")

Data <- read.delim("oyster_size_eggQual_VB_LC.txt", header = T)
plot(Data$FA, Data$Length)
abline(fit <- (lm(Data$Length ~ Data$FA)))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)
aov(fit)
FA.model = lmer(lglength ~ lgFA + (1|DamID) + (1|SireID), data = Data)
Anova(FA.model,test.statistic = "F")

Data <- read.delim("oyster_size_eggQual_VB_LC.txt", header = T)
plot(Data$CH, Data$Length)
abline(fit <- (lm(Data$Length ~ Data$CH)))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)
aov(fit)
CH.model = lmer(lglength ~ lgCH + (1|DamID) + (1|SireID), data = Data)
Anova(CH.model,test.statistic = "F")

Data <- read.delim("oyster_size_eggQual_VB_LC.txt", header = T)
plot(Data$PL, Data$Length)
abline(fit <- (lm(Data$Length ~ Data$PL)))
legend("topright", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))
summary(fit)
aov(fit)
PL.model = lmer(lglength ~ lgPL + (1|DamID) + (1|SireID), data = Data)
Anova(PL.model,test.statistic = "F")
