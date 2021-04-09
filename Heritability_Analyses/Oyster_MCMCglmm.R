
library(MCMCglmm)

setwd("~/Documents/Oyster Heritability Data")

Data<-as.data.frame(read.table(file="oyster_size.txt",header=TRUE))
Data$animal<-as.factor(Data$animal)
Data$Sire<-as.factor(Data$Sire)
Data$Dam<-as.factor(Data$Dam)
Data$SireID<-as.factor(Data$SireID)
Data$DamID<-as.factor(Data$DamID)
Data$Salinity<-as.factor(Data$Salinity)
Data$Pop<-as.factor(Data$Pop)

Data$Length2 <- Data$Length * 1000

#select either low or high salinity
Data <- subset(Data, Salinity=="L")
Data <- subset(Data, Salinity=="H")

Ped<-as.data.frame(read.table(file="oyster_ped.txt",header=TRUE)) 
Ped$animal<-as.factor(Ped$animal)
Ped$FATHER<-as.factor(Ped$FATHER)
Ped$MOTHER<-as.factor(Ped$MOTHER)
head(Ped)

# specify priors, but tell the model to pay little attention to it since we want priors to come from our data (n=1)
p.var<-var(Data$Length,na.rm=TRUE)
p.var<-var(Data$Length2,na.rm=TRUE)
prior1.1<-list(G=list(G1=list(V=matrix(p.var/2),n=1)), R=list(V=matrix(p.var/2),n=1))

#now we can run the model
model1.1<-MCMCglmm(Length~1,random=~animal, pedigree=Ped,data=Data,prior=prior1.1)

#now look at output to see how confident we are in results: The important point here is that a consistent amount of
#variation around a largely unchanging mean value of the intercept was obtained, and
#the posterior distribution of the intercept appears to be valid.The plot on the left shows a time series of
#the values of 1000 samples of the posterior distribution of the the model intercept (mean birthweight).
plot(model1.1$Sol)

# this plots distributions of the estimates of the additive genetic (animal) and residual (units) effects.
plot(model1.1$VCV)

#Because the values in the left-hand plots appear to have different values at the beginning of the run,
#we might suspect that a longer burn-in period might be required. We can reduce the
#autocorrelation by lengthening the rest of the run and sampling the chain less frequently.
model1.1<-MCMCglmm(Length~1,random=~animal,
                   pedigree=Ped,data=Data,
                   nitt=65000,thin=50,burnin=15000,
                   prior=prior1.1,verbose=FALSE)

# now notice in plots that the autocorrelation is reduced
plot(model1.1$Sol)
plot(model1.1$VCV)

#A more compact way to evaluate the validity of the posterior distributions is to calculate autocorrelation among samples, as follows:
#the autocorrelation for all lag values greater than zero should be near zero.
#*****Mine are not that close to 0, not sure what to do about this, maybe adding he fixed and random effects will lower this?
autocorr(model1.1$VCV)

#We can obtain estimates of the additive genetic and residual variance by calculating the modes of the posterior distributions:
posterior.mode(model1.1$VCV)

#We can obtain the Bayesian equivalent of confidence intervals by calculating the the
#values of the estimates that bound 95% (or any other proportion) of the posterior distributions:
HPDinterval(model1.1$VCV)

#We specified weak priors in this analyses. Now we will check whether or not these priors infuenced the results that we obtained. 
#The simplest way to do this is to rerun the model with different priors. Rather than splitting the observed phenotypic variance
#evenly between the priors for the genetic and residual effects, we will construct priors
#with the same degree of belief parameters, but we will specify a much larger proportion of genetic control:
prior1.1.2<-list(G=list(
  G1=list(V=matrix(p.var*0.95),n=1)),
  R=list(V=matrix(p.var*0.05),n=1))
model1.1.2<-MCMCglmm(Length~1,random=~animal,
                     pedigree=Ped,data=Data,prior=prior1.1.2,
                     nitt=65000,thin=50,burnin=15000,verbose=FALSE)
posterior.mode(model1.1$VCV)

posterior.mode(model1.1.2$VCV)
#if posterior results are similar, prior had little effect, but mine are more different than in the tutorial, how different is acceptable?

#obtain an estimate of the heritability by applying the
#basic formula h2=VA/VP to each sample of the posterior disribution:
posterior.heritability1.1<-model1.1$VCV[,"animal"]/(model1.1$VCV[,"animal"]+model1.1$VCV[,"units"])
HPDinterval(posterior.heritability1.1,0.95)
posterior.mode(posterior.heritability1.1)

#Generate a plot of the posterior distribution of this heritability estimate
plot(posterior.heritability1.1)


#dam and sire acclimation as fixed effects
prior1.1<-list(G=list(G1=list(V=matrix(p.var/2),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.1<-MCMCglmm(Length~Sire+Dam,random=~animal,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.1,verbose=FALSE)
posterior.mode(model1.1$VCV)

posterior.mode(model1.1$Sol[,"DamL"])
HPDinterval(model1.1$Sol[,"DamL"],0.95) #crosses zero for VB_LC low sal!  Does not cross zero for VB_LC high sal!
posterior.mode(model1.1$Sol[,"SireL"])
HPDinterval(model1.1$Sol[,"SireL"],0.95) #crosses zero for VB_LC low sal! Crosses zero for VB_LS high sal!

##just sire acclimation as fixed effect
prior1.2<-list(G=list(G1=list(V=matrix(p.var/2),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.2<-MCMCglmm(Length~Sire,random=~animal,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.2,verbose=FALSE)
posterior.mode(model1.2$VCV)

posterior.mode(model1.1$Sol[,"SireL"])
HPDinterval(model1.1$Sol[,"SireL"],0.95) #crosses zero for VB_LC at low sal! crosses zero for VB_LC at high sal! 

#just dam acclimation as fixed effect
prior1.3<-list(G=list(G1=list(V=matrix(p.var/2),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.3<-MCMCglmm(Length~Dam,random=~animal,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.3,verbose=FALSE)
posterior.mode(model1.3$VCV)

posterior.mode(model1.1$Sol[,"DamL"])
HPDinterval(model1.1$Sol[,"DamL"],0.95) #crosses zero for VB_LC low sal!  Does not cross zero for VB_LC at high sal!

#sire and dam acclimation as fixed effects, DamID as random 
prior1.4<-list(G=list(G1=list(V=matrix(p.var/3),n=1),
                      G2=list(V=matrix(p.var/3),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.4<-MCMCglmm(Length2~Sire+Dam,random=~animal+DamID,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.4,verbose=FALSE)
posterior.mode(model1.4$VCV)
model1.4new <- model1.4

#sire and dam acclimation as fixed effects, SireID as random
prior1.5<-list(G=list(G1=list(V=matrix(p.var/3),n=1),
                      G2=list(V=matrix(p.var/3),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.5<-MCMCglmm(Length~Sire+Dam,random=~animal+SireID,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.5,verbose=FALSE)
posterior.mode(model1.5$VCV)

#sire and dam acclimation as fixed effects, SireID and DamID as random
prior1.6<-list(G=list(G1=list(V=matrix(p.var/4),n=1),
                      G2=list(V=matrix(p.var/4),n=1),
                      G3=list(V=matrix(p.var/4),n=1)),
               R=list(V=matrix(p.var/4),n=1))
model1.6<-MCMCglmm(Length~Sire+Dam,random=~animal+SireID+DamID,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.6,verbose=FALSE)
posterior.mode(model1.6$VCV)

#dam acclimation as fixed effects, DamID as random 
prior1.7<-list(G=list(G1=list(V=matrix(p.var/3),n=1),
                      G2=list(V=matrix(p.var/3),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.7<-MCMCglmm(Length~Dam,random=~animal+DamID,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.7,verbose=FALSE)
posterior.mode(model1.7$VCV)

#dam acclimation as fixed effects, SireID as random
prior1.8<-list(G=list(G1=list(V=matrix(p.var/3),n=1),
                      G2=list(V=matrix(p.var/3),n=1)),
               R=list(V=matrix(p.var/3),n=1))
model1.8<-MCMCglmm(Length~Dam,random=~animal+SireID,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.8,verbose=FALSE)
posterior.mode(model1.8$VCV)

#dam acclimation as fixed effects, SireID and DamID as random
prior1.9<-list(G=list(G1=list(V=matrix(p.var/4),n=1),
                      G2=list(V=matrix(p.var/4),n=1),
                      G3=list(V=matrix(p.var/4),n=1)),
               R=list(V=matrix(p.var/4),n=1))
model1.9<-MCMCglmm(Length~Dam,random=~animal+SireID+DamID,
                   pedigree=Ped,data=Data,
                   nitt=1300000,thin=50,burnin=300000,
                   prior=prior1.9,verbose=FALSE)
posterior.mode(model1.9$VCV)

save(model1.1, model1.2, model1.3, model1.4, model1.4new, model1.7, model1.8, model1.9, file = "Oyster_HighSal_MCMCglmm_results.RData")
save(model1.1, model1.2, model1.3, model1.4, model1.4new, model1.7, model1.8, model1.9, file = "Oyster_LowSal_MCMCglmm_results.RData")

load("Oyster_LowSal_MCMCglmm_results.RData")
load("Oyster_HighSal_MCMCglmm_results.RData")

model1.1$DIC #dam+sire
model1.2$DIC #sire
model1.3$DIC #dam
model1.4$DIC #dam+sire, r=damID
model1.5$DIC #dam+sire, r=sireID
model1.6$DIC #dam+sire, r=damID+sireID
model1.7$DIC #dam, r=damID
model1.8$DIC #dam, r=sireID
model1.9$DIC #dam, r=damID+sireID 

posterior.heritability1.9 <- model1.9$VCV[,"animal"]/(model1.9$VCV[,"animal"]+model1.9$VCV[,"DamID"]+model1.9$VCV[,"SireID"]+model1.9$VCV[,"units"])
posterior.mode(posterior.heritability1.9)
HPDinterval(posterior.heritability1.9,0.95)

posterior.maternal1.9 <- model1.9$VCV[,"DamID"]/(model1.9$VCV[,"animal"]+model1.9$VCV[,"DamID"]+model1.9$VCV[,"SireID"]+model1.9$VCV[,"units"])
posterior.mode(posterior.maternal1.9)
HPDinterval(posterior.maternal1.9,0.95)


posterior.heritability1.8 <- model1.8$VCV[,"animal"]/(model1.8$VCV[,"animal"]+model1.8$VCV[,"SireID"]+model1.8$VCV[,"units"])
posterior.mode(posterior.heritability1.8)
HPDinterval(posterior.heritability1.8,0.95)

posterior.maternal1.8 <- model1.8$VCV[,"DamID"]/(model1.8$VCV[,"animal"]+model1.8$VCV[,"DamID"]+model1.8$VCV[,"SireID"]+model1.8$VCV[,"units"])
posterior.mode(posterior.maternal1.8)
HPDinterval(posterior.maternal1.8,0.95)

###Decided on model1.4
posterior.heritability1.4 <- model1.4$VCV[,"animal"]/(model1.4$VCV[,"animal"]+model1.4$VCV[,"DamID"]+model1.4$VCV[,"units"])
posterior.mode(posterior.heritability1.4)
HPDinterval(posterior.heritability1.4,0.95)

'''
low sal
0.6951874 
lower    upper
var1 0.4515146 0.865113

high sal
0.7111734 
lower    upper
0.4393718 0.8613924
'''

posterior.maternal1.4 <- model1.4$VCV[,"DamID"]/(model1.4$VCV[,"animal"]+model1.4$VCV[,"DamID"]+model1.4$VCV[,"units"])
posterior.mode(posterior.maternal1.4)
HPDinterval(posterior.maternal1.4,0.95)
'''
low sal
0.1001325
       lower     upper
var1 0.03283 0.3518339

high sal
0.1629563 
          lower     upper
var1 0.04245382 0.3563461
'''

#############determining variance components for table
posterior.mode(model1.4$VCV[,"animal"])
HPDinterval(model1.4$VCV[,"animal"],0.95)

posterior.mode(model1.4$VCV[,"DamID"])
HPDinterval(model1.4$VCV[,"DamID"],0.95)

posterior.mode(model1.4$Sol[,"SireL"])
HPDinterval(model1.4$Sol[,"SireL"],0.95)

posterior.mode(model1.4$Sol[,"DamL"])
HPDinterval(model1.4$Sol[,"DamL"],0.95)

posterior.mode((model1.4$VCV[,"animal"]+model1.4$VCV[,"DamID"]+model1.4$VCV[,"units"]))
HPDinterval((model1.4$VCV[,"animal"]+model1.4$VCV[,"DamID"]+model1.4$VCV[,"units"]),0.95)

posterior.mode(model1.4.L2$VCV[,"animal"])
HPDinterval(model1.4.L2$VCV[,"animal"],0.95)

posterior.mode(model1.4.L2$VCV[,"DamID"])
HPDinterval(model1.4.L2$VCV[,"DamID"],0.95)

posterior.mode(model1.4.L2$Sol[,"SireL"])
HPDinterval(model1.4.L2$Sol[,"SireL"],0.95)

posterior.mode(model1.4.L2$Sol[,"DamL"])
HPDinterval(model1.4.L2$Sol[,"DamL"],0.95)

posterior.mode((model1.4.L2$VCV[,"animal"]+model1.4.L2$VCV[,"DamID"]+model1.4.L2$VCV[,"units"]))
HPDinterval((model1.4.L2$VCV[,"animal"]+model1.4.L2$VCV[,"DamID"]+model1.4.L2$VCV[,"units"]),0.95)


############### breeder's equation for significance of heritability estimates
library("dplyr")

growth <- read.table("oyster_size.txt", header = T)
low_sal=subset(growth, growth$Salinity=="L")
low_sal$size <- low_sal$Length * 1000

low_upper = subset(low_sal, low_sal$Length > 0.065) #larvae smaller than 65 um would have to die to see shift in mean size at low salinity by 5um
mean(low_upper$Length) #0.0794

percentile <- ecdf(low_sal$Length)

percentile(0.065) #0.255, 25.5% smallest individuals die to evolve back to ambient salinity size
hist(low_sal$Length)

#plot histogram/density of larval sizes at low salinity
dat <- with(density(low_sal$size), data.frame(x, y))
windows()
ggplot(dat, mapping = aes(x = x, y = y)) + 
  geom_line(size = 1, color = "black")+
  geom_vline(aes(xintercept = 74.4),col='black',size=0.5)+
  geom_vline(aes(xintercept = 79.4),col='red',size=0.5)+
  geom_area(mapping = aes(x = ifelse(x>30 & x< 65 , x,0)), fill = "red") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA))+
  ylim(0,0.04)+
  scale_x_continuous(breaks=seq(20,150,10))
##########################

