library(MCMCglmm)

Ped<-as.data.frame(read.table(file="oyster_ped_covar_matrix.txt",header=TRUE)) 
Ped$animal<-as.factor(Ped$animal)
Ped$FATHER<-as.factor(Ped$FATHER)
Ped$MOTHER<-as.factor(Ped$MOTHER)
head(Ped)

family_meansv2 <- read.csv("family_means_v2.csv", header=T)
family_meansv2$animal<-as.factor(family_meansv2$Animal)
family_meansv2$Sire_Acclim<-as.factor(family_meansv2$Sire_Acclim)
family_meansv2$Dam_Acclim<-as.factor(family_meansv2$Dam_Acclim)
family_meansv2$Sire_ID<-as.factor(family_meansv2$Sire_ID)
family_meansv2$Dam_ID<-as.factor(family_meansv2$Dam_ID)
family_meansv2$H_Length <- family_meansv2$H_Length * 1000
family_meansv2$L_Length <- family_meansv2$L_Length * 1000


prior1.1<-list(G=list(G1=list(V=diag(2)*(0.002/1.002),n=1.002), #number in diag(X) is the number of response variables (ie Salinity is 2 because its either H or L)
                      G2=list(V=diag(2)*(0.002/1.002),n=1.002)),
               R=list(V=diag(2)*(0.002/1.002),n=1.002)) #the number in this list depends on how many variables in rcov

cov_mcmc.model <- MCMCglmm(fixed = cbind(L_Length,H_Length) ~ Dam_Acclim + Sire_Acclim + trait, 
                           random=~us(trait):animal + idh(trait):Dam_ID, 
                           nitt=1300000, thin=50, burnin=300000,
                           rcov = ~us(trait):units,
                           pedigree=Ped, data=family_meansv2, prior=prior1.1, 
                           family=c("gaussian", "gaussian"), verbose=FALSE)

posterior.mode(cov_mcmc.model$VCV)
HPDinterval(cov_mcmc.model$VCV[,"traitH_Length:traitL_Length.animal"],0.95)
HPDinterval(cov_mcmc.model$VCV[,"traitH_Length:traitH_Length.animal"],0.95)
HPDinterval(cov_mcmc.model$VCV[,"traitL_Length:traitL_Length.animal"],0.95)

