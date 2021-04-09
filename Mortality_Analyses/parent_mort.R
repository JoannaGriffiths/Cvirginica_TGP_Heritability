
setwd("~/LSU/Research/Oyster Heritability/Oyster_Heritability_Manuscript_boxfolder")
library(lme4)
library(car)


mort <- read.delim("parent_mort.txt", header=T)
mort.model = lm(Mort ~ Site, data = mort)
Anova(mort.model,test.statistic = "F") 

cumul_mort <- read.delim("cumulative_parent_mort.txt", header=T)

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = (sd(x[[col]]))/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df3 <- data_summary(cum_mort, varname="Mortality", 
                    groupnames=c("Site", "Cumul_date"))
head(df3)



windows()
ggplot(data=df3, aes(x=Cumul_date, y=Mortality, group=Site)) +
  geom_line(aes(color=Site))+
  geom_point(aes(color=Site)) +
  geom_errorbar(data=df3, aes(ymin=Mortality-se, ymax=Mortality+se), width=.2)+
  scale_color_grey() +  theme_classic() +
  labs(y="Mortality (%)", x="Date")+
  scale_x_discrete(breaks=c(seq(1, 8, by = 1)))
                     labels=c("1" = "Dose 0.5", "2" = "Dose 1","3" = "Dose 2", "4" = "Dose 0.5", "5" = "Dose 1","6" = "Dose 2", "7" = "Dose 0.5","8" = "Dose 0.5"))


ggplot(data=dat5, aes(x=Parental_sal, y=Length, fill= Larval_sal)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
  geom_errorbar(data=dat5, aes(ymin=Length-ster, ymax=Length+ster), width=.2, position=position_dodge(.7)) +
  scale_fill_manual(values=c("#3399FF", "#FF9900")) +
  labs(y=expression("Length (µm)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
  coord_cartesian(ylim = 65:81)