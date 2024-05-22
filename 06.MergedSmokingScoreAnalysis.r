# Emma Walker
# E.M.Walker@exeter.ac.uk
# 13/08/2020

#### Test changes in smoking score as factor of time on clozapine
####

library(lme4)
library(lmerTest)
library(dplyr)

print("loading pheno and meth data...")
setwd("")

#load in version of pheno file with smoking score
load("rdataFiles/MergedCohortPhenoWithSmoking.rdata")

print("changing model terms to be right type...")
pheno$days <- as.numeric(as.character(pheno$days))
pheno$Individual_ID <- as.factor(pheno$Individual_ID)
pheno$Age <- as.numeric(as.character(pheno$Age))
pheno$Institute <- as.factor(pheno$Institute)
pheno$SmokingScore <- as.numeric(pheno$SmokingScore)

#add in time squared object
timesq <- pheno$days*pheno$days

#run model
full.model<-lmer(pheno$SmokingScore ~ timesq + pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute),	REML = FALSE)
null.model<-lmer(pheno$SmokingScore ~ pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute),	REML = FALSE)




################## delta plot ====================================================

pdf("Plots/MergedPlotSmokingChangeWithClozapineExposure.pdf")

# remove samples with no a visit and recalculate timesq
pheno <- pheno[-c(which(pheno$Individual_ID == "0044")),] 

#vals for reg line
xvals <- seq(min(pheno$days), max(pheno$days), 1)


par(mar=c(5,5,4,1)+.1)

  
  #rerun the model 
mod<-lmer(pheno$SmokingScore ~ pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute),  REML = FALSE)
  
  var<-pheno$days*fixef(mod)["pheno$days"]+residuals(mod)
  names(var) <- pheno$Basename
  
  #calculate delta smoking score
  varcomps <- t(as.data.frame(var))
  deltacomps <- vector()
  for(each in unique(pheno$Individual_ID)){
    indPheno <- pheno[pheno$Individual_ID == each,]
    indComps <- varcomps[,colnames(varcomps) %in% indPheno$Basename]
    aVisit <- as.character(indPheno$Basename[indPheno$Visit == "a"])
    indComps <- indComps - indComps[aVisit]
    deltacomps <- c(deltacomps, indComps) 
  }  
  
  
  # order the pheno and betas files to be the same again (above caluculation of change in methylation puts them out of order)
  pheno$Basename <- as.character(pheno$Basename)
  pheno <- pheno[order(pheno$Basename),]
  deltacomps <- deltacomps[order(names(deltacomps))]
  print(identical(pheno$Basename, names(deltacomps)))
  
  #calculate yvals for trendline
  yvals<-xvals*fixef(mod)["pheno$days"]
  
  #create object for plot
  plotdf <- as.data.frame(cbind(pheno$days, pheno$Individual_ID, deltacomps))
  colnames(plotdf) <- c("days", "ID", "deltacomps")
  
  
  # and for the trendline (needs a seperate dataframe)
  linedf<-as.data.frame(cbind(xvals, yvals))
  
  
  P <- ggplot() + 
    geom_line(aes(days, deltacomps, group=ID), plotdf, colour = "grey")+
    geom_point(aes(days, deltacomps), plotdf)+
    geom_line(aes(xvals, yvals), linedf, size =2)+ # trendline
    xlab("Time (weeks)")+
    ylab("Smoking Score")+
    scale_x_continuous(breaks = seq(0,270,42),minor_breaks = seq(from = 0, to = 270, by = 7), labels = seq(0,39,6))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  

print(P)

dev.off()
