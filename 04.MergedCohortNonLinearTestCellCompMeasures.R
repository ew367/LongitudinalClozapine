# Emma Walker
# E.M.Walker@exeter.ac.uk


#### Test changes in cell composition as factor of time on clozapine
####

#a = baseline
#b = 6 weeks from clozapine initiation 
#c = 8 weeks from c i  (only taken if low clozapine levels at 6 weeks)
#d = 12 weeks
#e = 6 months
#Exact times vary and are recorded individually but that's the rough schedule.


# load packages
library(lme4)
library(lmrTest)

print("loading pheno and methylation data...")
setwd("")
load("rdataFiles/Merged_CLZ_data.rdat")

cellMeasures<-c("CD8T", "CD4T", "NK", "Bcell", "Gran", "Mono")

out<-matrix(data = NA, nrow = length(cellMeasures), ncol = 10)
colnames(out)<-c("linear_Beta", "linear_SE", "linear_P", "sq_Beta", "sq_SE", "sq_P", "sq_linear_Beta", "sq_linear_SE", "sq_linear_P", "ANOVA_Pr(>Chisq)")
rownames(out)<-cellMeasures

print("changing model terms to be right type...")
pheno$days <- as.numeric(as.character(pheno$days))
pheno$Individual_ID <- as.factor(pheno$Individual_ID)
pheno$Age <- as.numeric(as.character(pheno$Age))
pheno$Institute <- as.factor(pheno$Institute)

#add in time squared object
timesq <- pheno$days*pheno$days

# run regression analysis on each cell type
for(each in cellMeasures){
  model<-lmer(pheno[,each] ~ timesq + pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute), REML = FALSE)
  model.null<-lmer(pheno[,each] ~ pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute), REML = FALSE)
  P<-anova(model,model.null)["model","Pr(>Chisq)"]
  out[each,1:3]<-summary(model.null)$coefficients["pheno$days", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[each,4:6]<-summary(model)$coefficients["timesq", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[each,7:9]<-summary(model)$coefficients["pheno$days", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[each, "ANOVA_Pr(>Chisq)"] <- P
}

write.csv(out, "Results/MergedNonLinearCellCompositionMeasuresChangeWithClozapineUseNewStart.csv")



################## delta cell composition plots ====================================================

  # remove samples with no a visit and recalculate timesq
  pheno <- pheno[-c(which(pheno$Individual_ID == "0044")),] 
  betas <- betas[,which(colnames(betas) %in% pheno$Basename)]
  timesq <- pheno$days*pheno$days
  
  #vals for reg line
  xvals <- seq(min(pheno$days), max(pheno$days), 1)
  xvalsq <- xvals*xvals
  
  
  
  pdf("Plots/MergedNonLinearCellCompositionMeasuresChangeWithClozapineUseNewStartForPaper12thJune.pdf")
  
  # run just for Bcells as this is the only cell type that showed a significant change over time
  cell<- "Bcell"
  
  # option to loop over all cell types
  #par(mar=c(5,5,4,1)+.1)
  #for(cell in cellMeasures){
    
    #rerun the model 
    mod <- lmer(pheno[,cell] ~ timesq + pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute), REML = FALSE)
    
    # add in the residuals
    var<-pheno$days*fixef(mod)["pheno$days"]+timesq*fixef(mod)["timesq"]+residuals(mod)
    names(var) <- pheno$Basename
    
    #calculate delta composition
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
    yvals<-xvals*fixef(mod)["pheno$days"]+xvalsq*fixef(mod)["timesq"]
    
    # create dataframe for ggplot
    plotdf <- as.data.frame(cbind(pheno$days, pheno$Individual_ID, deltacomps))
    colnames(plotdf) <- c("days", "ID", "deltacomps")
    
    
    # and for the trendline (needs a seperate dataframe)
    linedf<-as.data.frame(cbind(xvals, yvals))
    
    # create plot of change in cell composition over time
    P <- ggplot() + 
      geom_line(aes(days, deltacomps, group=ID), plotdf, colour = "grey")+
      geom_point(aes(days, deltacomps), plotdf)+
      geom_line(aes(xvals, yvals), linedf, size =2)+ # trendline
      xlab("Time (weeks)")+
      ylab(paste0(cell, " (proportion)"))+
      scale_x_continuous(breaks = seq(0,270,42),minor_breaks = seq(from = 0, to = 270, by = 7), labels = seq(0,39,6))+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    
  #}
    

  
  print(P)
    
dev.off()  
  