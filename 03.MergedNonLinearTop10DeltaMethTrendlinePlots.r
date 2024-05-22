# Emma Walker
# E.M.Walker@exeter.ac.uk


##### Plot change in methylation over time and add trendline from regression analyses


library(lme4)

print("loading pheno and meth data...")
setwd("")
load("rdataFiles/Merged_CLZ_data.rdat")

#load in the results data (out2 = ordered by squared term for time pval, out = ordered by linear term p val)
out2 <- read.csv("Results/MergedNonLinearMethylationChangeWithClozapineUse_ControlForCellComp_NewStart_BatchRandomEffectLinearTerm.csv", row.names = 1)


#check they are in the same order 
pheno$Basename <- as.character(pheno$Basename)
identical(pheno$Basename, colnames(betas))


# remove samples with no a visit
pheno <- pheno[-c(which(pheno$Individual_ID == "0044")),] 
betas <- betas[,which(colnames(betas) %in% pheno$Basename)]

#make sure betas and pheno are in the same order:
print(identical(as.character(pheno$Basename), colnames(betas)))
#[1] TRUE

#changing model terms to be right type
pheno$days <- as.numeric(as.character(pheno$days)) # this was originally in date format
pheno$Individual_ID <- as.factor(pheno$Individual_ID)
pheno$Age <- as.numeric(as.character(pheno$Age))
pheno$Institute <- as.factor(pheno$Institute)

#add in time squared object
timesq <- pheno$days*pheno$days


#vals for reg line
xvals <- seq(min(pheno$days), max(pheno$days), 1)
xvalsq <- xvals*xvals


#### plots

#### Nonlinear - plot top 10 most signifcant probes

pdf("Plots/Merged_DeltaMethNonLinearTop10_WeeksTEST12thJune.pdf")
par(mar=c(5,5,4,1)+.1)
for(i in 1:10){
  
  
  #rerun the model 
  probe<-lmer(unlist(betas[rownames(out2)[i],]) ~ timesq + pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute) + pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Gran,	REML = FALSE)
  
  var<-pheno$days*fixef(probe)["pheno$days"]+timesq*fixef(probe)["timesq"]+residuals(probe)  
  
  #calculate deltabetas so that each visit is calculcated as change to the baseline (a visit = 0)
  varbetas <- t(as.data.frame(var))
  deltabetas <- vector()
  for(each in unique(pheno$Individual_ID)){
    indPheno <- pheno[pheno$Individual_ID == each,]
    indBetas <- varbetas[,colnames(varbetas) %in% indPheno$Basename]
    aVisit <- as.character(indPheno$Basename[indPheno$Visit == "a"])
    indBetas <- indBetas - indBetas[aVisit]
    deltabetas <- c(deltabetas, indBetas) 
  }  
  
  
  # order the pheno and betas files to be the same again (above code chunk puts them in different orders)
  deltabetas <- deltabetas[order(match(names(deltabetas), pheno$Basename))]
  print(identical(pheno$Basename, names(deltabetas))) # check they are in the same order
  
  #calculate yvals for trendline (multiply by 100 to convert to %)
  yvals<-(xvals*fixef(probe)["pheno$days"]+xvalsq*fixef(probe)["timesq"])*100
  
  
  #create dataframes of data for the plots
  plotdf<- as.data.frame(cbind(pheno$days, pheno$Individual_ID, pheno$Institute, deltabetas*100))
  colnames(plotdf) <- c("days", "ID", "batch", "deltabetas")
  plotdf$batch <- as.factor(plotdf$batch)
  
  # and for the trendline (needs a seperate dataframe)
  linedf<-as.data.frame(cbind(xvals, yvals))
  
  
  P <- ggplot() + 
    geom_line(aes(days, deltabetas, group=ID), plotdf, colour = "grey")+
    geom_point(aes(days, deltabetas), plotdf)+
    geom_line(aes(xvals, yvals), linedf, size =2)+ # trendline
    xlab("Time (weeks)")+
    ylab("% DNA Methylation")+
    ggtitle(paste(rownames(out2)[i], "-", out2$UCSC_REFGENE_NAME[i], "Effect size =", signif(out2$sq_Beta[i]*100*7, 2),"P value =", signif(out2$sq_P[i],2)))+
    scale_x_continuous(breaks = seq(0,270,42),minor_breaks = seq(from = 0, to = 270, by = 7), labels = seq(0,39,6))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
  
  print(P)
}



dev.off()
