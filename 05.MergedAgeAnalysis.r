# Emma Walker
# E.M.Walker@exeter.ac.uk

# epigentic age correlation and regression analysis

setwd("")

# load in age prediction files and combine
KCL <- read.csv("../agePred/KCL_age_pred.csv")
NL <- read.csv("../agePred/NL_age_pred.csv")

age <- rbind(KCL,NL)
age$Basename <- age$ID

cor.test(age$age, age$enpred)
# 0.9541674
cor.test(age$age, age$enpred)$p.value
# 8.069095e-67

# load pheno and methylation data and add in age prediction data
load("rdataFiles/Merged_CLZ_data.rdat")
pheno <- dplyr::left_join(pheno, age)


#### Test changes in predicted age as a factor of time on clozapine

#create object to store regression results
out<-matrix(data = NA, nrow = 1, ncol = 10)
colnames(out)<-c("linear_Beta", "linear_SE", "linear_P", "sq_Beta", "sq_SE", "sq_P", "sq_linear_Beta", "sq_linear_SE", "sq_linear_P", "ANOVA_Pr(>Chisq)")


print("changing model terms to be right type...")
pheno$days <- as.numeric(as.character(pheno$days))
pheno$Individual_ID <- as.factor(pheno$Individual_ID)
pheno$Age <- as.numeric(as.character(pheno$Age))
pheno$Institute <- as.factor(pheno$Institute)

#add in time squared object
timesq <- pheno$days*pheno$days

print("running models")
  model<-lmer(pheno$enpred ~ timesq + pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute), REML = FALSE)
  model.null<-lmer(pheno$enpred ~ pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute), REML = FALSE)
  P<-anova(model,model.null)["model","Pr(>Chisq)"]
  out[1,1:3]<-summary(model.null)$coefficients["pheno$days", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[1,4:6]<-summary(model)$coefficients["timesq", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[1,7:9]<-summary(model)$coefficients["pheno$days", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[1, "ANOVA_Pr(>Chisq)"] <- P



####  plots

# remove samples with no a visit and recalculate timesq
pheno <- pheno[-c(which(pheno$Individual_ID == "0044")),] 
#betas <- betas[,which(colnames(betas) %in% pheno$Basename)]
timesq <- pheno$days*pheno$days

#vals for reg line
xvals <- seq(min(pheno$days), max(pheno$days), 1)
xvalsq <- xvals*xvals


pdf("Plots/Merged_DeltaAgeLinearWeeks.pdf")
par(mar=c(5,5,4,1)+.1)

probe<-lmer(pheno$enpred ~ pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute), REML = FALSE)
  
var<-pheno$days*fixef(probe)["pheno$days"]+residuals(probe)
names(var) <- pheno$Basename
  
  #calculate deltabetas so that each visit is calculcated as change to the baseline (a visit = 0)
  deltabetas <- vector()
  for(each in unique(pheno$Individual_ID)){
    indPheno <- pheno[pheno$Individual_ID == each,]
    indVar <- var[names(var) %in% indPheno$Basename]
    aVisit <- as.character(indPheno$Basename[indPheno$Visit == "a"])
    indVar <- indVar - indVar[aVisit]
    deltabetas <- c(deltabetas, indVar) 
  }  
  
  
  # order the pheno and betas files to be the same again (above caluculation of change in methylation puts them out of order)
  deltabetas <- deltabetas[order(match(names(deltabetas), pheno$Basename))]
  print(identical(pheno$Basename, names(deltabetas))) # check they are in the same order
  
  #calculate yvals for trendline
  yvals<-(xvals*fixef(probe)["pheno$days"])
  
  #create dataframes of data for the plots
  plotdf<- as.data.frame(cbind(pheno$days, pheno$Individual_ID, pheno$Institute, deltabetas))
  colnames(plotdf) <- c("days", "ID", "batch", "deltabetas")
  plotdf$batch <- as.factor(plotdf$batch)
  
  # and for the trendline (needs a seperate dataframe)
  linedf<-as.data.frame(cbind(xvals, yvals))
  
  # plot
  P <- ggplot() + 
    geom_line(aes(days, deltabetas, group=ID), plotdf, colour = "grey")+
    geom_point(aes(days, deltabetas), plotdf)+
    geom_line(aes(xvals, yvals), linedf, size =2)+ # trendline
    xlab("Time (weeks)")+
    ylab("Epigenetic Age")+
    scale_x_continuous(breaks = seq(0,270,42),minor_breaks = seq(from = 0, to = 270, by = 7), labels = seq(0,39,6))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))

print(P)


dev.off()
