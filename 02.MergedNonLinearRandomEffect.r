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

print("loading lme4")
library(lme4)
print("loading lmertest")
library(lmerTest)
print("loading dplyr")
library(dplyr)

print("loading pheno and meth data...")
setwd("")
load("rdataFiles/Merged_CLZ_data.rdat")

#ensure betas and pheno are in the same order:
print(identical(as.character(pheno$Basename), colnames(betas)))
#[1] TRUE


#create object to store regression results
out<-matrix(data = NA, nrow = nrow(betas), ncol = 10)
colnames(out)<-c("linear_Beta", "linear_SE", "linear_P", "sq_Beta", "sq_SE", "sq_P", "sq_linear_Beta", "sq_linear_SE", "sq_linear_P", "ANOVA_Pr(>Chisq)")
rownames(out)<-rownames(betas)

print("changing model terms to be right type...")
pheno$days <- as.numeric(as.character(pheno$days))
pheno$Individual_ID <- as.factor(pheno$Individual_ID)
pheno$Age <- as.numeric(as.character(pheno$Age))
pheno$Institute <- as.factor(pheno$Institute)

#add in time squared object
timesq <- pheno$days*pheno$days

print("running models")
for(i in 1:nrow(betas)){
  full.model<-lmer(betas[i,] ~ timesq + pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute) + pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Gran,	REML = FALSE)
  null.model<-lmer(betas[i,] ~ pheno$days + pheno$Age + pheno$Sex + (1|pheno$Individual_ID) + (1|pheno$Institute) + pheno$CD8T + pheno$CD4T + pheno$NK + pheno$Bcell + pheno$Gran,	REML = FALSE)
  out[i,7:9]<-summary(full.model)$coefficients["pheno$days", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[i,4:6]<-summary(full.model)$coefficients["timesq", c("Estimate", "Std. Error", "Pr(>|t|)")]
  out[i,1:3]<-summary(null.model)$coefficients["pheno$days", c("Estimate", "Std. Error", "Pr(>|t|)")]
  P<-anova(full.model,null.model)["full.model","Pr(>Chisq)"]
  out[i, "ANOVA_Pr(>Chisq)"] <- P
}


write.csv(out, "Results/MergedNotAnnoNonLinear.csv")
out <- read.csv("Results/MergedNotAnnoNonLinear.csv", row.names = 1)

## annotate eith gene info etc. using illumina manifest
load("../../EWAS/AllProbeIlluminaAnno.Rdata")
probeAnnot<-probeAnnot[rownames(out),]
out<-cbind(out, probeAnnot[,c("CHR", "MAPINFO", "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP")])


out<-out[order(out[,3]),] # order by linear model P value
out2<-out[order(out[,6]),] # order by square term P value


write.csv(out, "Results/MergedMethylationChangeWithClozapineUse_ControlForCellComp_NewStart_BatchRandomEffectLinearTerm.csv")
write.csv(out2, "Results/MergedNonLinearMethylationChangeWithClozapineUse_ControlForCellComp_NewStart_BatchRandomEffectLinearTerm.csv")

