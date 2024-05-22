# Emma Walker
# E.M.Walker@exeter.ac.uk
# 22/01/2021

### Merge Cohorts  ====

library(dplyr)

setwd("")

## load in individual cohort objects ====

print("loading KCL pheno and meth data...")
load("../KCL_NormalisedBetas_NewStartDate.rdata")
kclbetas <- betas
kclpheno <- pheno
kclpheno$Institute <- rep("KCL")
rm(betas)
rm(pheno)

print("loading NL pheno and meth data...")
load("../NL_NormalisedBetas_NewStartDate.rdata")
nlbetas <- betas
nlpheno <- pheno
rm(betas)
rm(pheno)


## Join pheno data ====

nl <- nlpheno %>% dplyr::select(Basename, Individual_ID, ID.Number, TptLetters, Institute, Age, Sex, CD8T, CD4T, Bcell, Gran, Mono, NK, Collection.dates, Clozapine.start.date, Response.group, days, SmokingScore)
colnames(nl) <- c("Basename", "Individual_ID", "ID.Number", "Visit", "Institute", "Age", "Sex", "CD8T", "CD4T", "Bcell", "Gran", "Mono", "NK", "Collection.dates", "Clozapine.start.date", "Response.group", "days", "SmokingScore")

kcl <- kclpheno %>% dplyr::select(Basename, id, ID.Number., visit, Institute, Age, predictedSex, CD8T, CD4T, Bcell, Gran, Mono, NK,  Collection.dates, Clozapine.start.date, Response.group, days, SmokingScore)
colnames(kcl) <- c("Basename", "Individual_ID", "ID.Number", "Visit", "Institute", "Age", "Sex", "CD8T", "CD4T", "Bcell", "Gran", "Mono", "NK", "Collection.dates", "Clozapine.start.date", "Response.group", "days", "SmokingScore")



## Join the 2 cohorts pheno data together ====

pheno <- rbind(kcl, nl)



## Join the betas together ====

# keep only probes that are in both datasets
head(print("store as a list"))
data<-list("NL" = nlbetas, "KCL" = kclbetas) ## add in all cohorts
ncohorts<-length(data)

mincohorts<-ncohorts ## minimum number of cohorts probe present in to include in analysis

print("create vector of probes to include in meta-analysis")
probeCount<-table(unlist(lapply(data, rownames)))
probes<-names(probeCount[which(probeCount >= mincohorts)])

print("subset beta files to only contain probes in both")
nlbetas <- nlbetas[row.names(nlbetas) %in% probes,]
kclbetas <- kclbetas[row.names(kclbetas) %in% probes,]

# check they are in the same order then join together
all.equal(row.names(nlbetas), row.names(kclbetas))
#[1] TRUE

betas <- cbind(kclbetas, nlbetas)
dim(betas)
#[1] 355845    126
dim(pheno)
#[1] 126  18

## save the new objects to an rdata file ====

save(pheno, betas, file = "Merged_CLZ_data.rdat")


