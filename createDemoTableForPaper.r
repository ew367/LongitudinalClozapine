# Emma Walker
# E.M.Walker@exeter.ac.uk
# 10/07/2023

# make demo table for paper


# read in pheno file

setwd("")
load("rdataFiles/Merged_CLZ_data.rdat")

#number of samples
table(pheno$Institute)
#     KCL Utrecht 
#     92      34 

KCL <- pheno[which(pheno$Institute == "KCL"),]
NL <- pheno[which(pheno$Institute == "Utrecht"),]

# number of individuals
length(unique(KCL$Individual_ID))
# 26
length(unique(NL$Individual_ID))
# 12

# sex percentages
table(KCL$Sex)

#F  M 
#22 70 

table(NL$Sex)

# F  M 
# 6 28 

getSex <- KCL[!duplicated(KCL$Individual_ID), ] 
table(getSex$Sex)
#F  M 
#6 20 

getSex <- NL[!duplicated(NL$Individual_ID), ] 
table(getSex$Sex)
#F  M 
#2 10 


## age

kclBaseline <- KCL[which(KCL$Visit == "a"),]
mean(as.numeric(kclBaseline$Age))
sd(as.numeric(kclBaseline$Age))


nlBaseline <- NL[which(NL$Visit == "a"),]
mean(as.numeric(nlBaseline$Age))
sd(as.numeric(nlBaseline$Age))


# number of time points
table(KCL$Visit)
table(NL$Visit)
