# Emma Walker
# E.M.Walker@exeter.ac.uk
# 11/08/2020

#Zhang age predictor https://github.com/qzhang314/DNAm-based-age-predictor/blob/master/data.age

############################## format input files for pred.R


### format input_file

#input_file: a R object file which contains the DNA methylation information for 
#individuals (N * M matrix). N is the number of individuals and M is the number of 
#CpG sites. Beta value is used as DNA methylation measurement.


### format age files

#age_file: an input file which has two column: individual ID and real 
#chronological age. Please note, the first line should be the header.

library(dplyr)

# load in QC'd data and create age file
setwd("")

# load KCL betas
load("../KCL_NormalisedBetas_NewStartDate.rdata")
saveRDS(object = betas, file = "~/clozapine/Zhang_age_pred/AgePredAll3cohorts/KCL1_betas.rds") # save betas as RDS file
age_file <- pheno %>% dplyr::select(Basename, Actual.Age)
names(age_file) <-c("ID", "age") #change names to be the same as in example age file
write.table(age_file, "~/clozapine/Zhang_age_pred/AgePredAll3cohorts/KCL1.age") #needs to be saved as a table

rm(list=ls())


# load NL betas
load("../NL_NormalisedBetas_NewStartDate.rdata")
saveRDS(object = betas, file = "~/clozapine/Zhang_age_pred/AgePredAll3cohorts/NL_betas.rds") # save betas as RDS file
age_file <- pheno %>% dplyr::select(Basename, Age)
names(age_file) <-c("ID", "age") #change names to be the same as in example age file
write.table(age_file, "~/clozapine/Zhang_age_pred/AgePredAll3cohorts/NL.age") #needs to be saved as a table

rm(list=ls())


# use this on command line to run predictor
# module load R/3.6.2-foss-2019b
# cd ~/clozapine/Zhang_age_pred/AgePredAll3cohorts/
# Rscript pred.R -i KCL1_betas.rds -o KCL1.age.pred -a KCL1.age 
# Rscript pred.R -i NL_betas.rds -o NL.age.pred -a NL.age 


