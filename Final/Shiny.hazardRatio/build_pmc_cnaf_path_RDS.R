# build_pmc_cnaf_path_RDS.R

# Create pmc_cnaf_path.RDS for use with shiny app

library(haven)
(pmc_cnaf_path <- read_dta("pmc_cnaf_path.dta"))
saveRDS(pmc_cnaf_path, file = "pmc_cnaf_path.RDS")
pmc_cnaf_path <- readRDS("pmc_cnaf_path.RDS")
 pmc <- pmc_cnaf_path$PointMutationCount
 cnaf <- pmc_cnaf_path$CNAFractionGenomeAltered
 cancer <-pmc_cnaf_path$cancer
 age <-pmc_cnaf_path$DiagnosisAge
 stage <-pmc_cnaf_path$CancerStage
 grade <-pmc_cnaf_path$TumorGrade
 logPM <- log(pmc)
 x <- cbind(logPM,cnaf,cancer,age,stage, grade)
 head(x)
 nrow(x)

