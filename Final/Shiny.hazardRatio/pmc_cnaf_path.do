log using pmc_cnaf_path.log,  replace
* pmc_cnaf_path.log 

* Strip build.dta to the minimum size for a Shiny app.
* Eliminate records with missing values

*pause on
use ../../build.dta, clear
drop if OverallSurvivalMonths == . | overallfate == .  
gen followup = OverallSurvivalMonths/12

* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
* the haven library does not seem to be handling variables with value labels properly
* drop the value label for cancer
label drop value_lab
list if cancer == .

keep PointMutationCount CNAFractionGenomeAltered cancer  DiagnosisAge CancerStage TumorGrade fate15 follow15
drop if PointMutationCount ==  . | CNAFractionGenomeAltered == .
 order follow15 fate15, first
 
 #delimit ;
 label define cancer 
 1    "Adrenocortical Carcinoma"
 2    "Bladder Cancer"
 3    "Breast Cancer"
 4    "Cervical Squamous Cell Carcinoma"
 5    "Cholangiocarcinoma"
 6    "Colorectal Adenocarcinoma"
 7    "Diffuse Large B-Cell Lymphoma"
 8    "Esophageal Adenocarcinoma"
 9    "Glioblastoma Multiforme"
 10   "Head and Neck Squamous Cell Carcinoma"
 11   "Kidney Chromophobe"
 12   "Kidney Renal Clear Cell Carcinoma"
 13   "Kidney Renal Papillary Cell Carcinoma"
 14   "Acute Myeloid Leukemia"
 15   "Brain Lower Grade Glioma"
 16   "Liver Hepatocellular Carcinoma"
 17   "Lung Adenocarcinoma"
 18   "Lung Squamous Cell Carcinoma"
 19   "Mesothelioma"
 20   "Ovarian Serous Cystadenocarcinoma"
 21   "Pancreatic Adenocarcinoma"
 22   "Pheochromocytoma and Paraganglioma"
 23   "Prostate Adenocarcinoma"
 24   "Sarcoma"
 25   "Skin Cutaneous Melanoma"
 26   "Stomach Adenocarcinoma"
 27   "Testicular Germ Cell Tumors"
 28   "Thyroid Carcinoma"
 29   "Thymoma"
 30   "Uterine Corpus Endometrial Carcinoma"
 31   "Uterine Carcinosarcoma"
 32   "Uveal Melanoma" ;
 #delimit cr
 label values cancer cancer
   
 codebook
 
 * Add two extra rows that will be used to calculate hazard ratios
 * Place the values PointMutationCount = 1 and CNAFractionGenomeAltered =0
 * in the penultimate row.
 
 local np2 = _N +2
 set obs `np2'
 replace PointMutationCount       = 1 if _n == _N -1
 replace CNAFractionGenomeAltered = 0 if _n == _N -1
 
save pmc_cnaf_path.dta, replace
gen logPM = log(PointMutationCount)

stset follow15, failure(fate15) scale(1)
mkspline _Spm = logPM, cubic nknots(3) displayknots
mkspline _Scna = CNAFractionGenomeAltered, cubic nknots(4) displayknots

stcox _S*, nohr

stcox logPM CNAFractionGenomeAltered

* Estimate the hazard ratio associated with PMC = 10 and CNAF = .5
replace logPM = log(10) if _n==_N
replace CNAFractionGenomeAltered = .5 if _n==_N
predict loghr, xb
predict se, stdp
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
gen hr = exp(loghr)
list hr lb ub loghr logPM CNAFractionGenomeAltered if _n > _N-3

stcox _S* i.cancer, nohr

stcox _S* i.cancer DiagnosisAge, nohr

preserve
drop if cancer==14 | cancer == 32 
stcox _S* i.cancer DiagnosisAge, nohr
restore

sort cancer
by cancer: tabulate CancerStage fate15, missing
stcox _S* i.cancer  i.CancerStage, nohr

stcox _S*   i.cancer DiagnosisAge i.CancerStage, nohr

by cancer: tabulate TumorGrade CancerStage, missing
stcox _S*   i.cancer i.TumorGrade i.CancerStage, nohr

stcox _S*   i.cancer DiagnosisAge i.TumorGrade i.CancerStage, nohr

preserve
collapse (firstnm) PointMutationCount CNAFractionGenomeAltered DiagnosisAge CancerStage TumorGrade , by(cancer)
drop if cancer == .
list
foreach var in PointMutationCount CNAFractionGenomeAltered DiagnosisAge CancerStage TumorGrade  {
    gen ok`var' = "Yes"
    replace ok`var' = "No" if `var' == .
}
rename okPointMutationCount 		SMcount
rename okCNAFractionGenomeAltered	cnaf 
rename okDiagnosisAge 			age
rename okCancerStage 			stage
rename okTumorGrade			grade

list cancer SMcount cnaf age stage grade
restore

use ../../build.dta, clear
drop if OverallSurvivalMonths == . | overallfate == .

sort cancer
gen denominator = .
forvalues i = 1/32 {
    preserve
    keep if cancer == `i'
    local denom = _N
    restore
replace denominator = `denom' if cancer == `i'
}
*pause
* Create indicator valiables for patients in models A throgh G of Supplemental Table 3
gen inA = 100/denominator if PointMutationCount !=  . & CNAFractionGenomeAltered != . // Model with TMB
gen inB = 100/denominator if inA != . & cancer != .		// Model with TMB and cancer
gen inC = 100/denominator if inB != . & DiagnosisAge != .	// Model with TMB, cancer and age
gen inD = 100/denominator if inB != . & CancerStage != .		// Model with TMB, cancer and stage
gen inE = 100/denominator if inD != . & DiagnosisAge != .	// Model with TMB, cancer, age and stage
gen inF = 100/denominator if inD != . & TumorGrade != .		// Model with TMB, cancer, stage and grade
gen inG = 100/denominator if inF != . & DiagnosisAge != .	// Model with TMB, cancer, age, stage and grade

collapse (sum) inA inB inC inD inE inF inG , by(cancer)
format in* %5.1f
list
log close   
*translate hgrader_PMC_overall_surv.log hr_PMC_overall_surv.pdf, replace

