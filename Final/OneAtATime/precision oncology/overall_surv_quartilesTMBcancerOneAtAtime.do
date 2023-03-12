log using "overall_surv_quartilesTMBcancerOneAtAtime.log", replace
* overall_surv_quartilesTMBcancerOneAtAtime.log
* This program creates graphs for the supplement of the second submission to 
* Precision Oncology
* 

* Explore the effects of the cancer type and cancer type plus TMB
* for each cancer type evaluated separately.

* In these analyses, the TMB covariates consist of cubic spline covariates from the 
* 3-knot model of SMC and the 4-knot model of CNAF

use ../../../build.dta
drop if OverallSurvivalMonths == . | overallfate == .
gen followup = OverallSurvivalMonths/12

* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
stset follow15, failure(fate15) scale(1)

gen logPM = log(PointMutationCount)
mkspline _Spmc = logPM, cubic nknots(3) displayknots
mkspline _Scna = CNAFractionGenomeAltered, cubic nknots(4) displayknots

tabulate cancer

* Analyze the effect of TMB on survival in each cancer type separately

local figNo = 4
forvalues i = 1/32 {
    if `i' != 1 & `i' != 11 & `i' != 15 & `i' != 30  {
        preserve 
        local cancer_label: label value_lab `i'
        di "Analyze the effect of TMB on survival for `cancer_label'"
        keep if cancer == `i'
        tabulate cancer
    
        stcox _S*
        estimates stats
        di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
        predict loghr, xb
        centile loghr , centile(25 50 75)
        local c1 = r(c_1)
        local c2 = r(c_2)
        local c3 = r(c_3)
        sum loghr
        local max = r(max)
        gen loghr4 = recode(loghr, `c1',`c2',  `c3', `max')
        egen quartile = group(loghr4)
        tabulate quartile
        sts test quartile, logrank
        local p =  chi2tail(`r(df)', `r(chi2)')
        di "p = `p'"
        local p: di %8.1g = `p'
        di "p = `p' cancer i = `i'
        local pbf = `p'*32
        di "pbf = `pbf'
        local pbf: di %8.1g = `pbf'
        if `pbf' > 0.05 {
            local pbf  NS
            di "`pbf'"
        } 
        local figNo = `figNo'+1
        di "fig no = `figNo'"
* Draw a survival plot associated with cancer-type omitting the i-th cancer
    sts graph, by(quartile) ci  							///
        risktable(, rowtitle(1st quartile) 	group(#1) size(small)) 			///
        risktable(, rowtitle(2nd quartile)  group(#2) size(small)) 			///
        risktable(, rowtitle(3rd quartile)  group(#3) size(small))			///
        risktable(, rowtitle(4th quartile)  group(#4) size(small)) 			///
        plot1opts(lcolor(blue) lwidth(thick) ) 					///
        plot2opts(lcolor(red) lwidth(thick))					///
        plot3opts(lcolor(green)lwidth(thick)) 					///
        plot4opts(lcolor(black) lwidth(thick))  					///
        ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30)  lwidth(none)) 		///
        ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			///
        ytitle(Probability of overall survival) 	 				///
        ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
        title("", size(zero)) 							///
        subtitle("`cancer_label'" "P = `p'" "Bonferroni corrected P = `pbf'", 	///
            size(small) position(12) ring(1)) 					///
        legend(subtitle("Quartiles") order(9 "1st" 10 "2nd" 11 "3rd" 12 "4th") 	///	, size(small)
            symxsize(*0.5) cols(1) position(8) ring(1))  graphregion(color(white)) 	///
        xsize(7) ysize(7) scale(.7) name(quartile_cancer_`i')					///
        note(" "									///
        "{bf:Supplementary Figure S`figNo'.} Kaplan-Meier plots of overall survival for patients with" ///
            "`cancer_label'." ///
        "Patients are divided into quartiles based on their hazard ratio from the SM+CNA" "model." ///
            ,size(small) position(7) ring(3))					
        graph export quartile_graph_cancer_`i'.pdf, name(quartile_cancer_`i') replace  
        graph export quartile_graph_cancer_`i'.emf, name(quartile_cancer_`i') replace  
        drop loghr loghr4 quartile
        restore
    }
}
log close
translate "overall_surv_quartilesTMBcancerOneAtAtime.log" "overall_surv_quartilesTMBcancerOneAtAtime.pdf", replace
