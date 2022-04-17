log using overall_surv_quartilesTMBwithinStage.log, replace
* overall_surv_quartilesTMBwithinStage.log

* 

* Explore the effects of the TMB model within stage on survival. 

* In these analyses, the TMB covariates consist of cubic spline covariates from the 
* 3-knot model of SMC and the 4-knot model of CNAF

use ../build.dta
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


* Within subgroups defined by stage
*Draw a survival plot associated with SMC and CNA. This is based on the TMB model

preserve
keep if CancerStage == 1
di "Analysis restricted to cancer stage 1"
tabulate CancerStage
stcox _S*
estimates stats

di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))

test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 
di "P-value for PMC and CNA = "r(p)

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
di "p = " chi2tail(`r(df)', `r(chi2)')
local p = chi2tail(`r(df)', `r(chi2)')
local p: di %8.1g = `p'
if `p' == 0 {
    local p = "P < 1.0e-323"
}    
else {
    local p = "P = `p'"
}    

sts test quartile if quartile <3, logrank
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test quartile if quartile != 1 & quartile < 4, logrank
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test quartile if quartile > 2 & quartile != ., logrank
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test quartile if quartile ==2 | quartile == 4 , logrank
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test quartile if quartile !=3 & quartile != . , logrank
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test quartile if quartile ==1 | quartile ==4  , logrank
di "p = " chi2tail(`r(df)', `r(chi2)')

sts graph, by(quartile) ci  						///
    risktable(, rowtitle(1st quartile) 	group(#1) size(medsmall)) 		///
    risktable(, rowtitle(2nd quartile)  group(#2) size(medsmall)) 		///
    risktable(, rowtitle(3rd quartile)  group(#3) size(medsmall))		///
    risktable(, rowtitle(4th quartile)  group(#4) size(medsmall)) 		///
    plot1opts(lcolor(blue) lwidth(thick) ) 					///
    plot2opts(lcolor(red) lwidth(thick))					///
    plot3opts(lcolor(green)lwidth(thick)) 					///
    plot4opts(lcolor(black) lwidth(thick))  				///
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30)  lwidth(none)) 		///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			///
    ytitle(Probability of overall survival) 	 			///
    ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis)  		///
    xlabel(0(2)14) title("", size(zero)) 					///
    note("Shaded regions denote 95% CIs" ,ring(0) position(4) size(medsmall)) ///
    subtitle( "Stage 1 cancers, `p'",	///
        size(medsmall) position(2) ring(0)) 					///
    legend(subtitle("Quartiles") order(9 "1st" 10 "2nd" 11 "3rd" 12 "4th") 	///	
        symxsize(*0.5) cols(1) position(7) ring(0)) graphregion(color(white)) ///
    xsize(5) ysize(4) name(quartile_TMB_stage_1)
graph export quartile_graph_TMB_stage_1.pdf, name(quartile_TMB_stage_1) replace  
drop loghr loghr4 quartile
restore
forvalues i = 2/4 {
    preserve
    keep if CancerStage == `i'
    di "Analysis restricted to cancer stage `i'"
    tabulate CancerStage
    stcox _S*
    estimates stats
    
    di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
    
    test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 
    di "P-value for PMC and CNA = "r(p)
    
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
    di "p = " chi2tail(`r(df)', `r(chi2)')
    local p = chi2tail(`r(df)', `r(chi2)')
    local p: di %8.1g = `p'
    if `p' == 0 {
        local p = "P < 1.0e-323"
    }    
    else {
        local p = "P = `p'"
    }    

    sts test quartile if quartile <3, logrank
    di "p = " chi2tail(`r(df)', `r(chi2)')
    sts test quartile if quartile != 1 & quartile < 4, logrank
    di "p = " chi2tail(`r(df)', `r(chi2)')
    sts test quartile if quartile > 2 & quartile != ., logrank
    di "p = " chi2tail(`r(df)', `r(chi2)')
    sts test quartile if quartile ==2 | quartile == 4 , logrank
    di "p = " chi2tail(`r(df)', `r(chi2)')
    sts test quartile if quartile !=3 & quartile != . , logrank
    di "p = " chi2tail(`r(df)', `r(chi2)')
    sts test quartile if quartile ==1 | quartile ==4  , logrank
    di "p = " chi2tail(`r(df)', `r(chi2)')
    if `i' == 1 {
        local ci "Shaded regions denote 95% CIs."
    }
    else {
        local ci " "
    }
    if `i' == 4 {
        local pos = 1
    }
    else {
        local pos = 7
    }
    sts graph, by(quartile) ci  						///
        risktable(, rowtitle(1st quartile) 	group(#1) size(medsmall)) 	///
        risktable(, rowtitle(2nd quartile)  group(#2) size(medsmall)) 		///
        risktable(, rowtitle(3rd quartile)  group(#3) size(medsmall))		///
        risktable(, rowtitle(4th quartile)  group(#4) size(medsmall)) 		///
        plot1opts(lcolor(blue) lwidth(thick) ) 					///
        plot2opts(lcolor(red) lwidth(thick))					///
        plot3opts(lcolor(green)lwidth(thick)) 					///
        plot4opts(lcolor(black) lwidth(thick))  				///
        ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30)  lwidth(none)) 		///
        ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			///
        ytitle(Probability of overall survival) 	 			///
        ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis)  		///
        xlabel(0(2)14) title("", size(zero)) 					///
        subtitle( "Stage `i' cancers, `p'",					///
            size(medsmall) position(2) ring(0)) 				///
        legend(off) graphregion(color(white)) 					///
        xsize(5) ysize(4) name(quartile_TMB_stage_`i')
    graph export quartile_graph_TMB_stage_`i'.pdf, name(quartile_TMB_stage_`i') replace  
    drop loghr loghr4 quartile
    restore
}

log close
translate overall_surv_quartilesTMBwithinStage.log overall_surv_quartilesTMBwithinStage.pdf, replace
