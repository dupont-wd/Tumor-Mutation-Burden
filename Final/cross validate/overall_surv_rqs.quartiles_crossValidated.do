log using overall_surv_rqs_quartiles_crossValidated.log, replace
* overall_surv_rqs_quartiles_crossValidated.log
* Generate Figure 4A

* Determine the combined effects of a 3-knot model of log PMC 
* and a 4-knot model of CNA on survival. 
*set more on
*pause on
use ../../build.dta
drop if OverallSurvivalMonths == . | overallfate == .
drop if PointMutationCount == . | CNAFractionGenomeAltered == .
gen followup = OverallSurvivalMonths/12
set seed 6519841

* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
stset follow15, failure(fate15) scale(1)

gen logPM = log(PointMutationCount)
mkspline _Spmc = logPM, cubic nknots(3) displayknots
mkspline _Scna = CNAFractionGenomeAltered, cubic nknots(4) displayknots

          
stcox _Spmc* _Scna* 
estimates stats
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))

test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 
di "P-value for PMC and CNA = "r(p)
more 

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

* Draw a survival plot associated with PMC and CNA
sts graph, by(quartile) ci  							///
    risktable(, rowtitle(1st quartile) 	group(#1) size(small)) 			///
    risktable(, rowtitle(2nd quartile)  group(#2) size(small)) 			///
    risktable(, rowtitle(3rd quartile)  group(#3) size(small))			///
    risktable(, rowtitle(4th quartile)  group(#4) size(small)) 			///
    plot1opts(lcolor(blue) lwidth(medthick) ) 					///
    plot2opts(lcolor(red) lwidth(medthick))					///
    plot3opts(lcolor(green)lwidth(medthick)) 					///
    plot4opts(lcolor(black) lwidth(medthick))  					///
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30)  lwidth(none)) 		///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			///
    ytitle(Probability of overall survival) 	 				///
    ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
    title("", size(zero)) 							///
    subtitle("Shaded regions denote 95% CIs." 					///
        "CNA and log PMC modeled with 4 and 3 knots, respectivley", 		///
        size(small) position(2) ring(0)) 					///
    legend(subtitle("Quartiles") order(9 "1st" 10 "2nd" 11 "3rd" 12 "4th") 	///	
        cols(1) position(8) ring(0))  graphregion(color(white)) 		///
    xsize(5) ysize(4) name(quartile_rqs2)
graph export quartile_graph_CNA&PMC.pdf, name(quartile_rqs2) replace  

* Divide data into 5 equal groups stratified by the 32 different cancers

gen id = _n  // Create a numeric id
gen rand = uniform()
sort cancer rand
by cancer: gen n=_n
gen group = 5 // We want 5 equal groups
forvalues i = 1/32 {  // For all 32 values of cancer
    centile n if cancer==`i', centile(20 40 60 80)
    forvalues j= 1/4 {
        local num = 5-`j'
        replace group = `num' if n <= r(c_`num') & cancer== `i'
    }
}
sort id
tabulate cancer group

* Calculate 5-fold cross-validated hazard ratios 

gen loghr_cv = .
gen loghr_tmp = .
gen hr = exp(loghr)
gen hr_cv = .

* Calculate reference line
gen refx = 1 if _n==1
gen refy=1 if _n== 1
replace refx = 20 if _n == 2
replace refy = 20 if _n == 2

 
forvalues i = 1/5 {

di  "Select knots based on all of the data except patients in group `i'"
    drop _S*
    di "Calculate logPM knots omitting group `i'"
    mkspline _Spmc = logPM if group !=`i', cubic nknots(3) displayknots
    matrix knots = r(knots)
    matrix list knots
    drop _Spmc*
    local knotPM1 = knots[1,1]
    local knotPM2 = knots[1,2]
    local knotPM3 = knots[1,3]

di  "Define logPM spline covariates for all of the data using the knots from" 
di  "all patients except patients in group `i'"
    mkspline _Spmc = logPM , cubic knots(`knotPM1' `knotPM2' `knotPM3') displayknots

di  "Repeat the above for CNAF"
    mkspline _Scna = CNAFractionGenomeAltered if group != `i'	///
        , cubic nknots(4) displayknots
    matrix knots = r(knots)
    matrix list knots
    local knotCN1 = knots[1,1]
    local knotCN2 = knots[1,2]
    local knotCN3 = knots[1,3]
    local knotCN4 = knots[1,4]
    drop _Scna*
di  "Define CNAF spline covariates for all of the data using the knots from" 
di  "all patients except patients in group i"
    
    mkspline _Scna = CNAFractionGenomeAltered 			///
        , cubic knots(`knotCN1' `knotCN2' `knotCN3' `knotCN4') displayknots
    
di  "Regress death against spline covariates for patients from group `i'"    
    stcox _Spmc* _Scna* if group !=`i'
    drop loghr_tmp 
    
di  "Predict log hazard for patients in group `i' using model "
di  "generated from all other patients"
    predict loghr_tmp if group ==`i', xb
    replace loghr_cv = loghr_tmp if group == `i'
    replace hr_cv = exp(loghr_tmp) if group == `i'

* Make a matrix of the parameter estimates
    matrix PAR`i' = e(b)
* List parameter estimates for the current group
di "Parameter estimates for group `i' are as follows"
    matrix list PAR`i'
    frame create grid`i'
    frame change grid`i'
    local col = 101
    local N= `col'^2
    set obs `N'
    local inclogPM =0.1
    gen cna = int((_n-1)/`col')/(`col'-1)
    gen logPM = mod((_n-1),`col')*`inclogPM'
        foreach j of numlist  1 10 100 1000 10000 {
        forvalues k = 1/9 {
            local tmp = `j'*`k'
            local `tmp'  = log(`tmp')
        }
    }
    mkspline _Spmc = logPM , cubic knots(`knotPM1' `knotPM2' `knotPM3') displayknots
    mkspline _Scna = cna 			///
        , cubic knots(`knotCN1' `knotCN2' `knotCN3' `knotCN4') displayknots

    predict loghr_cv, xb
    sum loghr_cv
    gen hr_cv = exp(loghr_cv)
    sum hr_cv

    label variable cna "CNA fraction"
    label variable hr "Mortal hazard ratio"

    twoway contour hr_cv cna logPM,   ccuts(1(1.5)19.5)  		///
        ccolors(blue "0 59 255" "0 120 255" cyan*1.4 cyan cyan*.6 	///
            green*1.8 green*1.4 green  )	                	///
        xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 				///
            `10' "10"  `30' "30"  `60' "60"				///
            `100' "100"  `300' "300"  `600' "600"			///
            `1000' "1000"  `3000' "3000"  `6000' "6000" 		///
            `10000' "10000" `20000' "20000" , angle(45))	       	///
        xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 			///
            `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" ///
            `200' "200" `400' "400" `500' "500" `700' "700" `800' "800" ///
            `900' "900" `2000' "2000" `4000' "4000" `5000' "5000" 	///
            `7000' "7000" `8000' "8000" `9000' "9000") 			///
        note("Unadjusted,  " 						///
            "CNA and log PMC modeled with 4 and 3 knots, respectively.  " ///
            "Model parameters from group `i'  " " ", ring(0) position(4)) ///
        scolor(blue) ecolor(red) xtitle( Point mutation count) 		///
        ylabel(0(.1)1, angle(0)) graphregion(color(white)) ysize(3) 	///
        xsize(3.75) name(contour`i')
    graph export contourPMC_CNA`i'.pdf, name(contour`i') replace 
    frame change default
}
gen diff = hr_cv - hr
save cross_validated.dta, replace
twoway (scatter diff hr, mcolor(blue) msymbol(smcircle_hollow) )	///
  , yline(0) xtitle(unvalidated hazard ratio) name(diff)		///
    ytitle( "cross-validated - unvalidated hazard ratios") 

regress hr_cv hr
predict yhat, xb
sort hr
twoway (scatter hr_cv hr, msymbol(smcircle_hollow) mcolor(blue)) 	///
    (line refy refx, lwidth(medthick) lcolor(red)) 			///
    (line yhat hr, lwidth(medthick) lcolor(cyan))			///	
  , xtitle(Unvalidated hazard ratio) 				///
    ytitle( Cross-validated hazard ratio)   			///
    legend(ring(0) position(10) subtitle(All patients) 		///
        order(1 "Observed" 2 "Reference" 3 "Expected") col(1)) 	///
    graphregion(color(white)) name(regress)
graph export cross_validated_regression.pdf, name(regress) replace 

    
* Regenerate the preceding scatter plots separately by group   
forvalues i = 1/5  {
    preserve
    keep if group == `i'
    * Calculate reference line
    replace refx = 1 if _n==1
    replace  refy=1 if _n== 1
    replace refx = 20 if _n == 2
    replace refy = 20 if _n == 2
    di "Group `i'"
     gen diff`i' = hr_cv - hr
    twoway (scatter diff`i' hr, mcolor(blue) msymbol(smcircle_hollow)) 	///
      , yline(0) xtitle(unvalidated hazard ratio) name(diff`i')	///
        ytitle( "cross-validated - unvalidated hazard ratios") 		///
        subtitle("Group `i'", ring(0) position(10) )
    regress hr_cv hr
    drop yhat
    predict yhat, xb
    sort hr
    twoway (scatter hr_cv hr, msymbol(smcircle_hollow) mcolor(blue)) 	///
        (line refy refx, lwidth(medthick) lcolor(red)) 			///
        (line yhat hr, lwidth(medthick) lcolor(cyan))			///	
      , xtitle(unvalidated hazard ratio) 				///
        ytitle( cross-validated hazard ratio)   			///
        legend(ring(0) position(10) subtitle(Group `i')			///
        order(1 "observed" 2 "reference" 3 "expected") col(1)) name(regress`i')
    restore    
}
frame change grid1
forvalues i = 2/5 {
    frameappend grid`i'
}
preserve
collapse (min) hr_cv, by(cna logPM)
    label variable cna "CNA fraction"
    label variable hr "Minimum cross-validated mortal hazard ratio"

twoway contour hr_cv cna logPM,   ccuts(1(1.5)19.5)  			///
    ccolors(blue "0 59 255" "0 120 255" cyan*1.4 cyan cyan*.6 green*1.8 ///
        green*1.4 green  )	                			///
    xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 				///
        `10' "10"  `30' "30"  `60' "60"					///
        `100' "100"  `300' "300"  `600' "600"				///
        `1000' "1000"  `3000' "3000"  `6000' "6000" 			///
        `10000' "10000" `20000' "20000" , angle(45))	       		///
    xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 			///
        `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 	///
        `200' "200" `400' "400" `500' "500" `700' "700" `800' "800" 	///
        `900' "900" `2000' "2000" `4000' "4000" `5000' "5000" `7000' 	///
        "7000" `8000' "8000" `9000' "9000") 				///
    note("Unadjusted,  " 						///
        "CNA and log PMC modeled with 4 and 3 knots, respectively.  " 	///
        "Minimum hazard from 5 cross-validation models  " 		///
        " ", ring(0) position(4)) 					///
    scolor(blue) ecolor(red) xtitle( Point mutation count) 		///
    ylabel(0(.1)1, angle(0)) graphregion(color(white)) ysize(3) 	///
    xsize(3.75) name(contour_min)
graph export contourPMC_CNA_min.pdf, name(contour_min) replace 
restore
collapse (max) hr_cv, by(cna logPM)
    label variable cna "CNA fraction"
    label variable hr "Maximum cross-validated mortal hazard ratio"

twoway contour hr_cv cna logPM,   ccuts(1(1.5)19.5)  			///
    ccolors(blue "0 59 255" "0 120 255" cyan*1.4 cyan cyan*.6 green*1.8 ///
        green*1.4 green  )	                			///
    xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 				///
        `10' "10"  `30' "30"  `60' "60"					///
        `100' "100"  `300' "300"  `600' "600"				///
        `1000' "1000"  `3000' "3000"  `6000' "6000" 			///
        `10000' "10000" `20000' "20000" , angle(45))	       		///
    xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 			///
        `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 	///
        `200' "200" `400' "400" `500' "500" `700' "700" `800' "800" 	///
        `900' "900" `2000' "2000" `4000' "4000" `5000' "5000" `7000' 	///
        "7000" `8000' "8000" `9000' "9000") 				///
    note("Unadjusted,  " 						///
        "CNA and log PMC modeled with 4 and 3 knots, respectively.  " 	///
        "Maximum hazard from 5 cross-validation models  " 		///
        " ", ring(0) position(4)) 					///
    scolor(blue) ecolor(red) xtitle( Point mutation count) 		///
    ylabel(0(.1)1, angle(0)) graphregion(color(white)) ysize(3) 	///
    xsize(3.75) name(contour_max)
graph export contourPMC_CNA_max.pdf, name(contour_max) replace 

log close
translate overall_surv_rqs_quartiles_crossValidated.log 		///
          overall_surv_rqs_quartiles_crossValidated.pdf, replace
