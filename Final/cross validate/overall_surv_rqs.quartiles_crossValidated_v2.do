log using overall_surv_rqs_quartiles_crossValidated_v2.log, replace
* overall_surv_rqs_quartiles_crossValidated_v2.log

* Perform a cross-validation analysis of the TMB model, refitting the 
* model from scratch for each of the 5 cross-validation data sets

* First, regenerate the TMB model on the complete data set
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
gen hr_unvalidated = exp(loghr)
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
    subtitle("CNA fraction and log SM count modeled with 4 and 3 knots, respectivley" ///
        "Shaded regions denote 95% CIs", size(small) position(2) ring(0)) 	///
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
gen se = .
gen lb = .
gen ub = .
* Calculate reference line
gen refx = 1 if _n==1
gen refy=1 if _n== 1
replace refx = 20 if _n == 2
replace refy = 20 if _n == 2

 
forvalues j = 1/5 {
*!!!!!!!!!! Start cross-validation analysis !!!!!!!!!!!!

    di " ""!!!!!!!!!! start of group `j' "
    di "Generate optimal model based on all of the data except patients in group `j'"
    drop _S* loghr loghr_tmp hr se lb ub
    * First, fit log PMC (log SMC) to the survival data
    
    * Regress survival against logPM -- linear model
    stcox logPM if group != `j'
    estimates stats
    matrix STAT = r(S)
    matrix list STAT
    local old_aic =STAT[1,5]
    di "AIC from model regressing survival against logPM = `old_aic'"
    
    * Find optimal cubic spline model for the effect of logPM on survival
    forvalues i = 3/7 {
        di "Regress survival against logPM using `i' knots"
        mkspline _Slpm = logPM if group != `j', cubic nknots(`i') displayknots
        stcox _Slpm* if group != `j' 
        di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
        estimates stats
        matrix STAT = r(S)
        matrix list STAT
        local new_aic =STAT[1,5]
        di "AIC from model regressing survival against logPM = `new_aic' using `i' knots"
        if `new_aic' > `old_aic' {
            local final_i = `i' - 1
            local n_coef_lpm = `final_i' -1
            continue, break
        }
        local old_aic =`new_aic'
        drop _Slpm*
    }
    * Regress logPM against survival using the optimal number of knots
    drop _Slpm*
    mkspline _Slpm = logPM if group != `j', cubic nknots(`final_i') displayknots
    matrix lpm_knots = r(knots)
    matrix list lpm_knots

    stcox _Slpm*   if group != `j'
    di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
    estimates stats
        
    
    predict loghr if group != `j'  , xb
    sum loghr
    gen hr = exp(loghr)
    sort logPM
    predict se if group != `j', stdp
    gen lb = exp(loghr - 1.96*se)
    gen ub = exp(loghr + 1.96*se)
    list hr lb ub if PointMutationCount==100 
    foreach i of numlist  1 10 100 1000 10000 {
        forvalues jj = 1/9 {
            local tmp = `i'*`jj'
            local `tmp'  = log(`tmp')
        }
    }
    twoway (rarea lb ub logPM, fcolor(yellow) lcolor(yellow) 			///
            lwidth(none)) (line hr logPM, lcolor(red) lwidth(medium)) 		///
        , ytitle(Mortal hazard relative to patients with SM count = 1)		///
          ylabel(0(5)30, angle(zero)) xtitle(SM count)      	///
          xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 				///
              `10' "10" `20' "20" `30' "30"  `60' "60"				///
              `100' "100"  `200' "200" `300' "300"  `600' "600"			///
              `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
              `10000' "10000" `20000' "20000" , angle(45))	       		///
          xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 			///
              `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 		///
              `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 	///
              `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
          subtitle("`final_i'-knot model of log SM count" 		///
              "Cross-validation sample # `j'",   position(5) ring(0)) 		///
          legend(order(1  "95% CI" 2 "Expected") cols(1) position(1) ring(0)) 	///
          graphregion(color(white)) ysize(3) xsize(3.75) name(hrPMC`j')
    graph export hr_graph_by_PMC_v2_cv`j'.pdf, name(hrPMC`j') replace 
    
    * Second, fit CNAF to survival data
    
    drop loghr hr se lb ub
    
    * Regress survival against CNAFractionGenomeAltered --linear model
    stcox CNAFractionGenomeAltered if group != `j'   
    estimates stats
    matrix STAT = r(S)
    matrix list STAT
    local old_aic =STAT[1,5]
    di "AIC from model regressing survival against CNAFractionGenomeAltered = `old_aic'"
    
    * Find optimal cubic spline model for the effect of CNAFractionGenomeAltered on survival
    forvalues i = 3/7 {
        di "Regress survival against CNAFractionGenomeAltered using `i' knots"
        mkspline _Scna = CNAFractionGenomeAltered if group != `j', cubic nknots(`i') displayknots
        stcox _Scna*  if group != `j' 
        di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
        estimates stats
        matrix STAT = r(S)
        matrix list STAT
        local new_aic =STAT[1,5]
        di "AIC from model regressing survival against CNAFractionGenomeAltered = `new_aic' using `i' knots"
        if `new_aic' > `old_aic' {
            local final_i = `i' - 1
            local n_coef_cna = `final_i' -1
            continue, break
        }
        local old_aic =`new_aic'
        drop _Scna*
    }
    * Regress CNAFractionGenomeAltered against survival using the optimal number of knots
    drop _Scna*
    mkspline _Scna = CNAFractionGenomeAltered if group != `j', cubic nknots(`final_i') displayknots
    matrix cna_knots = r(knots)
    matrix list cna_knots
    
    stcox _Scna*   if group != `j'
    di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
    estimates stats
            
    predict loghr  if group != `j', xb
    sum loghr
    gen hr = exp(loghr)
    sort CNAFractionGenomeAltered
    predict se if group != `j', stdp
    gen lb = exp(loghr - 1.96*se)
    gen ub = exp(loghr + 1.96*se)
    list hr lb ub if CNAFractionGenomeAltered == .4
    
    twoway (rarea lb ub CNAFractionGenomeAltered , fcolor(yellow) lcolor(yellow)	///
    lwidth(none)) (line hr CNAFractionGenomeAltered , lcolor(red) lwidth(medium))     	///
        , ytitle(Mortal hazard relative to patients with CNA fraction = 0)      	///
        ylabel(1(.2)2.6, angle(zero)) xtitle(CNA fraction) xlabel(0(.1)1)     		///
        subtitle("`final_i'-knot model of CNA fraction"					///
            "Cross-validation sample # `j'", position(4) ring(0)) 			///
        legend(order(1  "95% CI" 2 "Expected") cols(1) position(1) ring(0)) 		///
         graphregion(color(white)) ysize(3) xsize(3.75) name(hrCNA`j')
    graph export hr_graph_by_CNA_v2_cv`j'.pdf, name(hrCNA`j') replace 
    
    * Third, fit fusion pair data to the survival data
    
    drop loghr hr se lb ub
    
    * Regress survival against FusionPairs --linear model
    stcox FusionPairs  if group != `j'
    estimates stats
    matrix STAT = r(S)
    matrix list STAT
    local old_aic =STAT[1,5]
    di "AIC from model regressing survival against FusionPairs = `old_aic'"
    
    * Find optimal cubic spline model for the effect of FusionPairs on survival
    forvalues i = 3/7 {
        di "Regress survival against FusionPairs using `i' knots"
        capture mkspline _Sfp = FusionPairs if group != `j', cubic nknots(`i') displayknots

        if  _rc == 0 {
            drop _Sfp*
            mkspline _Sfp = FusionPairs if group != `j', cubic nknots(`i') displayknots
        }
        else {
            local final_i = `i' - 1
            local n_coef_fp = `final_i' -1
            continue, break
        }        
        stcox _Sfp* if group != `j'
        di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
        estimates stats
        matrix STAT = r(S)
        matrix list STAT
        local new_aic =STAT[1,5]
        di "AIC from model regressing survival against FusionPairs = `new_aic' using `i' knots"
        if `new_aic' > `old_aic' {
            local final_i = `i' - 1
            continue, break
        }
        local old_aic =`new_aic'
        drop _Sfp*
    }
    * Regress FusionPairs against survival using the optimal number of knots
    capture drop _Sfp*
    mkspline _Sfp = FusionPairs if group != `j', cubic nknots(`final_i') displayknots
    matrix fp_knots = r(knots)
    matrix list fp_knots

    stcox _Sfp* if group != `j'
    di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
    estimates stats
        
    predict loghr if group != `j', xb
    sum loghr
    gen hr = exp(loghr)
    sort FusionPairs
    predict se if group != `j', stdp
    gen lb = exp(loghr - 1.96*se)
    gen ub = exp(loghr + 1.96*se)
    
    twoway (rarea lb ub FusionPairs , fcolor(yellow) lcolor(yellow)			///
    lwidth(none)) (line hr FusionPairs , lcolor(red) lwidth(medium))     		///
        , ytitle(Mortal hazard relative to patients with FP count = 0)      		///
        ylabel(1(.2)2.6, angle(zero)) xtitle(FP count) xlabel(0(10)60)     		///
        subtitle("`final_i'-knot model of FP count"					///
            "Cross-validation sample # `j'", position(4) ring(0))  			///
        legend(order(1  "95% CI" 2 "Expected") cols(1) position(1) ring(0)) 		///
         graphregion(color(white)) ysize(3) xsize(3.75) name(hrFP`j')
    graph export hr_graph_by_FP_v2_cv`j'.pdf, name(hrFP`j') replace 
    
    di " "
    di "Combine logPM, CNAF and FPC in a multiplicative model "
    
    stcox _S* if group != `j'
    local tvalue _Slpm1
    forvalues i = 2/`n_coef_lpm' {
        local nextcoef _Slpm`i'
        local tvalue `tvalue' `nextcoef'
    }    
    test `tvalue'
    
    local no_lpm = 0
    if r(p) >0.05 {
        drop _Slpm*
        local no_lpm =1
    }  
    
    local tvalue _Scna1
    forvalues i = 2/`n_coef_cna' {
        local nextcoef _Scna`i'
        local tvalue `tvalue' `nextcoef'
   }    
    local no_cna = 0
    test `tvalue'
    if r(p) >0.05 {
        drop _Scna*
        local no_cna = 1
    }
    
    local tvalue _Sfp1
    forvalues i = 2/`n_coef_fp' {
        local nextcoef _Sfp`i'
        local tvalue `tvalue' `nextcoef'
    }    
    
    local no_fp = 0
    test `tvalue'
    if r(p) >0.05 {
        drop _Sfp*
        local no_fp = 1
    }

di "Run the final model for cross-validated dataset `j' !!!!!!!!!!!!!"

    stcox _S* if group != `j'
*   Make a matrix of the parameter estimates
    matrix PAR`j' = e(b)
    
*   List parameter estimates for the current group
    di "Parameter estimates for group `j' are as follows"
    matrix list PAR`j'
    drop _S*
    if `no_lpm' == 0 {
        local cols = colsof(lpm_knots)
        local knots = lpm_knots[1,1]
        local knots `knots'
        di "`knots'"
        forvalues i = 2/`cols' {
	   local next =lpm_knots[1,`i']
	   local next `next'
	   local knots `knots' `next'
	   di "`knots'"
        }
        mkspline _Slpm = logPM , cubic knots(`knots') displayknots
    }

    if `no_cna' == 0 {
        local cols = colsof(cna_knots)
        local knots = cna_knots[1,1]
        local knots `knots'
        di "`knots'"
        forvalues i = 2/`cols' {
	   local next =cna_knots[1,`i']
	   local next `next'
	   local knots `knots' `next'
	   di "`knots'"
            }
            mkspline _Scna = CNAFractionGenomeAltered, cubic knots(`knots') displayknots
        }
        if `no_fp' == 0 {
            di "Fusion pairs are in the model"
        }
             
    di  "Predict log hazard for patients in group `j' using model "
    di  "generated from all other patients"
    predict loghr_tmp if group ==`j', xb
    replace loghr_cv = loghr_tmp if group == `j'
    replace hr_cv = exp(loghr_tmp) if group == `j'
    codebook hr_cv
    
    di "Start of generating grid for contour plot from model # `j'
    frame create grid`j'
    frame change grid`j'
    local col = 101
    local N= `col'^2
    set obs `N'
    local inclogPM =0.1
    gen cna = int((_n-1)/`col')/(`col'-1)
    gen logPM = mod((_n-1),`col')*`inclogPM'
        foreach jj of numlist  1 10 100 1000 10000 {
        forvalues k = 1/9 {
            local tmp = `jj'*`k'
            local `tmp'  = log(`tmp')
        }
    }
*   The following will not work if the model includes FPC or excludes
*   any of the other TMB variables. It should be robust to changes in 
*   the number of knots although it may not work if the linear
*   model is best

    if `no_lpm' == 0 {
        local cols_lpm = colsof(lpm_knots)
        local knots = lpm_knots[1,1]
        local knots `knots'
        di "`knots'"
        forvalues i = 2/`cols_lpm' {
	   local next =lpm_knots[1,`i']
	   local next `next'
	   local knots `knots' `next'
	   di "`knots'"
        }
        mkspline _Slpm = logPM , cubic knots(`knots') displayknots
    }

    if `no_cna' == 0 {
        local cols_cna = colsof(cna_knots)
        local knots = cna_knots[1,1]
        local knots `knots'
        di "`knots'"
        forvalues i = 2/`cols_cna' {
	   local next =cna_knots[1,`i']
	   local next `next'
	   local knots `knots' `next'
	   di "`knots'"
        }
        mkspline _Scna = cna, cubic knots(`knots') displayknots
    }
    if `no_fp' == 0 {
        di "Fusion pairs are in the model"
    }
      
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
        note("CNA fraction and log SM count modeled with `cols_cna' and `cols_lpm' knots, respectively.  " ///
            "Model parameters from group `j'  " , ring(1) position(12)) ///
        scolor(blue) ecolor(red) xtitle(SM count) 		///
        ylabel(0(.1)1, angle(0)) graphregion(color(white)) ysize(3) 	///
        xsize(3.75) name(contour`j')
    graph export contourPMC_CNA`j'.pdf, name(contour`j') replace 
    frame change default
}
gen diff = hr_cv - hr_unvalidated
save cross_validated.dta, replace
twoway (scatter diff hr_unvalidated, mcolor(blue) msymbol(smcircle_hollow) )	///
  , yline(0) xtitle(unvalidated hazard ratio) name(diff)		///
    ytitle( "cross-validated - unvalidated hazard ratios") 
graph export cross_validated_differences.pdf, name(diff) replace 
    

regress hr_cv hr_unvalidated
predict yhat, xb
sort hr
twoway (scatter hr_cv hr_unvalidated, msymbol(smcircle_hollow) mcolor(blue)) 	///
    (line refy refx, lwidth(medthick) lcolor(red)) 			///
    (line yhat hr_unvalidated, lwidth(medthick) lcolor(cyan))		///	
  , xtitle(Unvalidated hazard ratio) 				///
    ytitle( Cross-validated hazard ratio)   			///
    legend(ring(0) position(10) subtitle(All patients) 		///
        order(1 "Observed" 2 "Reference" 3 "Expected") col(1)) 	///
    graphregion(color(white)) name(regress)
graph export cross_validated_regression.pdf, name(regress) replace 

    
* Regenerate the preceding scatter plots separately by group   
forvalues j = 1/5  {
    preserve
    keep if group == `j'
    * Calculate reference line
    replace refx = 1 if _n==1
    replace  refy=1 if _n== 1
    replace refx = 20 if _n == 2
    replace refy = 20 if _n == 2
    di "Group `j'"
     gen diff`j' = hr_cv - hr_unvalidated
    twoway (scatter diff`j' hr_unvalidated, mcolor(blue) msymbol(smcircle_hollow)) 	///
      , yline(0) xtitle(unvalidated hazard ratio) name(diff`j')	///
        ytitle( "cross-validated - unvalidated hazard ratios") 		///
        subtitle("Group `j'", ring(0) position(10) )
    graph export cross_validated_differences.pdf, name(diff`j') replace 

    regress hr_cv hr_unvalidated
    drop yhat
    predict yhat, xb
    sort hr_unvalidated
    twoway (scatter hr_cv hr_unvalidated, msymbol(smcircle_hollow) mcolor(blue)) 	///
        (line refy refx, lwidth(medthick) lcolor(red)) 			///
        (line yhat hr_unvalidated, lwidth(medthick) lcolor(cyan))			///	
      , xtitle(unvalidated hazard ratio) 				///
        ytitle( cross-validated hazard ratio)   			///
        legend(ring(0) position(10) subtitle(Group `j')			///
        order(1 "observed" 2 "reference" 3 "expected") col(1)) name(regress`j')
    graph export cross_validated_regression`j'.pdf, name(regress`j') replace 
    restore    
}
frame change grid1
forvalues j = 2/5 {
    frameappend grid`j'
}
preserve
collapse (min) hr_cv, by(cna logPM)
    label variable cna "CNA fraction"
    label variable hr_cv "Minimum cross-validated mortal hazard ratio"

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
    subtitle("CNA fraction and log SM count modeled with optimal number of knots" ///
        "Minimum hazard from 5 cross-validation models" 		///
        , ring(1) position(12) size(small)) 			///
    scolor(blue) ecolor(red) xtitle( SM count) 		///
    ylabel(0(.1)1, angle(0)) graphregion(color(white)) ysize(3) 	///
    xsize(3.75) name(contour_min)
graph export contourPMC_CNA_min.emf, name(contour_min) replace 
restore
collapse (max) hr_cv, by(cna logPM)
    label variable cna "CNA fraction"
    label variable hr_cv "Maximum cross-validated mortal hazard ratio"

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
    subtitle("CNA fraction and log SM count modeled with optimal number of knots " ///
        "Maximum hazard from 5 cross-validation models" 		///
        , ring(1) position(12) size(small)) 				///
    scolor(blue) ecolor(red) xtitle( SM count) 		///
    ylabel(0(.1)1, angle(0)) graphregion(color(white)) ysize(3) 	///
    xsize(3.75) name(contour_max)
graph export contourPMC_CNA_max.emf, name(contour_max) replace 

* use Snagit to import emf files into Photoshop with high resolution
log close
translate overall_surv_rqs_quartiles_crossValidated_v2.log 		///
          overall_surv_rqs_quartiles_crossValidated_v2.pdf, replace
