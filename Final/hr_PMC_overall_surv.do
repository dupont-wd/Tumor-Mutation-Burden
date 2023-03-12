log using hr_PMC_overall_surv.log,  replace
* hr_PMC_overall_surv.log 

* Initial analysis of the effects of the point mutation count (PMC), 
* on hazard ratios for overall survival.
* Explore how model fit can be imporved with restricted cubic splines

* This program generated Figures 3A and 3B

use ../build.dta, clear
drop if OverallSurvivalMonths == . | overallfate == .
gen followup = OverallSurvivalMonths/12

* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
stset follow15, failure(fate15) scale(1)

* Regress survival against log PMC using a 3 knot model
gen logPM = log(PointMutationCount)

mkspline _S = logPM, cubic nknots(3) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
di "test of linearity p = " 2*normal(-abs(_b[_S2]/_se[_S2]))
estimates stats

* Try 4-knot models

* Regress survival against PMC using a 4-knot model
drop _S*
mkspline _S = logPM, cubic nknots(4) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
test _S2 _S3
di "test of linearity p = " r(p)

estimates stats

* N.B. The 4-knot model does not fit as well as the 3-knot model. We will use the
* 3-knot model going forward.

drop _S*
mkspline _S = logPM, cubic nknots(3) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(2, e(chi2))
di "test of linearity p = " 2*normal(-abs(_b[_S2]/_se[_S2]))
estimates stats
* Test proportional hazards assumption
 estat phtest


predict loghr, xb
sum loghr
gen hr = exp(loghr)
sort logPM
predict se, stdp
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
list hr lb ub if PointMutationCount==100 
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp')
    }
}
twoway (rarea lb ub logPM, fcolor(yellow) lcolor(yellow) 			///
        lwidth(none)) (line hr logPM, lcolor(red) lwidth(thick)) 		///
    , ytitle(Mortal hazard relative to patients with SM count = 1)			///
      ylabel(0(5)30, angle(zero)) xtitle(SM count)      		///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45))	       			///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      subtitle("  Unadjusted", position(10) ring(0)) 	///
      legend(order(1  "95% CI" 2 "Expected") cols(1) position(1) ring(0)) 	///
      graphregion(color(white)) ysize(3) xsize(3.75) name(hrPMC)
graph export hr_graph_by_PMC.pdf, name(hrPMC) replace 

* List hr values associated with the median SM counts for each sextile
* (See overall_surv_6tiles.log)
preserve
sort PointMutationCount
by   PointMutationCount: drop if _n > 1
list PointMutationCount logPM hr if PointMutationCount == 12 | PointMutationCount== 28 ///
                                  | PointMutationCount== 47  | PointMutationCount== 75 ///
                                  | PointMutationCount== 134  | PointMutationCount== 435
restore
* Create horizontal line color coded by sextile
gen y = 0
gen x1 =     log(1) 	if _n==1
replace x1=  log(20) 	if _n == 2
gen x2 =     log(21)	if _n==3
replace x2=  log(36)	if _n == 4
gen x3 =     log(37)	if _n==5
replace x3 = log(59)	if _n == 6
gen x4 =     log(60)	if _n==7
replace x4 = log(97)	if _n == 8
gen x5 =     log(98)	if _n==9
replace x5 = log(211)	if _n == 10
gen x6 =     log(212)	if _n==11
replace x6 = log(25711)	if _n == 12
twoway (rarea lb ub logPM, fcolor(yellow) lcolor(yellow) 			///
        lwidth(none)) (line hr logPM, lcolor(red) lwidth(thick)) 		///
        (line y x1, lcolor(blue) lwidth(thick))				///
        (line y x2, lcolor(red) lwidth(thick))					///
        (line y x3, lcolor(green) lwidth(thick))				///
        (line y x4, lcolor(cyan) lwidth(thick))				///
        (line y x5, lcolor(magenta) lwidth(thick))				///
        (line y x6, lcolor(black) lwidth(thick))				///
    , ytitle(Mortal hazard relative to patients with SM count = 1)		///
      subtitle("  Unadjusted", position(10) ring(0)) 	///
      ylabel(0(5)30, angle(zero)) xtitle(SM count)      			///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45))	       			///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      legend(order(1  "95% CI" 2 "Expected") cols(1) position(1) ring(0)) 	///
      graphregion(color(white)) ysize(3) xsize(3.75) name(hrPMC2)

preserve

* estimate HR and 95% CI when PMC = 10
* In general, we do not know that any patient has PMC = 10. Let's not depend on this.

local Np1 = _N+1
set obs `Np1'
replace PointMutationCount = 10 if _n ==`Np1'
drop logPM loghr hr se lb ub _S*
gen logPM = log(PointMutationCount) 
mkspline _S = logPM, cubic nknots(3) displayknots
predict loghr, xb
sum loghr
gen hr = exp(loghr)
predict se, stdp
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
di "HR when PMC = " PointMutationCount[`Np1'] " is " hr[`Np1'] " 95% CI = " lb[`Np1'] ", " ub[`Np1']

restore
* Regress survival against cancer and PMC with a 3-knot model
drop loghr hr se lb ub
stcox _S* ib3.cancer
test _S1 _S2
di "P = " r(p)
estimates stats
* Test proportional hazards assumption
 estat phtest

predictnl loghr = _S1*_b[_S1] + _S2*_b[_S2] ,se(se)
gen hr = exp(loghr)
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)

twoway (rarea lb ub logPM, fcolor(yellow) lcolor(yellow) 			///
        lwidth(none)) (line hr logPM, lcolor(red) lwidth(thick)) 		///
    , ytitle(Mortal hazard relative to patients with SM count = 1)			///
      ylabel(0(1)5, angle(zero)) xtitle(SM count)      		///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45))	       			///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      subtitle("  Adjusted for" "  cancer-type", position(10) ///
          ring(0)) legend(off) graphregion(color(white)) ysize(3) xsize(3.75) name(hrPMCandCA)
graph export hr_graph_by_PMCandCA.pdf, name(hrPMCandCA) replace

twoway (rarea lb ub logPM, fcolor(yellow) lcolor(yellow) 			///
        lwidth(none)) (line hr logPM, lcolor(red) lwidth(thick)) 		///
        (line y x1, lcolor(blue) lwidth(thick))				///
        (line y x2, lcolor(red) lwidth(thick))					///
        (line y x3, lcolor(green) lwidth(thick))				///
        (line y x4, lcolor(cyan) lwidth(thick))				///
        (line y x5, lcolor(magenta) lwidth(thick))				///
        (line y x6, lcolor(black) lwidth(thick))				///
    , ytitle(Mortal hazard relative to patients with SM count = 1)			///
      ylabel(0(1)5, angle(zero)) xtitle(SM count)      		///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45))	       			///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      subtitle("  Adjusted for" "  cancer-type", position(10) ///
          ring(0)) legend(off) graphregion(color(white)) ysize(3) xsize(3.75) name(hrPMCandCA2)
graph export hr_graph_by_PMCandCA2.pdf, name(hrPMCandCA2) replace
log close
translate hr_PMC_overall_surv.log hr_PMC_overall_surv.pdf, replace

