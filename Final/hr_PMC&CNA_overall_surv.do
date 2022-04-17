log using hr_PMC&CNA_overall_surv.log, replace
* hr_PMC&CNA_overall_surv.log

* Contour plot of mortal hazard ratio as a function of point mutation count (PMC) 
* and copy number aberration (CNA) fraction 

* This program generates Figures 3E and 3F 
pause on
use ../build.dta, clear
drop if OverallSurvivalMonths == . | overallfate == .
gen followup = OverallSurvivalMonths/12
* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
stset follow15, failure(fate15) scale(1) 
di "sample size = "_N
preserve
keep if PointMutationCount >3000 & CNAFractionGenomeAltered > .7
list PointMutationCount  CNAFractionGenomeAltered
*pause
restore

* Regress survival against log PMC and CNA using a 3 and 4 knots, respectively
gen logPM = log(PointMutationCount)
mkspline _Spm = logPM, cubic nknots(3) displayknots
* mkspline generates default knot locations from the real data. We need
* to save these locations as local macros so that we can generate
* log hazard ratios for the lattice of values needed for the contour plot
matrix knotspm = r(knots)
 local ktpm1 = knotspm[1,1]
 local ktpm2 = knotspm[1,2]
 local ktpm3 = knotspm[1,3]
 

mkspline _Scna = CNAFractionGenomeAltered, cubic nknots(4) displayknots
matrix knotscna = r(knots)
 local ktcna1 = knotscna[1,1]
 local ktcna2 = knotscna[1,2]
 local ktcna3 = knotscna[1,3]
 local ktcna4 = knotscna[1,4]

*pause
di "sample size = "_N

stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
test _Scna1 _Scna2 _Scna3
sum logPM CNAFractionGenomeAltered

* Let's look at the Wald test
test _Scna1 _Scna2 _Scna3 _Spm1 _Spm2

estimates stats
predict loghr, xb
sum loghr
gen hr = exp(loghr)
sum hr

preserve
di "sample size = "_N


* estimate HR and 95% CI when PMC = 100 & CNAF = .4
* estimate HR and 95% CI when PMC = 50  & CNAF = .3
*                             PMC = 300 & CNAF = .3
*                             PMC = 50  & CNAF = .6
*                             PMC = 300 & CNAF = .6
* In general, we do not know that any patient has these values Let's not depend on this.

local Np4 = _N+5
set obs `Np4'
replace PointMutationCount = 100      if _n ==`Np4' -4
replace CNAFractionGenomeAltered = .4 if _n ==`Np4' -4
replace PointMutationCount = 50       if _n ==`Np4' -3
replace CNAFractionGenomeAltered = .3 if _n ==`Np4' -3
replace PointMutationCount = 300      if _n ==`Np4' -2
replace CNAFractionGenomeAltered = .3 if _n ==`Np4' -2
replace PointMutationCount = 50       if _n ==`Np4' -1
replace CNAFractionGenomeAltered = .6 if _n ==`Np4' -1
replace PointMutationCount = 300      if _n ==`Np4' 
replace CNAFractionGenomeAltered = .6 if _n ==`Np4' 

sum hr
drop logPM loghr hr  _S*
gen logPM = log(PointMutationCount)
mkspline _Spm = logPM, cubic nknots(3) displayknots
mkspline _Scna = CNAFractionGenomeAltered, cubic nknots(4) displayknots

predict loghr, xb
sum loghr
gen hr = exp(loghr)
sum hr
predict se, stdp
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
di "HR when PMC = " PointMutationCount[`Np4'-4] ", & CNAF = " CNAFractionGenomeAltered[`Np4'-4] ///
    " is " hr[`Np4'-4] " 95% CI = " lb[`Np4'-4] ", " ub[`Np4'-4]
di "HR when PMC = " PointMutationCount[`Np4'-3] ", & CNAF = " CNAFractionGenomeAltered[`Np4'-3] ///
    " is " hr[`Np4'-3] " 95% CI = " lb[`Np4'-3] ", " ub[`Np4'-3]
di "HR when PMC = " PointMutationCount[`Np4'-2] ", & CNAF = " CNAFractionGenomeAltered[`Np4'-2] ///
    " is " hr[`Np4'-2] " 95% CI = " lb[`Np4'-2] ", " ub[`Np4'-2]
di "HR when PMC = " PointMutationCount[`Np4'-1] ", & CNAF = " CNAFractionGenomeAltered[`Np4'-1] ///
    " is " hr[`Np4'-1] " 95% CI = " lb[`Np4'-1] ", " ub[`Np4'-1]
di "HR when PMC = " PointMutationCount[`Np4'] ", & CNAF = " CNAFractionGenomeAltered[`Np4'] ///
    " is " hr[`Np4']   " 95% CI = " lb[`Np4'] ",   " ub[`Np4']
*pause
restore
preserve
* calculate a lattice of values of logPM and CNAFractionGenomeAltered for drawing 
* a contour plot of hr from this model

local delta = 0.01
sum logPM
local maxlogPM = r(max)
local inclogPM = `maxlogPM'*`delta'
local maxcna = 1
local inccna =`delta'
clear
local col = round(1/`delta')+1
local N = `col'^2
set obs `N'
di "_N = " _N
gen cna = int((_n-1)/`col')/(`col'-1)
gen logPM = mod((_n-1),`col')*`inclogPM'
pause
frame copy default grapharray
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp')
    }
}
* N.B.  We need to use the default knots generated from the real data
mkspline _Spm = logPM, cubic knots(`ktpm1'  `ktpm2' `ktpm3') displayknots
mkspline _Scna = cna, cubic knots(`ktcna1' `ktcna2' `ktcna3' `ktcna4') displayknots

predict loghr, xb
sum loghr
gen hr = exp(loghr)
sum hr

label variable cna "CNA fraction"
label variable hr "Mortal hazard ratio"

twoway contour hr cna logPM,   ccuts(1(1.5)19.5)  		///
    ccolors(blue "0 59 255" "0 120 255" cyan*1.4 cyan cyan*.6 green*1.8 green*1.4 ///
        green  )	                ///
    xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
        `10' "10"  `30' "30"  `60' "60"				///
        `100' "100"  `300' "300"  `600' "600"			///
        `1000' "1000"  `3000' "3000"  `6000' "6000" 		///
        `10000' "10000" `20000' "20000" , angle(45))	       			///
    xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
        `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
        `200' "200" `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
        `2000' "2000" `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
    subtitle("Unadjusted" "CNA fraction and log SM count modeled with 4 and 3 knots, respectively" ///
        , ring(1) position(12) size(medsmall)) scolor(blue) ecolor(red)			///
    xtitle(SM count) ylabel(0(.1)1, angle(0)) 			///
    graphregion(color(white)) ysize(3) xsize(3.75) name(contour1)
graph export contourPMC_CNA1.emf, name(contour1) replace 
*pause

* Adjust for cancer-type

restore
stcox _S* ib3.cancer
test _Spm1 _Spm2  _Scna1 _Scna2 _Scna3
di "P = " r(p)
test _Scna1 _Scna2 _Scna3
di "P = " r(p)
estimates stats
frame change grapharray
* N.B.  We need to use the default knots generated from the real data
mkspline _Spm = logPM, cubic knots(`ktpm1'  `ktpm2' `ktpm3') displayknots
mkspline _Scna = cna, cubic knots(`ktcna1' `ktcna2' `ktcna3' `ktcna4') displayknots

*drop loghr hr cna logpmr
gen loghr = _Spm1*_b[_Spm1] + _Spm2*_b[_Spm2] + _Scna1*_b[_Scna1] + _Scna2*_b[_Scna2]+ _Scna3*_b[_Scna3] 
gen hr = exp(loghr)
sum hr
*gen cna = round( CNAFractionGenomeAltered, .003)
*gen logpmr = round(logPM,.03)
label variable cna "CNA fraction"
label variable hr "Mortal hazard ratio"

* delete duplicate values
*sort cna logpmr
*drop if cna==cna[_n-1] & logpmr == logpmr[_n-1] 
*drop if hr ==. | cna ==. |logpmr==.

twoway contour hr cna logPM,   ccuts(1(.2)3.8)	ccolors(blue "0 59 255" 	///
    "0 120 255" cyan*1.4 cyan cyan*.6 green*1.8 green*1.4 green  )	        ///
    xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
        `10' "10"  `30' "30"  `60' "60"				///
        `100' "100"  `300' "300"  `600' "600"			///
        `1000' "1000"  `3000' "3000"  `6000' "6000" 		///
        `10000' "10000" `20000' "20000" , angle(45))	       			///
    xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
        `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
        `200' "200" `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
        `2000' "2000" `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
    subtitle("Adjusted for cancer-type" "CNA fraction and log SM count modeled with 4 and 3 knots, respectively" ///
        , ring(1) position(12)size( medsmall)) scolor(blue) ecolor(red)			///
     xtitle(SM count) ylabel(0(.1)1, angle(0)) 			///
     graphregion(color(white)) ysize(3) xsize(3.75) name(contour2)
     
graph export contourPMC_CNA2.emf, name(contour2) replace 
frame change default
* N.B. To maximize graph resolution, save the graph as an .emf file, open it in
* Snagit, and then copy and paste it into Photoshop.
log close
translate hr_PMC&CNA_overall_surv.log hr_PMC&CNA_overall_surv.pdf, replace

