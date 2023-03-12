log using hr_CNA_overall_surv.log, replace
* hr_CNA_overall_surv.log

* Initial analysis of the effects of the copy number aberration fraction score
* (CNA) on overall survival. Explore how model fit can be imporved with 
* restricted cubic splines

* This program generated Figures 3C and 3D 

use ../build.dta 
drop if OverallSurvivalMonths == . | overallfate == .
gen followup = OverallSurvivalMonths/12

* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
stset follow15, failure(fate15) scale(1)
     
* Regress survival against CNA using a 3 knot model
mkspline _S = CNAFractionGenomeAltered , cubic nknots(3) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
di "test of linearity p = " 2*normal(-abs(_b[_S2]/_se[_S2]))
estimates stats

* Regress survival against CNA using a 4-knot model
drop  _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(4) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
test _S2 _S3
di "test of linearity p = " r(p)
estimates stats

* Regress survival against CNA using a 5-knot model
drop  _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(5) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
test  _S2 _S3 _S4
di "test of linearity P = " r(p)
estimates stats

* N.B. The 5-knot model does not fit as well as the 4-knot model. We will use the
* 4-knot model going forward.

* Regress survival against CNA using a 4-knot model
drop  _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(4) displayknots
stcox _S*
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
test _S2 _S3
di "test of linearity p = " r(p)
estimates stats
* Test proportional hazards assumption
 estat phtest


predict loghr, xb
sum loghr
gen hr = exp(loghr)
predict se, stdp
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
sort CNAFractionGenomeAltered
list hr lb ub if CNAFractionGenomeAltered == .4

twoway (rarea lb ub CNAFractionGenomeAltered , fcolor(yellow) lcolor(yellow)		///
lwidth(none)) (line hr CNAFractionGenomeAltered , lcolor(red) lwidth(thick))     	///
    , ytitle(Mortal hazard relative to patients with CNA fraction = 0)      		///
    ylabel(1(.2)2.6, angle(zero)) xtitle(CNA fraction) xlabel(0(.1)1)     		///
    legend(off) 		///
     graphregion(color(white)) ysize(3) xsize(3.75) name(hrCNA)
graph export hr_graph_by_CNA.pdf, name(hrCNA) replace 

* Create horizontal line color coded by sextile
gen y = 1
gen x1 =     0 		if _n==1
replace x1=  0.033 	if _n == 2
gen x2 =     0.033	if _n==3
replace x2=  0.112	if _n == 4
gen x3 =     0.112	if _n==5
replace x3 = 0.194	if _n == 6
gen x4 =     0.194	if _n==7
replace x4 = 0.308	if _n == 8
gen x5 =     0.308	if _n==9
replace x5 = 0.480	if _n == 10
gen x6 =     0.480	if _n==11
replace x6 = 1.0	if _n == 12
twoway (rarea lb ub CNAFractionGenomeAltered , fcolor(yellow) lcolor(yellow)		///
lwidth(none)) (line hr CNAFractionGenomeAltered , lcolor(red) lwidth(thick))     	///
        (line y x1, lcolor(blue) lwidth(thick))						///
        (line y x2, lcolor(red) lwidth(thick))						///
        (line y x3, lcolor(green) lwidth(thick))					///
        (line y x4, lcolor(cyan) lwidth(thick))						///
        (line y x5, lcolor(magenta) lwidth(thick))					///
        (line y x6, lcolor(black) lwidth(thick))					///
    , ytitle(Mortal hazard relative to patients with CNA fraction = 0)      		///
    ylabel(1(.2)2.6, angle(zero)) xtitle(CNA fraction) xlabel(0(.1)1)     		///
    legend(off) 		///
     graphregion(color(white)) ysize(3) xsize(3.75) name(hrCNA2)


preserve

* estimate HR and 95% CI when CNAF = 0.1
* In general, we do not know that any patient has CNAF = 0.1. Let's not depend on this.

local Np1 = _N+1
set obs `Np1'
replace CNAFractionGenomeAltered = 0.1 if _n ==`Np1'
drop  loghr hr se lb ub _S*
mkspline _S = CNAFractionGenomeAltered , cubic nknots(4) displayknots
predict loghr, xb
sum loghr
gen hr = exp(loghr)
predict se, stdp
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
di "HR when CNAF = " CNAFractionGenomeAltered[`Np1'] " is " hr[`Np1'] " 95% CI = " lb[`Np1'] ", " ub[`Np1']

restore
drop loghr hr se lb ub 

* Regress survival against CNA and cancer type using a 4-knot model
stcox _S* ib3.cancer
test _S1 _S2 _S3
di "P = " r(p)
test  _S2 _S3
di "test of linearity P = " r(p)

estimates stats
* Test proportional hazards assumption
 estat phtest

predictnl loghr = _S1*_b[_S1] + _S2*_b[_S2]  + _S3*_b[_S3],se(se)
gen hr = exp(loghr)
gen lb = exp(loghr - 1.96*se)
gen ub = exp(loghr + 1.96*se)
twoway (rarea lb ub CNAFractionGenomeAltered , fcolor(yellow) lcolor(yellow) lwidth(none)) ///
    (line hr CNAFractionGenomeAltered , lcolor(red) lwidth(thick)) 			   ///
  , ytitle(Mortal hazard relative to patients with CNA fraction = 0) ylabel(1(.1)1.8, 	   ///
        angle(zero)) xtitle(CNA fraction) xlabel(0(.1)1) 				   ///
    legend(off) 		   ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(hrCNAandCA)
graph export hr_graph_by_CNAandCA.pdf, name(hrCNAandCA) replace





log close
translate hr_CNA_overall_surv.log hr_CNA_overall_surv.pdf, replace
