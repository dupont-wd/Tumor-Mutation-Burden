log using overall_surv_rqs_quartiles.log, replace
* overall_surv_rqs_quartiles.log
* Generate Figure 4A

* Explore the combined effects of a 3-knot model of log PMC 
* and a 4-knot model of CNA on survival. Also, look at fusion pairs
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

stcox _S*
estimates stats
* Determine the best RQS model for FusionPairs
* Note that attempts to define 4 and 5 knot splines for FusionPairs failed

stcox FusionPairs
mkspline _S = FusionPairs, cubic nknots(3) displayknots
stcox _S1 _S2
estimates stats
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))

stcox _S1 _S2 i.cancer
test _S1 _S2
di "P-value for fusion pairs adjusted for cancer = "r(p)

drop _S*
*mkspline _S = FusionPairs, cubic nknots(4) displayknots
*stcox _S1 _S2 _S3
* Analysis with 4 knots fails "Sample size too small for this many knots."

* Combine point mutations and CNA fraction and fusion pairs in the same model
mkspline _Spmc = logPM, cubic nknots(3) displayknots
mkspline _Scna = CNAFractionGenomeAltered, cubic nknots(4) displayknots
mkspline _Sfpc = FusionPairs, cubic nknots(3) displayknots

stcox _Spmc* _Scna* FusionPairs
estimates stats
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))
test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 FusionPairs
di "P-value for PMC, CNA and fusion pairs = "r(p)

test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 
di "P-value for PMC & CNA adjusted for fusion pairs = "r(p)

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
sts test quartile if quartile <3, logrank
sts test quartile if quartile != 1 & quartile < 4, logrank
sts test quartile if quartile > 2 & quartile != ., logrank

* Draw a survival plot associated with log PMC, CNA and fusion pairs
sts graph, by(quartile) ci risktable risktable(, rowtitle(1st quartile) 	///
        size(small) group(#1)) 							///
    risktable(, rowtitle(2nd quartile) size(small) group(#2)) 			///
    risktable(, rowtitle(3rd quartile) size(small) group(#3))			///
    risktable(, rowtitle(4th quartile) size(small) group(#4)) 			///
    plot1opts(lcolor(blue)) plot2opts(lcolor(red)) 				///
    plot3opts(lcolor(forest_green)) plot4opts(lcolor(gold))  			///
    ci1opts(fcolor(blue%40)) ci2opts(fcolor(red%40) lwidth(none))		///
    ci3opts(fcolor(forest_green%40)) ci4opts(fcolor(gold%40)) 			///
    ytitle(Probability of overall survival) 	 				///
    ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
    title("", size(zero)) 							///
    note("Shaded regions denote 95% CIs." 					///
        "Survival associated with log PMC (3-knots), CNA (4-knots) & fusion pairs"	///
        ,size(small) position(2) ring(0)) 					///
    legend(title(Quartiles) order(9 "1st" 10 "2nd" 11 "3rd" 12 "4th") 		///
    symxsize(*.75) textwidth(*.75) cols(1) position(8) ring(0)) name(quartile_rqs1)
*graph export quartile_graph_rqs1.pdf, name(quartile_rqs1) replace     

* Fusion pairs are not significant when included in a model with PMC and CNAF
* Let's repeat the above using a 3-knot model for FPC

drop loghr loghr4 quartile
stcox _Spmc* _Scna* _Sfpc*
estimates stats
di "Chi^2 = " e(chi2) " p = " chi2tail(e(df_m), e(chi2))

test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 _Sfpc1 _Sfpc2
di "P-value for PMC, CNA and fusion pairs = "r(p)

test  _Sfpc1 _Sfpc2
di "P-value for FPC adjusted for PMC, CNA  = "r(p)

test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 
di "P-value for PMC & CNA adjusted for fusion pairs = "r(p)


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
sts test quartile if quartile <3, logrank
sts test quartile if quartile != 1 & quartile < 4, logrank
sts test quartile if quartile > 2 & quartile != ., logrank

* Draw a survival plot associated with log PMC, CNA and fusion pairs
sts graph, by(quartile) ci risktable risktable(, rowtitle(1st quartile) 	///
        size(small) group(#1)) 							///
    risktable(, rowtitle(2nd quartile) size(small) group(#2)) 			///
    risktable(, rowtitle(3rd quartile) size(small) group(#3))			///
    risktable(, rowtitle(4th quartile) size(small) group(#4)) 			///
    plot1opts(lcolor(blue)) plot2opts(lcolor(red)) 				///
    plot3opts(lcolor(forest_green)) plot4opts(lcolor(gold))  			///
    ci1opts(fcolor(blue%40)) ci2opts(fcolor(red%40) lwidth(none))		///
    ci3opts(fcolor(forest_green%40)) ci4opts(fcolor(gold%40)) 			///
    ytitle(Probability of overall survival) 	 				///
    ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
    title("", size(zero)) 							///
    note("Shaded regions denote 95% CIs." 					///
        "Survival associated with log PMC (3-knots), CNA (4-knots) & FPC (3-knots)"	///
        ,size(small) position(2) ring(0)) 					///
    legend(title(Quartiles) order(9 "1st" 10 "2nd" 11 "3rd" 12 "4th") 		///
    symxsize(*.75) textwidth(*.75) cols(1) position(8) ring(0)) name(quartile_rqs3)



* Let's drop fusion pairs from the model


            
drop loghr loghr4 quartile
stcox _Spmc* _Scna* 
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
    risktable(, rowtitle(1st quartile) 	group(#1) size(medium)) 		///
    risktable(, rowtitle(2nd quartile)  group(#2) size(medium))			///
    risktable(, rowtitle(3rd quartile)  group(#3) size(medium))			///
    risktable(, rowtitle(4th quartile)  group(#4) size(medium))			///
    plot1opts(lcolor(blue) lwidth(thick) ) 					///
    plot2opts(lcolor(red) lwidth(thick))					///
    plot3opts(lcolor(green)lwidth(thick)) 					///
    plot4opts(lcolor(black) lwidth(thick))  					///
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30)  lwidth(none)) 		///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			///
    ytitle(Probability of overall survival) 	 				///
    ylabel(0(.2)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
    title("", size(zero)) 							///
    subtitle("Shaded regions denote 95% CIs." 					///
        "CNA count and log SM count modeled" "with 4 and 3 knots, respectivley", ///
        size(medium) position(2) ring(0)) 					///
    legend(subtitle("Quartiles") order(9 "1st" 10 "2nd" 11 "3rd" 12 "4th") 	///	
        cols(1) position(8) ring(0))  graphregion(color(white)) 		///
    xsize(3.25) ysize(2.6) name(quartile_rqs2)
graph export quartile_graph_CNA&PMC.pdf, name(quartile_rqs2) replace  

* Repeat the preceding analysis adjusted for cancer type

stcox _Spmc* _Scna* i.cancer
estimates stats

test _Spmc1 _Spmc2  _Scna1 _Scna2 _Scna3 
di "P-value for PMC and CNA = "r(p)


log close
translate overall_surv_rqs_quartiles.log overall_surv_rqs_quartiles.pdf, replace
