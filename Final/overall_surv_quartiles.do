log using overall_surv_quartiles.log, replace
* overall_surv_quartiles.log
*This program generates panels A, B and C of Figure 2 

* Draw survival plots of SM count, CNA fraction, and FG count by quartiles

use ../build.dta, clear
drop if OverallSurvivalMonths == . | overallfate == .
gen followup = OverallSurvivalMonths/12

* Censor survival at 15 years
gen follow15 = followup
replace follow15 = 15 if followup > 15
gen fate15 = overallfate
replace fate15= 0 if followup > 15
stset follow15, failure(fate15) scale(1)

* Explore effect of PMC on survival
centile PointMutationCount, centile(25(25)75)
forvalues i = 1/3 {
    local c`i' = r(c_`i')
}
sum PointMutationCount
local max = r(max)
gen pmc = recode( PointMutationCount, `c1' , `c2' , `c3'   , `max')
tabulate pmc, missing
sts test pmc
di "p = " chi2tail(`r(df)', `r(chi2)')

sts test pmc if pmc > `c1' & pmc != .
di "p = " chi2tail(`r(df)', `r(chi2)')

gen first_vs_rest = recode(pmc, `c1', `max')
tabulate first_vs_rest
sts test first_vs_rest
di "p = " chi2tail(`r(df)', `r(chi2)')

* Draw Kaplan-Meier survival curvers by quartile
sts graph, by(pmc) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by SM count quartile" "Shaded regions denote 95% CIs." , ///
    position(1) ring(0) size(medium)) 	 					 ///
    plot1opts(lcolor(blue) lwidth(thick)) 					 ///
    plot2opts(lcolor(red) lwidth(thick)) 					 ///
    plot3opts(lcolor(green) lwidth(thick)) 					 ///
    plot4opts(lcolor(black) lwidth(thick)) 					 ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			 ///
    legend(order( 9 "1     {&le} SM count {&le} 28" 10 "29   {&le} SM count {&le} 59" 	 ///
        11 "60   {&le} SM count {&le} 133" 12 "134 {&le} SM count {&le} 25,711 " ) cols(1) ///
            position(7) ring(0)) 						 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(quartilePMC)
graph export survival_quartile_PMC.pdf, name(quartilePMC) replace 

* Regress survival against PMC quartile
stcox i.pmc

* Regress survival against PMC quartile and cancer-type
stcox i.pmc ib3.cancer
sum PointMutationCount, detail

* Explore effect of CNA on survival

centile CNAFractionGenomeAltered, centile(25(25)75)
forvalues i = 1/3 {
    local c`i' = r(c_`i')
}
sum CNAFractionGenomeAltered
local max = r(max)
gen cna = recode( CNAFractionGenomeAltered, `c1' , `c2' , `c3' ,  `max')
tabulate cna, missing
sts test cna
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test cna if cna < `c2' + .0001
di "p = " chi2tail(`r(df)', `r(chi2)')
gen threeor4vs2 = 0 if cna > `c2'-.0001 & cna < `c2' + .0001
replace threeor4vs2 = 1 if cna > `c2' & cna != .
tabulate threeor4vs2
sts test threeor4vs2
di "p = " chi2tail(`r(df)', `r(chi2)')
sts test cna if cna > `c2' + .0001 & cna <= `max'



* Draw Kaplan-Meier survival curvers by quartile
sts graph, by(cna) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)   ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival)  ///
    subtitle("Survival curves by CNA fraction quartile", position(1) ring(0) size(medium)) ///
    plot1opts(lcolor(blue) lwidth(thick)) plot2opts(lcolor(red) lwidth(thick))	  ///
    plot3opts(lcolor(green) lwidth(thick)) plot4opts(lcolor(black) lwidth(thick)) ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none))  ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			  ///
    legend( order(9 "0        {&le} CNA fraction {&le} 0.070" 				  ///
            10 "0.070 < CNA fraction {&le} 0.194" 11 "0.194 < CNA fraction {&le} 0.385" 	  ///
            12 "0.385 < CNA fraction {&le} 1 ") cols(1) position(7) ring(0)) 		  ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(quartileCNA) 
graph export survival_quartile_CNA.pdf, name(quartileCNA) replace 

* Regress survival against CNA quartile
sort cna
egen cna_quartile = group(cna)
tabulate cna cna_quartile

stcox i.cna_quartile

* Regress survival against CNA quartile and cancer-type
stcox i.cna_quartile ib3.cancer
sum CNAFractionGenomeAltered, detail

* Explore effect of FPC on survival

centile FusionPairs, centile(25(25)75)
forvalues i = 1/3 {
    local c`i' = r(c_`i')
}
sum FusionPairs
local max = r(max)
gen fpc = recode( FusionPairs, `c1' , `c2' , `c3'  , `max')
tabulate fpc
sts test fpc
di "p = " chi2tail(`r(df)', `r(chi2)')

* Draw Kaplan-Meier survival curvers by quartile
sts graph, by(fpc) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by FG count quartile", position(1) size(medium) ring(0)) ///
    plot1opts(lcolor(blue) lwidth(thick)) plot2opts(lcolor(red) lwidth(thick))	 ///
    plot3opts(lcolor(green) lwidth(thick)) plot4opts(lcolor(black) lwidth(thick)) ///	 
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(black%30)) 			 ///
    legend( order(9 "FG count = 0" 10 "FG count = 1" 11 "FG count = 2 or 3" 			 ///
        12 "4 {&le} FG count {&le} 60") cols(1) position(7) ring(0))	 		 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(quartileFPC) 
graph export survival_quartile_FPC.pdf, name(quartileFPC) replace 

* Regress survival against FPC quartile
sort fpc
egen fpc_quartile = group(fpc)
tabulate fpc fpc_quartile

stcox i.fpc_quartile

* Run logrank tests on FPC
sts test fpc
gen fpc1 = fpc
replace fpc1 = 1 if fpc == 0
sts test fpc1
gen fpc2 = fpc1
replace fpc2 = 60 if fpc1 == 4
sts test fpc2

* Regress survival against CNA quartile and cancer-type
stcox i.fpc_quartile ib3.cancer
sum FusionPairs, detail

log close
translate overall_surv_quartiles.log overall_surv_quartiles.pdf, replace

