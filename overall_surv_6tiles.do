log using overall_surv_6tiles.log, replace
* overall_surv_6tiles.log

* Draw survival plots of SM count and CNA fraction by sextiles
* This program generates Figures 2A and 2B

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
centile PointMutationCount, centile(16.6667(16.6667)83.3335 )
forvalues i = 1/5 {
    local c`i' = r(c_`i')
}
sum PointMutationCount
local max = r(max)
local c6 = `max'
local c0 = 0

* Calculate median PointMutationCount for patients in  sextiles 1 -- 6
forvalues i = 1/6 {
    di "Summarize PointMutationCount in sextile `i'"
    local im1 = `i'-1
    di "Sextile `i' goes from `c`im1'' to `c`i''"
    sum PointMutationCount if (PointMutationCount <= `c`i'' & PointMutationCount > `c`im1''), detail
    di "Median for sextile `i' = " r(p50)
}


gen pmc = recode( PointMutationCount, `c1' , `c2' , `c3'  ,`c4' , `c5', `max')
di "local values `c1' `c2' `c3' `c4' `c5' `max'"
codebook PointMutationCount
tabulate pmc , missing
sts test pmc
di "p = " chi2tail(`r(df)', `r(chi2)')

sts test pmc if pmc > `c1' & pmc != .
di "p = " chi2tail(`r(df)', `r(chi2)')

forvalues i=1/5 {
    di " "
    di "Sextile `i' vs sextile 6"
    sts test pmc if pmc == `c`i'' | pmc == `max'
    di "p = " chi2tail(`r(df)', `r(chi2)')
}
sts test pmc if pmc == `c3' | pmc == `c4' | pmc == `c5'
di "p = " chi2tail(`r(df)', `r(chi2)')

sts test pmc if pmc == `c1' | pmc == `c2' 
di "p = " chi2tail(`r(df)', `r(chi2)')

sts test pmc if pmc == `c2' | pmc == `c3' 
di "p = " chi2tail(`r(df)', `r(chi2)')

sts test pmc if pmc == `c1' | pmc == `c2' 
di "p = " chi2tail(`r(df)', `r(chi2)')


gen first_vs_rest = recode(pmc, `c1', `max')
tabulate first_vs_rest
sts test first_vs_rest
di "p = " chi2tail(`r(df)', `r(chi2)')
preserve
drop if pmc < `c3' | pmc == .

gen last_vs_s3to5 = pmc == `max'
tabulate last_vs_s3to5
sts test last_vs_s3to5
di "p = " chi2tail(`r(df)', `r(chi2)')
restore

*gen last_vs_3rd_or_4th = 2 if pmc== `max'
*replace last_vs_3rd_or_4th = 1 if pmc==`c3' | pmc == `c4'
*tabulate last_vs_3rd_or_4th
*sts test last_vs_3rd_or_4th

* Draw Kaplan-Meier survival curvers by 6tile of PM count
sts graph, by(pmc) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by SM count sextile" "Shaded regions denote 95% CIs." , ///
    position(1) ring(0) size(medium)) 	 					 ///
    plot1opts(lcolor(blue) lwidth(thick)) 					 ///
    plot2opts(lcolor(red) lwidth(thick)) 					 ///
    plot3opts(lcolor(green) lwidth(thick)) 					 ///
    plot4opts(lcolor(cyan) lwidth(thick)) 					 ///	
    plot5opts(lcolor(magenta) lwidth(thick)) 					 ///	
    plot6opts(lcolor(black) lwidth(thick)) 					 ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(cyan%30)) 				 ///
    ci5opts(fcolor(magenta%30)) ci6opts(fcolor(black%30)) ///
    legend(order(13 "    1 {&le} SM count {&le} 20"  14 "  21 {&le} SM count {&le} 36" ///
                  15 "  37 {&le} SM count {&le} 59"  16 "  60 {&le} SM count {&le} 97" ///
                  17 "  98 {&le} SM count {&le} 211" 18  "212 {&le} SM count {&le} 25,711 " ///
       ) cols(1)  position(7) ring(0) symxsize(*0.5))	 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(sextilePMC)     
    
graph export survival_6tile_PMC.pdf, name(sextilePMC) replace 

sts graph, by(pmc) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by SM count sextile" "Shaded regions denote 95% CIs." , ///
    position(1) ring(0) size(medium)) 	 					 ///
    plot1opts(lcolor(blue) lwidth(thick)) 					 ///
    plot2opts(lcolor(red) lwidth(thick)) 					 ///
    plot3opts(lcolor(green) lwidth(thick)) 					 ///
    plot4opts(lcolor(cyan) lwidth(thick)) 					 ///	
    plot5opts(lcolor(magenta) lwidth(thick)) 					 ///	
    plot6opts(lcolor(black) lwidth(thick)) 					 ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(cyan%30)) 				 ///
    ci5opts(fcolor(magenta%30)) ci6opts(fcolor(black%30)) ///
    legend(off)	 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(sextilePMCnoleg)     
    
graph export survival_6tile_PMCnoleg.pdf, name(sextilePMCnoleg) replace 

* Regress survival against PMC quartile
stcox i.pmc

stcox ib`max'.pmc 
sum PointMutationCount, detail

* Explore the effect of CNA fraction on survival @@@@@@@@@@@@@@@@@@@@@@@

* We need to convert CNA fraction to a variable with integer values
gen cnaf = round(CNAFractionGenomeAltered*1000)
centile cnaf, centile(16.6667(16.6667)83.3335 )
forvalues i = 1/5 {
    local c`i' = round(r(c_`i'))
}
sum CNAFractionGenomeAltered cnaf
local max = r(max)
/*
* Calculate median CNAFractionGenomeAltered for patients in the first quartile
sum CNAFractionGenomeAltered if CNAFractionGenomeAltered <= `c1', detail
di "Median for first quartile = " r(p50)

* Calculate median CNAFractionGenomeAltered for patients in the second quartile
sum CNAFractionGenomeAltered if CNAFractionGenomeAltered > `c1' & CNAFractionGenomeAltered <= `c2', detail
di "Median for 2nd quartile = " r(p50)

* Calculate median CNAFractionGenomeAltered for patients in the third quartile
sum CNAFractionGenomeAltered if CNAFractionGenomeAltered > `c2' & CNAFractionGenomeAltered <= `c3', detail
di "Median for 3rd quartile = " r(p50)

* Calculate median CNAFractionGenomeAltered for patients in the fourth quartile
sum CNAFractionGenomeAltered if CNAFractionGenomeAltered > `c3' & CNAFractionGenomeAltered !=., detail
di "Median for 4th quartile = " r(p50)
*/


gen cna = recode( cnaf, `c1' , `c2' , `c3'  ,`c4' , `c5', `max')
di "local values `c1' `c2' `c3' `c4' `c5' `max'"
codebook CNAFractionGenomeAltered cnaf
tabulate cna , missing
sts test cna
di "p = " chi2tail(`r(df)', `r(chi2)')

sts test cna if cna > `c1' & cna != .
di "p = " chi2tail(`r(df)', `r(chi2)')

drop first_vs_rest
gen first_vs_rest = recode(cna, `c1', `max')
tabulate first_vs_rest
sts test first_vs_rest
di "p = " chi2tail(`r(df)', `r(chi2)')

forvalues i=1/5 {
    di " "
    di "Sextile `i' vs sextile 6"
    sts test cna if cna == `c`i'' | cna == `max'
    di "p = " chi2tail(`r(df)', `r(chi2)')
}


*gen last_vs_3rd_or_4th = 2 if cna== `max'
*replace last_vs_3rd_or_4th = 1 if cna==`c3' | cna == `c4'
*tabulate last_vs_3rd_or_4th
*sts test last_vs_3rd_or_4th

* Draw Kaplan-Meier survival curvers by 6tile of CNA fraction
sts graph, by(cna) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by CNA fraction sextile"  , ///
    position(1) ring(0) size(medium)) 	 					 ///
    plot1opts(lcolor(blue) lwidth(thick)) 					 ///
    plot2opts(lcolor(red) lwidth(thick)) 					 ///
    plot3opts(lcolor(green) lwidth(thick)) 					 ///
    plot4opts(lcolor(cyan) lwidth(thick)) 					 ///	
    plot5opts(lcolor(magenta) lwidth(thick)) 					 ///	
    plot6opts(lcolor(black) lwidth(thick)) 					 ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(cyan%30)) ci5opts(fcolor(magenta%30)) ci6opts(fcolor(black%30)) ///
    legend(order(13 "0.000 {&le} CNA fraction {&le} 0.033"  14 "0.033 < CNA fraction {&le} 0.112" ///
                 15 "0.112 < CNA fraction {&le} 0.194"  16 "0.194 < CNA fraction {&le} 0.308" ///
                 17 "0.308 < CNA fraction {&le} 0.480" 18  "0.480 < CNA fraction {&le} 1.000 " ///
       ) cols(1)  position(7) ring(0) symxsize(*0.5))	 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(sextileCNA)     ///
    
graph export survival_6tile_CNA.pdf, name(sextileCNA) replace 

sts graph, by(cna) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by CNA fraction sextile" , ///
    position(1) ring(0) size(medium)) 	 					 ///
    plot1opts(lcolor(blue) lwidth(thick)) 					 ///
    plot2opts(lcolor(red) lwidth(thick)) 					 ///
    plot3opts(lcolor(green) lwidth(thick)) 					 ///
    plot4opts(lcolor(cyan) lwidth(thick)) 					 ///	
    plot5opts(lcolor(magenta) lwidth(thick)) 					 ///	
    plot6opts(lcolor(black) lwidth(thick)) 					 ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(cyan%30)) ci5opts(fcolor(magenta%30)) ci6opts(fcolor(black%30)) ///
    legend(off)	 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(sextileCNAnoleg)     ///
    
graph export survival_6tile_CNAnoleg.pdf, name(sextileCNAnoleg) replace 


* Regress survival against PMC quartile
stcox i.cna

stcox ib`max'.cna 
sum CNAFractionGenomeAltered, detail


* Explore effect of FPC on survival by sextiles
* Note that the first two sextiles have an FG  count of zero

centile FusionPairs, centile(16.6667(16.6667)83.3335 )
forvalues i = 1/5 {
    local c`i' = r(c_`i')
}
sum FusionPairs
local max = r(max)
gen fpc = recode( FusionPairs, `c1', `c2', `c3', `c4', `c5', `max')
tabulate fpc
sts test fpc
di "p = " chi2tail(`r(df)', `r(chi2)')

* Draw Kaplan-Meier survival curvers by sextile
sts graph, by(fpc) ci ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis)  ///
    xlabel(0(1)15) title("", size(zero)) ytitle(Probability of overall survival) ///
    subtitle("Survival curves by FG count sextile", position(1) size(medium) ring(0)) ///
    plot1opts(lcolor(blue) lwidth(thick)) 					 ///
    plot2opts(lcolor(red) lwidth(thick)) 					 ///
    plot3opts(lcolor(green) lwidth(thick)) 					 ///
    plot4opts(lcolor(cyan) lwidth(thick)) 					 ///	
    plot5opts(lcolor(black) lwidth(thick)) 					 ///	
    ci1opts(fcolor(blue%30)) ci2opts(fcolor(red%30) lcolor(red%30) lwidth(none)) ///
    ci3opts(fcolor(green%30)) ci4opts(fcolor(cyan%30)) ci5opts(fcolor(black%30))  ///
    legend( order(11 "FG count = 0" 12 "FG count = 1" 13 "FG count = 2 " 14 "FG count = 3 or 4" ///
        15 "5 {&le} FG count {&le} 60") cols(1) position(7) ring(0) region(color(white)))	 ///
    graphregion(color(white)) ysize(3) xsize(3.75) name(quartileFPC) 
graph export survival_quartile_FPC.pdf, name(quartileFPC) replace 

* Regress survival against FPC quartile
sort fpc
egen fpc_quartile = group(fpc)
tabulate fpc fpc_quartile

stcox i.fpc_quartile

* Run logrank tests on FPC
sts test fpc
di "p = " chi2tail(`r(df)', `r(chi2)')

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
translate overall_surv_6tiles.log overall_surv_6tiles.pdf, replace

