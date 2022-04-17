log using fpc_vs_cna.log , replace
* fpc_vs_cna.log

* This program generates  Figure 2C

use ../build.dta ,clear
regress FusionPairs CNAFractionGenomeAltered 
twoway (scatter FusionPairs CNAFractionGenomeAltered ), name(fpccna)
pause on
gen lfpc = log(FusionPairs +1)
twoway (scatter lf CNAFractionGenomeAltered ), name(lfpccna)
regress lfpc CNAFractionGenomeAltered 
estimates stats

* Find the optimal number of cubic spline knots
mkspline _S = CNAFractionGenomeAltered, cubic nknots(3) displayknots
regress lfpc _S*
estimates stats
drop _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(4) displayknots
regress lfpc _S*
test _S2 _S3 
estimates stats
drop _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(5) displayknots
regress lfpc _S*
test _S2 _S3 _S4
estimates stats
drop _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(6) displayknots
regress lfpc _S*
estimates stats
drop _S*
mkspline _S = CNAFractionGenomeAltered, cubic nknots(7) displayknots
regress lfpc _S*
estimates stats
test _S2 _S3 _S4 _S5 _S6

* Stata does not give default knot locations for more than 7 knots
* Let's try 8 knots
centile CNAFractionGenomeAltered, centile(1.25, 98.75)
local c1 = r(c_1)
local c8 = r(c_2)
preserve
drop if CNAFractionGenomeAltered < `c1' | CNAFractionGenomeAltered > `c8'
local pct =100/7
forvalues i = 1/6 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAFractionGenomeAltered, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6')
forvalues i = 2/7 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAFractionGenomeAltered , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8')
regress lfpc _S*
estimates stats

* Let's try 9 knots
centile CNAFractionGenomeAltered, centile(1, 99)
local c1 = r(c_1)
local c9 = r(c_2)
preserve
drop if CNAFractionGenomeAltered <= `c1' | CNAFractionGenomeAltered >= `c9'
local pct =100/8
forvalues i = 1/7 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAFractionGenomeAltered, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7')
forvalues i = 2/8 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAFractionGenomeAltered , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9')
regress lfpc _S*
estimates stats

* Let's try 10 knots
centile CNAFractionGenomeAltered, centile(1, 99)
local c1 = r(c_1)
local c10 = r(c_2)
preserve
drop if CNAFractionGenomeAltered <= `c1' | CNAFractionGenomeAltered >= `c10'
local pct =100/9
forvalues i = 1/8 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAFractionGenomeAltered, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7' , `pct8')
forvalues i = 2/9 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAFractionGenomeAltered , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9' `c10')
regress lfpc _S*
estimates stats

* Let's try 11 knots
centile CNAFractionGenomeAltered, centile(1, 99)
local c1 = r(c_1)
local c11 = r(c_2)
preserve
drop if CNAFractionGenomeAltered <= `c1' | CNAFractionGenomeAltered >= `c11'
local pct =100/10
forvalues i = 1/9 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAFractionGenomeAltered, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7' , `pct8' , `pct9')
forvalues i = 2/10 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAFractionGenomeAltered , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9' `c10' `c11')
regress lfpc _S*
estimates stats


* 10 knots gives the best fit. Let's run it again
centile CNAFractionGenomeAltered, centile(1, 99)
local c1 = r(c_1)
local c10 = r(c_2)
preserve
drop if CNAFractionGenomeAltered <= `c1' | CNAFractionGenomeAltered >= `c10'
local pct =100/9
forvalues i = 1/8 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAFractionGenomeAltered, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7' , `pct8')
forvalues i = 2/9 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAFractionGenomeAltered , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9' `c10')
regress lfpc _S*
estimates stats

predict yhat, xb
replace yhat = . if yhat < 0

predict se, stdp
gen lb = max(yhat - se*1.96, 0)
gen ub = yhat + se*1.96
sum PointMutationCount
local 0 = log(1 + 0)
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(1+`tmp')
    }
}

* Plot regression curve of log FPC vs. CNA fraction
twoway (scatter lfpc CNAFractionGenomeAltered, msymbol(circle_hollow) mcolor(blue)) ///
    (rarea lb ub CNAFractionGenomeAltered, sort fcolor(yellow) lwidth(none)) 	    ///
    (line yhat CNAFractionGenomeAltered, sort lcolor(red) lwidth(medthick))	    ///
    , ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  			    ///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero)) 	    ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") 			    ///
      xtitle(CNA fraction) xlabel(0(.1)1) ytitle( Fusion pair count) 		    ///
      legend(ring(0) position(1) col(1) 					    ///
          order(1 "95% CI" 2 "Observed" 3 "Expected") symxsize(*.3)) name(fpccna1) 
graph export fpc_vs_cna_graph.pdf, name(fpccna1) replace            

* The preceding scatter plot obscures the number of patients in high density areas
* Let's try a density distrubution plot. We draw this plot twice with two different
* legend components that will be spliced together with photoshop
sunflower lfpc CNAFractionGenomeAltered, 					     ///
    addplot((rarea lb ub CNAFractionGenomeAltered, sort fcolor(yellow) lwidth(none)) ///
    (line yhat CNAFractionGenomeAltered, sort lcolor(red) lwidth(medthick)))	     ///
     ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  			     ///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero)) 	     ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") xtitle(CNA fraction) 	     ///
      xlabel(0(.1)1) ytitle( FG count) legend(col(1) ring(0) position(1)    ///
          order(1 "Observed" 2 3 ) symysize(*.9) symxsize(*.9) size(*.9) )           ///
      graphregion(color(white)) ysize(3) xsize(4.5) name(fpccna_sun1) 
graph export fpc_vs_cna_sunflower_tmp1.pdf, name(fpccna_sun1) replace            

sunflower lfpc CNAFractionGenomeAltered, 					     ///
    addplot((rarea lb ub CNAFractionGenomeAltered, sort fcolor(yellow) lwidth(none)) ///
    (line yhat CNAFractionGenomeAltered, sort lcolor(red) lwidth(medthick)))	     ///
     ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  			     ///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero))           ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") xtitle(CNA fraction)	     ///
      xlabel(0(.1)1) ytitle( FG count) legend(col(1) ring(0) position(10)   ///
          order( 4 "95% CI" 5 "Expected" ) symysize(*.9)  size(*.9) ) 		     ///
      graphregion(color(white)) ysize(3) xsize(4.5) name(fpccna_sun2)
graph export fpc_vs_cna_sunflower_tmp2.pdf, name(fpccna_sun2) replace            

log close
translate fpc_vs_cna.log fpc_vs_cna.pdf, replace
*/
