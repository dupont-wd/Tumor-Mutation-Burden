log using fusionpairs_vs_pointmutations.log , replace
* fusionpairs_vs_pointmutations.log

* This program generated Supplemental Figure 1B

use ../build.dta ,clear
regress PointMutationCount FusionPairs
twoway (scatter PointMutationCount FusionPairs)
pause on
regress  FusionPairs PointMutationCount
twoway (scatter FusionPairs PointMutationCount)
gen lf = log(FusionPairs +1)
twoway (scatter lf PointMutationCount)
gen lpmc = log(1+ PointMutationCount)
twoway (scatter lf lpmc)

* Find the optimal number of cubic spline knots
mkspline _S = lpmc, cubic nknots(3)
regress lf _S*
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(4)
regress lf _S*
test _S2 _S3 
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(5)
regress lf _S*
test _S2 _S3 _S4
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(6)
regress lf _S*
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(7)
regress lf _S*
estimates stats
test _S2 _S3 _S4 _S5 _S6

* Stata does not give default knot locations for more than 7 knots
* Let's try 8 knots
centile lpmc, centile(1.25, 98.75)
local c1 = r(c_1)
local c8 = r(c_2)
preserve
drop if lpmc < `c1' | lpmc > `c8'
local pct =100/7
forvalues i = 1/6 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile lpmc, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6')
forvalues i = 2/7 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = lpmc , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8')
regress lf _S*
estimates stats

* Let's try 9 knots
centile lpmc, centile(1, 99)
local c1 = r(c_1)
local c9 = r(c_2)
preserve
drop if lpmc < `c1' | lpmc > `c9'
local pct =100/8
forvalues i = 1/7 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile lpmc, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7')
forvalues i = 2/8 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = lpmc , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9')
regress lf _S*
estimates stats

* 8 knots gives the best fit. Let's run it again
centile lpmc, centile(1.25, 98.75)
local c1 = r(c_1)
local c8 = r(c_2)
preserve
drop if lpmc < `c1' | lpmc > `c8'
local pct =100/7
forvalues i = 1/6 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile lpmc, centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6')
forvalues i = 2/7 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = lpmc , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8')
regress lf _S*

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

* Plot regression curve of log fusion pair count against log point mutation count
twoway (scatter lf lpmc, msymbol(circle_hollow) mcolor(blue))	///
    (rarea lb ub lpmc, sort fcolor(yellow) lwidth(none)) 	///
    (line yhat lpmc, sort lcolor(red) lwidth(medthick))		///
    , ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  	///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero)) ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") 	///
      ytitle(Fusion pair count)	name(fpvspmc)			///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6"  			///
          `10' "10" `20' "20" `30' "30"  `60' "60"		///
          `100' "100"  `200' "200" `300' "300"  `600' "600"	///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" ///
          `10000' "10000" `20000' "20000" , angle(45))	       	///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 		///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 	///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" ///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      xtitle( Point mutation count) legend(ring(0) position(1) col(1) ///
          order(1 "95% CI" 2 "Observed" 3 "Expected") symxsize(*.75)) 
graph export fusionpairs_vs_pointmutations_graph.pdf, name(fpvspmc) replace            

* The preceding scatter plot obscures the number of patients in high density areas
* Let's try a density distrubution plot
sunflower lf lpmc, addplot((rarea lb ub lpmc, sort fcolor(yellow) lwidth(none)) ///
    (line yhat lpmc, sort lcolor(red) lwidth(medthick)))	///
     ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  	///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero)) ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") 	///
      ytitle(FG count)	name(fpvspmc_sun)		///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6"  			///
          `10' "10" `20' "20" `30' "30"  `60' "60"		///
          `100' "100"  `200' "200" `300' "300"  `600' "600"	///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" ///
          `10000' "10000" `20000' "20000" , angle(45))	       	///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 		///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 	///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" ///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      xtitle(SM count) legend(ring(0) position(1) col(1) ///
          order(1 "Observed" 2 3 4 "95% CI" 5 "Expected" ) 	///
              symxsize(*.25)) graphregion(color(white)) ysize(3) xsize(4.5) 
graph export fusionpairs_vs_pointmutations__sunflower_graph.pdf, name(fpvspmc_sun) replace            


log close
translate fusionpairs_vs_pointmutations.log fusionpairs_vs_pointmutations.pdf, replace

