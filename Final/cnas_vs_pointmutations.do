log using cnas_vs_pointmutations.log , replace
* cnas_vs_pointmutations.log

*This program generates supplemental figure 1A

use ../build.dta ,clear
regress PointMutationCount CNAFractionGenomeAltered
twoway (scatter PointMutationCount FusionPairs)
pause on
regress  CNAFractionGenomeAltered PointMutationCount
twoway (scatter FusionPairs PointMutationCount)
gen lpmc = log(1+ PointMutationCount)
twoway (scatter CNAFractionGenomeAltered lpmc)

regress CNAFractionGenomeAltered lpmc
estimates stats

* Find the optimal number of cubic spline knots
mkspline _S = lpmc, cubic nknots(3)
regress CNAFractionGenomeAltered _S*
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(4)
regress CNAFractionGenomeAltered _S*
test _S2 _S3 
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(5)
regress CNAFractionGenomeAltered _S*
test _S2 _S3 _S4
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(6)
regress CNAFractionGenomeAltered _S*
test _S2 _S3 _S4 _S5
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(7)
regress CNAFractionGenomeAltered _S*
estimates stats
test _S2 _S3 _S4 _S5 _S6

* 6 knots give best fit.
drop _S*
mkspline _S = lpmc, cubic nknots(6)
regress CNAFractionGenomeAltered _S*
test _S2 _S3 _S4 _S5

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

* Plot regression curve of CNA fraction against 
* log point mutation count
twoway (scatter CNAFractionGenomeAltered lpmc, 			///
        msymbol(circle_hollow) mcolor(blue)) 			///
       (rarea lb ub lpmc, sort fcolor(yellow) lwidth(none)) 	///    
       (line yhat lpmc, sort lcolor(red) lwidth(medthick))	///
       , ylabel(0 (.1)1 , angle(zero)) 				///
      ytitle(CNA fraction)	name(cnapmc1)			///
     xlabel(`1' "1" `2' "2" `3' "3"`6' "6"  			///
          `10' "10" `20' "20" `30' "30"  `60' "60"		///
          `100' "100"  `200' "200" `300' "300"  `600' "600"	///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' 	///
          "6000" `10000' "10000" `20000' "20000" , angle(45))	///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 		///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 	///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' ///
          "900" `4000' "4000" `5000' "5000" `7000' "7000" 	///
          `8000' "8000" `9000' "9000" ) 			///
      xtitle( Point mutation count) legend(ring(0) position(1)  ///
          col(1) order(1 "95% CI" 2 "Observed" 3 "Expected") 	///
          symxsize(*.75)) 
          
graph export cna_vs_pointmutations_graph.pdf, name(cnapmc1) replace            

* The preceding scatter plot obscures the number of patients in 
* high density areas. Let's try a density distrubution plot
sort lpmc
sunflower CNAFractionGenomeAltered lpmc, addplot((rarea lb ub 	///
        lpmc, fcolor(yellow) lwidth(none)) 			///
    (line yhat lpmc,  lcolor(red) lwidth(medthick)))		///
    ylabel(0(.1)1 , angle(zero)) ytitle(CNA fraction)		///
     name(cnapmc_sun) xlabel(`1' "1" `2' "2" `3' "3"`6' "6"  	///
          `10' "10" `20' "20" `30' "30"  `60' "60"		///
          `100' "100"  `200' "200" `300' "300"  `600' "600"	///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' 	///
          "6000" `10000' "10000" `20000' "20000" , angle(45))   ///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 		///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 	///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' ///
          "900" `4000' "4000" `5000' "5000" `7000' "7000" 	///
          `8000' "8000" `9000' "9000") 				///
      xtitle(SM count) legend(ring(0) position(1) 	///
          col(1) order(1 "Observed" 2 3 4 "95% CI" 		///
          5 "Expected") symxsize(*.5)) graphregion(color(white)) ysize(3) xsize(4.5) 
graph export cna_vs_pointmutations__sunflower_graph.pdf, name(cnapmc_sun) replace            


log close
translate cnas_vs_pointmutations.log cnas_vs_pointmutations.pdf, replace
