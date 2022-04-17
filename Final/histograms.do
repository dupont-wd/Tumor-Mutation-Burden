log using histograms.log, replace
* histograms.log

* Draw histograms of SM count, CNA fraction and FG count
* See also overall_surv_CNA.do and overall_surv_PMC.do
* This program creates panels A, B and C of Figure 1

use ../build.dta 
drop if OverallSurvivalMonths == . | overallfate == .
        
histogram CNAFractionGenomeAltered, width(.02) start(-.01) frequency 	///
    ytitle(Number of patients) ylabel(0(100)1300, angle(zero)) 		///
    xtitle(CNA fraction) 						///
    graphregion(color(white)) ysize(3) xsize(4.5) name(histogramCNA)    
graph export histogramCNA.pdf, name(histogramCNA) replace 
codebook CNAFractionGenomeAltered PointMutationCount FusionPairs

* Plot a histogram of PMC
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp')
    }
}
local start = log(0.63)
gen logPM = log(PointMutationCount)
histogram logPM, frequency xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 		///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45)) start(`start')		///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
      xtitle( SM count) ytitle(Number of patients) 			///
      ylabel(0(100)1000, angle(zero)) graphregion(color(white)) ysize(3) 	///
      xsize(4.5) name(histogramPMC)
graph export histogramPMC.pdf, name(histogramPMC) replace 

*Plot a histogram of FPC
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp'+1)
    }
}
local 0 = 0
local start = log(0.64)

gen logFPC = log( FusionPairs)    
histogram logFPC, frequency xlabel(`0' "0" `1' "1" `2' "2" `3' "3"`6' "6" 	///
          `10' "10" `20' "20" `30' "30"  `60' "60") start(`start')	///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" `40' "40" `50' "50" )	///
      xtitle( FG count) ytitle(Number of patients) 			///
      ylabel(0(200)1600, angle(zero)) ymtick(100(200)1500) 			///
      graphregion(color(white)) ysize(3) xsize(4.5) name(histogramFPC)
graph export histogramFPC.pdf, name(histogramFPC) replace 

log close
translate histograms.log histograms.pdf, replace
