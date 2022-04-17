log using tcga_tmb.log, replace
* tcga_tmb.log

* Analyze the effect of TCGA TMB on survival

use "C:\Users\dupontwd\Dropbox\OFFICE_dropbox\Parl\Mutational burden & cancer survival\build.dta", clear
sort PatientID
by PatientID: gen duplicates= _N
codebook duplicates
merge 1:1 PatientID SampleID using ..\tcga_tmb
* The following patients are either missing follow-up time or they are from patients
* with multiple samples where the listed sample is not the primary sample
list PatientID SampleID  OverallSurvivalStatus OverallSurvivalMonths_string  if _merge == 2

* keep patients in build.dta

keep if _merge == 3

scatter PointMutationCount tcga_tmb, msymbol(Oh) name(SMvsTCGA)
correlate PointMutationCount tcga_tmb
list PointMutationCount tcga_tmb if tcga_tmb > 800 & tcga_tmb != .
local adj = 25711/ 861.86667
di "SM count / TMBcol ratio = `adj'"
local low_slope =  15825  /    1052.2
scatter PointMutationCount tcga_tmb if PointMutationCount/ tcga_tmb > (`adj' + `low_slope')/2, msymbol(Oh) name(bigslope)
correlate PointMutationCount tcga_tmb if PointMutationCount/ tcga_tmb > (`adj' + `low_slope')/2
scatter PointMutationCount tcga_tmb if PointMutationCount/ tcga_tmb < (`adj' + `low_slope')/2, msymbol(Oh) name(smallslope)
correlate PointMutationCount tcga_tmb if PointMutationCount/ tcga_tmb < (`adj' + `low_slope')/2
gen bigslope =  PointMutationCount/ tcga_tmb >  (`adj' + `low_slope')/2

regress PointMutationCount tcga_tmb if bigslope
local adjslope = _b[tcga_tmb]
di " adj slope = `adjslope'"
tabulate bigslope
tabulate cancer bigslope, col
tabulate CancerStage bigslope, col
tabulate TumorGrade bigslope, col
ttest PointMutationCount, by(bigslope) unequal
ttest CNAFractionGenomeAltered , by(bigslope) unequal
ttest FusionPairs , by(bigslope) unequal
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
predict loghr, xb
sum loghr
gen hr = exp(loghr)
line hr logPM , sort lcolor(red) lwidth(medium) name(logPMhr)
sum tcga_tmb
gen log_tcga_tmb = log( tcga_tmb)
codebook tcga_tmb
drop _S*
mkspline _S = log_tcga_tmb , cubic nknots(3) displayknots
stcox _S*
predict loghr_tmb, xb
gen hr_tmb = exp(loghr_tmb)

foreach i of numlist  2 20  200 2000 {
    forvalues j = 1/9 {
* local names may not contain periods    
        local tmp = `i'*`j'
        di "tmp = `tmp'"
        local `tmp'  = log(`tmp'/100)
        di "log tmp/100 = ``tmp''"      
    }
}

line hr_tmb log_tcga_tmb ,sort lcolor(blue) lwidth(medium) name(tcga_hr1)	///
     xlabel(`4' "0.04" `20' "0.2" `40' "0.4" , angle(45) )/*`3' "3"`6' "6" 					///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45) labcolor(red))	       			///
      xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 				///
          `40' "40" `50' "50" `70' "70" `80' "80" `90' "90" 			///
          `400' "400" `500' "500" `700' "700" `800' "800" `900' "900" 		///
          `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" `9000' "9000") ///
*/
sum hr_tmb
* When tcga_tmb = 0.0333, hr_tmb = .0763806
gen hr_tmb_adj = hr_tmb/ .0763806
twoway (line hr_tmb_adj log_tcga_tmb ,sort lcolor(blue) lwidth(medium)) 	///
    (line hr logPM , sort lcolor(red) lwidth(medium)) 			///
  , name(tcga_hr2)
gen tcga_tmb_adj = tcga_tmb *`adjslope'
sum tcga_tmb_adj PointMutationCount if bigslope
gen log_tcga_tmb_adj = log(tcga_tmb_adj)

foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp')
    }
}

* Draw graph with xaxis labels for SM count
twoway (line hr_tmb_adj log_tcga_tmb_adj ,sort lcolor(blue) lwidth(medium) )	///
    (line hr logPM , sort lcolor(red) lwidth(medium))			///
  , name(tcga_hr3) ytitle(Mortal hazard ratio) 				///
      ylabel(0(2)20, angle(zero)) xtitle("SM count", color(red))      		///
      xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 					///
          `10' "10" `20' "20" `30' "30"  `60' "60"				///
          `100' "100"  `200' "200" `300' "300"  `600' "600"			///
          `1000' "1000" `2000' "2000" `3000' "3000"  `6000' "6000" 		///
          `10000' "10000" `20000' "20000" , angle(45) labcolor(red))	       			///
      xmtick(`4' `5'  `7' `8'  `9' `40'  `50'  `70'  `80' `90'  			///
          `400'  `500'  `700'  `800'  `900' `4000' `5000' `7000' `8000' `9000') ///
    legend(ring(0) position(10) cols(1) order(2 "SM count" 1  	///
    "TMB{sub:RATE}" ) symxsize(*.4)) ///
    graphregion(color(white)) ysize(3) xsize(3.75)
graph export TMBcol1.pdf, name(tcga_hr3) replace  

foreach i of numlist  2 20  200 2000 20000{
    forvalues j = 1/9 {
* local names may not contain periods    
        local tmp = `i'*`j'
        di "tmp = `tmp'"
        local `tmp'  = log(`adjslope'*`tmp'/100)
        di "log adjslope*tmp/100 = ``tmp''"      
    }
}
* Draw graph with xaxis labels for TCGA TMB
twoway (line hr_tmb_adj log_tcga_tmb_adj ,sort lcolor(blue) lwidth(medium) )	///
    (line hr logPM , sort lcolor(red) lwidth(medium))			///
  , name(tcga_hr4) ytitle(Mortal hazard ratio) 				///
      ylabel(0(2)20, angle(zero)) xtitle("TMB{sub:RATE}", color(blue)) ///
      xlabel(`4' "0.04" `8'  "0.08"    `20' "0.2"    `40' "0.4" ///
                       `80'   "0.8"   `200' "2.0"   `400' "4.0"		///
                      `800'   "8.0"  `2000'  "20"  `4000'  "40"		///
                     `8000'    "80" `20000' "200" `40000' "400" `80000' "800" ///
                     , angle(45) labcolor(blue) ) ///
      xmtick(`12'  `16'  `120' `160'  `1200'   `1600'    `12000'   `16000' )   ///  
    legend(off)									///
    graphregion(color(white)) ysize(3) xsize(3.75)
graph export TMBcol2.pdf, name(tcga_hr4) replace  
   
    */
    
log close
