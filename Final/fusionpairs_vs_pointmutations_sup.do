log using fusionpairs_vs_pointmutations_sup.log , replace
* fusionpairs_vs_pointmutations_sup.log

* This program generated Supplemental Figure S2

use ../build.dta ,clear
gen lf = log(FusionPairs +1)
twoway (scatter lf PointMutationCount)
gen lpmc = log(1+ PointMutationCount)

sum PointMutationCount
local 0 = log(1 + 0)
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(1+`tmp')
    }
}

* The preceding scatter plot obscures the number of patients in high density areas
* Let's try a density distrubution plot
sunflower lf lpmc,  ///
     ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  	///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero)) ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") 	///
      ytitle(FG count)	name(fpvspmc_sun_sup)		///
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
          order(1 "Observed" 2 3  )) 	///
               graphregion(color(white)) ysize(4.879) xsize(6.5) 
graph export fusionpairs_vs_pointmutations__sunflower_graph_sup.pdf, name(fpvspmc_sun_sup) replace            


log close
translate fusionpairs_vs_pointmutations.log fusionpairs_vs_pointmutations.pdf, replace

