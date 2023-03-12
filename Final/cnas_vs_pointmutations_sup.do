log using cnas_vs_pointmutations_sup.log , replace
* cnas_vs_pointmutations_sup.log

*This program generates supplemental figure S1

use ../build.dta ,clear
gen lpmc = log(1+ PointMutationCount)

sort lpmc

sum PointMutationCount
local 0 = log(1 + 0)
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(1+`tmp')
    }
}

sunflower CNAFractionGenomeAltered lpmc ,			///
    ylabel(0(.1)1 , angle(zero)) ytitle(CNA fraction)		///
     name(cnapmc_sun_sup) xlabel(`1' "1" `2' "2" `3' "3"`6' "6" ///
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
          col(1) order(1 "Observed" 2 3)) graphregion(color(white)) ysize(4.875) xsize(6.5) 
graph export cna_vs_pointmutations__sunflower_graph_sup.pdf, name(cnapmc_sun_sup) replace            


log close
translate cnas_vs_pointmutations_sup.log cnas_vs_pointmutations_sup.pdf, replace
