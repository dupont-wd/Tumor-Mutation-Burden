log using fusionpairs_vs_CNAfraction_sup.log , replace
* fusionpairs_vs_CNAfraction_sup.log

* This program generates Supplemental Figure S3

use ../build.dta ,clear
gen lfpc = log(FusionPairs +1)
sum PointMutationCount
local 0 = log(1 + 0)
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(1+`tmp')
    }
}

sunflower lfpc CNAFractionGenomeAltered, 					     ///
     ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  			     ///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero))           ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") xtitle(CNA fraction)	     ///
      xlabel(0(.1)1) ytitle( FG count) legend(col(1) ring(0) position(1) 		///
      order(1 "Observed" 2 3)  symysize(*.9) symxsize(*.9) size(*.9))  			///
      graphregion(color(white)) ysize(4.333) xsize(6.5) name(f1)		
      
graph export fpc_vs_cna_sunflower_graph_sup.pdf, name(f1) replace            
/*sunflower lfpc CNAFractionGenomeAltered, 					     ///
     ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  			     ///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero))           ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") xtitle(CNA fraction)	     ///
      xlabel(0(.1)1) ytitle( FG count) legend(off)  					///
      graphregion(color(white)) ysize(4.333) xsize(6.5) name(f2)
graph export fpc_vs_cna_sunflower_graph_sup2.pdf, name(fpccna_sup2) replace            
*/
log close
translate fpc_vs_cna.log fpc_vs_cna.pdf, replace
*/
