log using fusionpairs_vs__CNAfraction_heatmap.log , replace
* fusionpairs_vs_pointmutations_heatmap.log
pause on
*This program generates  figure 1F

use ../build.dta ,clear
gen CNAfraction = CNAFractionGenomeAltered
gen FPcount = FusionPairs
keep CNAfraction FPcount
drop if CNAfraction ==.
drop if FPcount ==.
gen logFPcount = log(FPcount +1)
codebook  CNAfraction  FPcount logFPcount

* Select evenly spaced bins for CNAfraction
sum CNAfraction, detail
local maxx = r(max)
local binsx = 80
local binwidthCNA = round(`maxx'/`binsx',0.001)
gen CNAfraction_binned = `binwidthCNA'/2
local binsxm1= `binsx'-1
forvalues i = 1/`binsxm1' {
    replace CNAfraction_binned = `binwidthCNA'/2 + `i'*`binwidthCNA' if CNAfraction >`i'*`binwidthCNA' 
}
tabulate  CNAfraction_binned
histogram CNAfraction_binned, freq width(`binwidth') name(CNA)

* Select evenly spaced bins for logFPcount
sum logFPcount, detail
local maxy = r(max)
local binsy = 60
local binwidthFP = round(`maxy'/`binsy',0.001)
gen logFPcount_binned = `binwidthFP'/2
local binsym1= `binsy'-1
forvalues i = 1/`binsym1' {
    replace logFPcount_binned = `binwidthFP'/2 + `i'*`binwidthFP' if logFPcount >`i'*`binwidthFP' 
}

tabulate logFPcount_binned
histogram logFPcount_binned, freq width(`binwidth') name(FPc)
sort CNAfraction_binned logFPcount_binned
save tmp, replace
local obs = `binsx' * `binsy'

clear
di "`binsx' `binsy' `obs'"
set obs `obs'
gen CNAfraction_binned = .
gen logFPcount_binned = .
forvalues i =1/`binsx' {
    quietly replace CNAfraction_binned = -`binwidthCNA'/2 + `i'*`binwidthCNA' if `i' == 1 + int((_n-1)/`binsy')
    forvalues j = 1/`binsy' {
        quietly replace logFPcount_binned = -`binwidthFP'/2 + `binwidthFP'*`j' if (`j' == mod(_n,`binsy') | mod(_n,`binsy') ==0)
    }
}
merge 1:m CNAfraction_binned logFPcount_binned using tmp
gen observation = _merge==3
collapse (sum) count = observation , by( CNAfraction_binned logFPcount_binned)

gen logcount = log(count)
sum logcount
replace logcount = -0.5 if count == 0
twoway (contour logcount logFPcount_binned CNAfraction_binned, heatmap ccuts(-.5 (.5) 7.5)  ccolors(blue "0 59 255" 	///
    "0 120 255" cyan*1.4 cyan cyan*.6 green*1.8 green*1.4 green  )   scolor(blue) ecolor(red)	name(heat1))        ///

sum logcount
local inc = r(max)/15
di `inc'
replace logcount = -`inc' if count == 0


local ivalues "`inc' "
forvalues i =1/15 {
    local nexti = `i'*`inc' + `inc'
    local ivalues "`ivalues' `nexti'"
}
di "`ivalues'"

local j = 0
local allcuts "-`inc' 0 "
di "`allcuts'"
foreach i in  `ivalues' {
    di "original cut = " `i'
    local cut = exp(`i')
    local j = `j'+1
    di "original cut exponentiated " `cut'
    local cut`j' `i'
    local allcuts "`allcuts' `cut`j'' "
}
forvalues j=1/15 {
    di "`cut`j''"
}  
di "`cut15'"
di "`cut15'"
* increase cut15 slightly to account to the fact that there is a cell with 950 patients
local cut15 = log(951)
di "`allcuts'"
di "`cut1'"

* Define x axis labels
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp')
        local x`tmp' = log(`tmp')
*       di "tmp `tmp' xtmp `x`tmp''"
    }
}

* Define y axis labels
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp' + 1)
        local y`tmp' = log(`tmp' + 1)
*       di "tmp `tmp' xtmp `x`tmp''"
    }
}

*di "x10 `x10' x600 `x600'"
local y0 = 0
twoway (contour logcount logFPcount_binned CNAfraction_binned 	///
  , heatmap ccuts(`allcuts')  ccolors(white white blue blue*.5 	///
    yellow yellow cyan cyan*1.5 green green*1.5    	///
     magenta*.5 magenta magenta*1.5  yellow*1.5 orange red red*1.5)	///   
    zlabel( -`inc' "    0" 0 "    1" ///
        `cut1'  "    2"   ///
        `cut2'  "    3" ///  
        `cut3'  "    4" ///  
        `cut4'  "    7" ///  
        `cut5'  "  10" ///  
        `cut6'  "  16" ///  
        `cut7'  "  25" ///  
        `cut8'  "  39" ///  
        `cut9'  "  62" ///  
        `cut10' "  97" ///  
        `cut11' "153" ///  
        `cut12' "242" ///  
        `cut13' "381" ///  
        `cut14' "602"  ///
        `cut15' "950" )  ///
     ylabel(`y0' "0"  `y1' "1" `y2' "2" `y3' "3"`y4' "4" `y6' "6"  ///
         `y10' "10" `y20' "20" `y30' "30" `y40' "40" `y60' "60" , angle(zero)) ///
     ymtick(`y5'  `y7'  `y8' `y9' `y50')			///
     xlabel(0 (.1) 1) xmtick(.05 (.1) .9) 			///
    ztitle(Number of patients) xtitle(CNA fraction) ytitle(FG count) ///
    name(heat2) graphregion(color(white)) ysize(2.8) xsize(3.75) ) 

use tmp, clear 
twoway (scatter FPcount CNAfraction), name(fpccna)

gen lfpc = log(FPcount +1)
twoway (scatter lfpc CNAfraction  ), name(lfpccna)
regress lfpc CNAfraction  
estimates stats

* Find the optimal number of cubic spline knots
mkspline _S = CNAfraction , cubic nknots(3) displayknots
regress lfpc _S*
estimates stats
drop _S*
mkspline _S = CNAfraction , cubic nknots(4) displayknots
regress lfpc _S*
test _S2 _S3 
estimates stats
drop _S*
mkspline _S = CNAfraction , cubic nknots(5) displayknots
regress lfpc _S*
test _S2 _S3 _S4
estimates stats
drop _S*
mkspline _S = CNAfraction , cubic nknots(6) displayknots
regress lfpc _S*
estimates stats
drop _S*
mkspline _S = CNAfraction , cubic nknots(7) displayknots
regress lfpc _S*
estimates stats
test _S2 _S3 _S4 _S5 _S6

* Stata does not give default knot locations for more than 7 knots
* Let's try 8 knots
centile CNAfraction , centile(1.25, 98.75)
local c1 = r(c_1)
local c8 = r(c_2)
preserve
drop if CNAfraction  < `c1' | CNAfraction  > `c8'
local pct =100/7
forvalues i = 1/6 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAfraction , centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6')
forvalues i = 2/7 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAfraction  , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8')
regress lfpc _S*
estimates stats

* Let's try 9 knots
centile CNAfraction , centile(1, 99)
local c1 = r(c_1)
local c9 = r(c_2)
preserve
drop if CNAfraction  <= `c1' | CNAfraction  >= `c9'
local pct =100/8
forvalues i = 1/7 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAfraction , centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7')
forvalues i = 2/8 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAfraction  , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9')
regress lfpc _S*
estimates stats

* Let's try 10 knots
centile CNAfraction , centile(1, 99)
local c1 = r(c_1)
local c10 = r(c_2)
preserve
drop if CNAfraction  <= `c1' | CNAfraction  >= `c10'
local pct =100/9
forvalues i = 1/8 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAfraction , centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7' , `pct8')
forvalues i = 2/9 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAfraction  , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9' `c10')
regress lfpc _S*
estimates stats

* Let's try 11 knots
centile CNAfraction , centile(1, 99)
local c1 = r(c_1)
local c11 = r(c_2)
preserve
drop if CNAfraction  <= `c1' | CNAfraction  >= `c11'
local pct =100/10
forvalues i = 1/9 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAfraction , centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7' , `pct8' , `pct9')
forvalues i = 2/10 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAfraction  , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9' `c10' `c11')
regress lfpc _S*
estimates stats


* 10 knots gives the best fit. Let's run it again
centile CNAfraction , centile(1, 99)
local c1 = r(c_1)
local c10 = r(c_2)
preserve
drop if CNAfraction  <= `c1' | CNAfraction  >= `c10'
local pct =100/9
forvalues i = 1/8 {
    local pct`i' = `pct' +(`i'-1)*`pct'
}
centile CNAfraction , centile(`pct1', `pct2', `pct3', `pct4', `pct5', `pct6', `pct7' , `pct8')
forvalues i = 2/9 {
local j =`i'-1
    local c`i' =r(c_`j')
}
restore
drop _S*
mkspline _S = CNAfraction  , cubic displayknots knots(`c1'  `c2'  `c3'  `c4'  `c5'  `c6'  `c7'  `c8' `c9' `c10')
regress lfpc _S*
estimates stats

predict yhat, xb
replace yhat = . if yhat < 0

predict se, stdp
gen lb = max(yhat - se*1.96, 0)
gen ub = yhat + se*1.96
sum CNAfraction
local 0 = 0
foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(1+`tmp')
    }
}

* Plot regression curve of log FPC vs. CNA fraction
twoway (rarea lb ub CNAfraction , sort fcolor(yellow) lwidth(none))     ///
    (line yhat CNAfraction , sort lcolor(red) lwidth(medthick))	    ///
    , ylabel(`0' "0"  `1' "1" `2' "2" `3' "3"`4' "4" `6' "6"  			    ///
          `10' "10" `20' "20" `30' "30" `40' "40" `60' "60" , angle(zero)) 	    ///
      ymtick(`5' "5" `7' "7" `8' "8" `9' "9" `50' "50") 			    ///
      xtitle(CNA fraction) xlabel(0(.1)1) ytitle( FG count) 		    ///
      legend(off) name(fpccna1)  graphregion(color(white)) ysize(2.8) xsize(2.99) 
graph export fpc_vs_cna_graph.pdf, name(fpccna1) replace            

    
log close
translate fusionpairs_vs__CNAfraction_heatmap.log fusionpairs_vs__CNAfraction_heatmap.pdf, replace
