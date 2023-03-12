log using fusionpairs_vs_pointmutations_heatmap.log , replace
* fusionpairs_vs_pointmutations_heatmap.log
pause on
*This program generates  figure 1D

use ../build.dta ,clear
gen SMcount= PointMutationCount
gen FPcount = FusionPairs
keep SMcount FPcount
drop if SMcount ==.
drop if FPcount ==.
gen logSMcount = log(SMcount)
gen logFPcount = log(FPcount +1)
codebook  SMcount logSMcount FPcount logFPcount

* Select evenly spaced bins for logSMcount
sum logSMcount, detail
local maxx = r(max)
local binsx = 80
local binwidthSM = round(`maxx'/`binsx',0.001)
gen logSMcount_binned = `binwidthSM'/2
local binsxm1= `binsx'-1
forvalues i = 1/`binsxm1' {
    replace logSMcount_binned = `binwidthSM'/2 + `i'*`binwidthSM' if logSMcount >`i'*`binwidthSM' 
}
tabulate logSMcount_binned
histogram logSMcount_binned, freq width(`binwidth') name(SM)

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
histogram logFPcount_binned, freq width(`binwidth') name(CNA)
sort logSMcount_binned logFPcount_binned
save tmp, replace
local obs = `binsx' * `binsy'

clear
di "`binsx' `binsy' `obs'"
set obs `obs'
gen logSMcount_binned = .
gen logFPcount_binned = .
forvalues i =1/`binsx' {
    quietly replace logSMcount_binned = -`binwidthSM'/2 + `i'*`binwidthSM' if `i' == 1 + int((_n-1)/`binsy')
    forvalues j = 1/`binsy' {
        quietly replace logFPcount_binned = -`binwidthFP'/2 + `binwidthFP'*`j' if (`j' == mod(_n,`binsy') | mod(_n,`binsy') ==0)
    }
}
merge 1:m logSMcount_binned logFPcount_binned using tmp
gen observation = _merge==3
collapse (sum) count = observation , by ( logSMcount_binned logFPcount_binned)

gen logcount = log(count)
sum logcount
replace logcount = -0.5 if count == 0
twoway (contour logcount logFPcount_binned logSMcount_binned, heatmap ccuts(-.5 (.5) 7.5)  ccolors(blue "0 59 255" 	///
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
* increase cut15 slightly to account to the fact that there is a cell with 232 patients
local cut15 = log(233)
di "`cut15'"
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
twoway (contour logcount logFPcount_binned logSMcount_binned 	///
  , heatmap ccuts(`allcuts')  ccolors(white white blue blue*.5 	///
    yellow yellow cyan cyan*1.5 green green*1.5    	///
     magenta*.5 magenta magenta*1.5  yellow*1.5 orange red red*1.5)	///   
    zlabel( -`inc' "    0" 0 "    1" ///
        `cut1'  "    2"   ///
        `cut2'  "    3" ///  
        `cut3'  "    3" ///  
        `cut4'  "    5" ///  
        `cut5'  "    7" ///  
        `cut6'  "    9" ///  
        `cut7'  "  13" ///  
        `cut8'  "  19" ///  
        `cut9'  "  27" ///  
        `cut10' "  38" ///  
        `cut11' "  55" ///  
        `cut12' "  79" ///  
        `cut13' " 113" ///  
        `cut14' "162"  ///
        `cut15' "232" )  ///
     ylabel(`y0' "0"  `y1' "1" `y2' "2" `y3' "3"`y4' "4" `y6' "6"  ///
         `y10' "10" `y20' "20" `y30' "30" `y40' "40" `y60' "60" , angle(zero)) ///
     ymtick(`y5'  `y7'  `y8' `y9' `y50')			///
     xlabel(`x1' "1" `x2' "2" `x3' "3"`x6' "6" 			///
        `x10' "10"  `x30' "30"  `x60' "60"				///
        `x100' "100"  `x300' "300"  `x600' "600"			///
        `x1000' "1000"  `x3000' "3000"  `x6000' "6000" 		///
        `x10000' "10,000" `x20000' "20,000" , angle(45))	       	///
    xmtick(`x4'  `x5'  `x7'  `x8'  `x9'  		///
           `x20' `x40'   `x50'   `x70'   `x80'  	///
           `x90' `x200'  `x400'  `x500'  		///
          `x700' `x800'  `x900'  `x2000'  	///
        `x4000'  `x5000' `x7000' `x8000'  ///
        `x9000'   ) 						///
    ztitle(Number of patients) xtitle(SM count) ytitle(FG count) ///
    name(heat2) graphregion(color(white)) ysize(2.8) xsize(3.75) ) 

use tmp, clear 
gen lpmc = log(SMcount)
*twoway (scatter CNAFractionGenomeAltered lpmc)

regress logFPcount lpmc
estimates stats

* Find the optimal number of cubic spline knots
mkspline _S = lpmc, cubic nknots(3)
regress logFPcount _S*
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(4)
regress logFPcount _S*
test _S2 _S3 
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(5)
regress logFPcount _S*
test _S2 _S3 _S4
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(6)
regress logFPcount _S*
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(7)
regress logFPcount _S*
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
regress logFPcount _S*
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
regress logFPcount _S*
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
regress logFPcount _S*

predict yhat if lpmc >= 0 & lpmc != ., xb
replace yhat = . if yhat < 0
predict se, stdp
gen lb = max(yhat - se*1.96, 0) if SMcount < 20000
gen ub = yhat + se*1.96 if SMcount < 20000
sum SMcount
 
twoway    (rarea lb ub lpmc, sort fcolor(yellow) lwidth(none)) 	///    
    (line yhat lpmc, sort lcolor(red) lwidth(medthick)),	///
    xlabel(`x1' "1" `x2' "2" `x3' "3"`x6' "6" 			///
        `x10' "10"  `x30' "30"  `x60' "60"				///
        `x100' "100"  `x300' "300"  `x600' "600"			///
        `x1000' "1000"  `x3000' "3000"  `x6000' "6000" 		///
        `x10000' "10,000" `x20000' "20,000" , angle(45))	       	///
    xmtick(`x4'  `x5'  `x7'  `x8'  `x9'  		///
        `x20' `x40'   `x50'   `x70'   `x80'  	///
        `x90' `x200'  `x400'  `x500'  		///
        `x700' `x800'  `x900'  `x2000'  	///
        `x4000'  `x5000' `x7000' `x8000'  ///
        `x9000'   ) 						///
      ylabel(`y0' "0"  `y1' "1" `y2' "2" `y3' "3"`y4' "4" `y6' "6"  	///
        `y10' "10" `y20' "20" `y30' "30" `y40' "40" `y60' "60" , angle(zero)) ///  
    ymtick(`y5'  `y7'  `y8' `y9' `y50')			///
    xtitle(SM count) ytitle(FG count) name(regression) legend(off) ///
    graphregion(color(white)) ysize(2.8) xsize(2.99) 
 
    
log close
translate fusionpairs_vs_pointmutations_heatmap.log fusionpairs_vs_pointmutations_heatmap.pdf, replace
