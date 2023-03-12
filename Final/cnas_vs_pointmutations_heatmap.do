log using cnas_vs_pointmutations_heatmap.log , replace
* cnas_vs_pointmutations_heatmap.log
pause on
*This program generates  figure 1B

use ../build.dta ,clear
gen SMcount= PointMutationCount
gen CNAfraction = CNAFractionGenomeAltered
keep SMcount CNAfraction
drop if SMcount ==.
drop if CNAfraction ==.
gen logSMcount = log(SMcount)
codebook  SMcount logSMcount CNAfraction

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

* Select evenly spaced bins for CNAfraction
sum CNAfraction, detail
local maxy = r(max)
local binsy = 60
local binwidthCNA = round(`maxy'/`binsy',0.001)
gen CNAfraction_binned = `binwidthCNA'/2
local binsym1= `binsy'-1
forvalues i = 1/`binsym1' {
    replace CNAfraction_binned = `binwidthCNA'/2 + `i'*`binwidthCNA' if CNAfraction >`i'*`binwidthCNA' 
}

tabulate CNAfraction_binned
histogram CNAfraction_binned, freq width(`binwidth') name(CNA)
sort logSMcount_binned CNAfraction_binned
save tmp, replace
local obs = `binsx' * `binsy'

clear
di "`binsx' `binsy' `obs'"
set obs `obs'
gen logSMcount_binned = .
gen CNAfraction_binned = .
forvalues i =1/`binsx' {
    quietly replace logSMcount_binned = -`binwidthSM'/2 + `i'*`binwidthSM' if `i' == 1 + int((_n-1)/`binsy')
    forvalues j = 1/`binsy' {
        quietly replace CNAfraction_binned = -`binwidthCNA'/2 + `binwidthCNA'*`j' if (`j' == mod(_n,`binsy') | mod(_n,`binsy') ==0)
    }
}
merge 1:m logSMcount_binned CNAfraction_binned using tmp
gen observation = _merge==3
collapse (sum) count = observation , by ( logSMcount_binned CNAfraction_binned)



gen logcount = log(count)
sum logcount
replace logcount = -0.5 if count == 0
twoway (contour logcount CNAfraction_binned logSMcount_binned, heatmap ccuts(-.5 (.5) 7.5)  ccolors(blue "0 59 255" 	///
    "0 120 255" cyan*1.4 cyan cyan*.6 green*1.8 green*1.4 green  )   scolor(blue) ecolor(red)	name(heat1))        ///

sum logcount
local inc = r(max)/15
di `inc'
replace logcount = -`inc' if count == 0


local ivalues "`inc' "
* N.B. Skip the third cut value because the 3rd and 4th value corespond to a count of 2
forvalues i =2/15 {
    local nexti = `i'*`inc' + `inc'
    local ivalues "`ivalues' `nexti'"
}
di "`ivalues'"

* But we want the key heights to be equal
local ivalues2 "`inc' "
forvalues i =1/14 {
    local nexti = `i'*`inc' + `inc'
    local ivalues2 "`ivalues2' `nexti'"
}


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
di "`allcuts'"
di "`cut1'"

foreach i of numlist  1 10 100 1000 10000 {
    forvalues j = 1/9 {
        local tmp = `i'*`j'
        local `tmp'  = log(`tmp')
    }
}

twoway (contour logcount CNAfraction_binned logSMcount_binned, 	///
    heatmap ccuts(`allcuts')  ccolors(white white blue blue*.5 	///
    yellow cyan*.5 cyan cyan*1.5 green green*1.5    	///
     magenta*.5 magenta magenta*1.5  yellow*1.5 orange red red*1.5)	///   
    zlabel( -`inc' "    0" 0 "    1" ///
        `cut1'  "    2"   ///
        `cut2'  "    3" ///  
        `cut3'  "    4" ///  
        `cut4'  "    5" ///  
        `cut5'  "    7" ///  
        `cut6'  "    9" ///  
        `cut7'  "  12" ///  
        `cut8'  "  16" ///  
        `cut9'  "  22" ///  
        `cut10' "  30" ///  
        `cut11' "  41" ///  
        `cut12' "  55" ///  
        `cut13' "  75" ///  
        `cut14' "100"  ///
        `cut15' "138" )  ///
     xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 			///
        `10' "10"  `30' "30"  `60' "60"				///
        `100' "100"  `300' "300"  `600' "600"			///
        `1000' "1000"  `3000' "3000"  `6000' "6000" 		///
        `10000' "10,000" `20000' "20,000" , angle(45))	       	///
    xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 		///
        `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" 	///
        `90' "90" `200' "200" `400' "400" `500' "500" 		///
        `700' "700" `800' "800" `900' "900" `2000' "2000" 	///
        `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" ///
        `9000' "9000") 						///
    ztitle(Number of patients) xtitle(SM count) ytitle(CNA fraction) ///
    ylabel(0(.1)1, angle(0)) name(heat2) graphregion(color(white)) ysize(2.8) xsize(3.75) ) 

use tmp, clear 
gen lpmc = log(SMcount)
*twoway (scatter CNAFractionGenomeAltered lpmc)

regress CNAfraction lpmc
estimates stats

* Find the optimal number of cubic spline knots
mkspline _S = lpmc, cubic nknots(3)
regress CNAfraction _S*
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(4)
regress CNAfraction _S*
test _S2 _S3 
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(5)
regress CNAfraction _S*
test _S2 _S3 _S4
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(6)
regress CNAfraction _S*
test _S2 _S3 _S4 _S5
estimates stats
drop _S*
mkspline _S = lpmc, cubic nknots(7)
regress CNAfraction _S*
estimates stats
test _S2 _S3 _S4 _S5 _S6

* 6 knots give best fit.
drop _S*
mkspline _S = lpmc, cubic nknots(6)
regress CNAfraction _S*
test _S2 _S3 _S4 _S5

predict yhat, xb
replace yhat = . if yhat <0 & lpmc > log(10000)
replace yhat = 0 if yhat < 0
predict se, stdp
gen lb = max(yhat - se*1.96, 0)
gen ub = yhat + se*1.96
sum SMcount
 
twoway    (rarea lb ub lpmc, sort fcolor(yellow) lwidth(none)) 	///    
    (line yhat lpmc, sort lcolor(red) lwidth(medthick)),	///
     xlabel(`1' "1" `2' "2" `3' "3"`6' "6" 			///
        `10' "10"  `30' "30"  `60' "60"				///
        `100' "100"  `300' "300"  `600' "600"			///
        `1000' "1000"  `3000' "3000"  `6000' "6000" 		///
        `10000' "10,000" `20000' "20,000" , angle(45))	       	///
    xmtick(`4' "4" `5' "5" `7' "7" `8' "8" `9' "9" 		///
        `20' "20" `40' "40" `50' "50" `70' "70" `80' "80" 	///
        `90' "90" `200' "200" `400' "400" `500' "500" 		///
        `700' "700" `800' "800" `900' "900" `2000' "2000" 	///
        `4000' "4000" `5000' "5000" `7000' "7000" `8000' "8000" ///
        `9000' "9000") 						///
    xtitle(SM count) ytitle(CNA fraction) 			///
    ylabel(0(.1)1, angle(0)) name(regression) 			///
    legend(order(1 "95% CI" 2 "Expected") ring(0) position(1) col(1)) 	///
    graphregion(color(white)) ysize(2.8) xsize(3.1) 
 
    
log close
translate cnas_vs_pointmutations_heatmap.log cnas_vs_pointmutations_heatmap.pdf, replace
