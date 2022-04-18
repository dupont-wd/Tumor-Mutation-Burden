log using modelFit.c2.log, replace
* ModelFit.log

use cross_validated.dta
stset follow15, failure(fate15) scale(1)

*stcox loghr, nohr
stcox loghr_cv, nohr


centile loghr_cv , centile(10(10) 90)
local c1 = r(c_1)
local c2 = r(c_2)
local c3 = r(c_3)
local c4 = r(c_4)
local c5 = r(c_5)
local c6 = r(c_6)
local c7 = r(c_7)
local c8 = r(c_8)
local c9 = r(c_9)
di exp(`c1')
di exp(`c2')
di exp(`c3')
di exp(`c4')
di exp(`c5')
di exp(`c6')
di exp(`c7')
di exp(`c8')
di exp(`c9')

sum loghr_cv
local max = r(max)
di exp(`max')

gen loghr_cv10 = recode(loghr_cv, `c1', `c2', `c3', `c4' , `c5', `c6', 	///
    `c7', `c8', `c9', `max')
egen decile = group(loghr_cv10)
tabulate decile
*sts test decile, logrank
*sts test quartile if quartile <3, logrank
*sts test quartile if quartile != 1 & quartile < 4, logrank
*sts test quartile if quartile > 2 & quartile != ., logrank

* Draw a survival plot associated with loghr_cv deciles with legend
sts graph, by(decile) 							///
    plot1opts(lcolor(blue) lwidth(medthick)) plot2opts(lcolor(red) lwidth(medthick)) ///
    plot3opts(lcolor(dkgreen) lwidth(medthick))				/// 
    plot4opts(lcolor(cyan) lwidth(medthick)) 				///
    plot5opts(lcolor(green) lwidth(medthick)) 				///
    plot6opts(lcolor(gold) lwidth(medthick)) 				///
    plot7opts(lcolor(magenta) lwidth(medthick))				///
    plot8opts(lcolor(orange) lwidth(medthick)) 				///
    plot9opts(lcolor(brown) lwidth(medthick)) 				///
    plot10opts(lcolor(black) lwidth(medthick)) 				///
    ytitle(Probability of overall survival) 	 			///
    ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
    title("", size(zero)) 							///
    legend( order(1 "1st     1.00 {&minus} 5.93 " 2 "2nd    5.94 {&minus} 8.60" ///
        3 "3rd     8.61 {&minus} 10.4" 4 "4th     10.5 {&minus} 12.1" 		///
        5 "5th     12.2 {&minus} 13.4" ) subtitle(Decile            Hazard ratio)	///	
        cols(1) position(5) ring(1)  ) graphregion(color(white)) 		///symxsize(*.5)
    subtitle("Kaplan-Meier curves by decile of cross-validated hazard ratios"	///
        , ring(0) position(1)) name(decile_rqs1)
graph export overall_survival_deciles_legend.pdf, name(decile_rqs1) replace     

* Draw a survival plot associated with loghr_cv deciles with legend
sts graph, by(decile) 							///
    plot1opts(lcolor(blue) lwidth(medthick)) plot2opts(lcolor(red) lwidth(medthick)) ///
    plot3opts(lcolor(dkgreen) lwidth(medthick))				/// 
    plot4opts(lcolor(cyan) lwidth(medthick)) 				///
    plot5opts(lcolor(green) lwidth(medthick)) 				///
    plot6opts(lcolor(gold) lwidth(medthick)) 				///
    plot7opts(lcolor(magenta) lwidth(medthick))				///
    plot8opts(lcolor(orange) lwidth(medthick)) 				///
    plot9opts(lcolor(brown) lwidth(medthick)) 				///
    plot10opts(lcolor(black) lwidth(medthick)) 				///
    ytitle(Probability of overall survival) 	 			///
    ylabel(0(.1)1, angle(zero)) xtitle(Years since diagnosis) xlabel(0(2)14) 	///
    title("", size(zero)) 							///
    legend( off ) graphregion(color(white)) 	///
    subtitle("Kaplan-Meier curves by decile of cross-validated hazard ratios"	///
        , ring(0) position(1)) name(decile_rqs2)
graph export overall_survival_deciles_nolegend.pdf, name(decile_rqs2) replace     
        
    
centile loghr_cv, centile(5(10) 95)
local c1 = r(c_1)
local c2 = r(c_2)
local c3 = r(c_3)
local c4 = r(c_4)
local c5 = r(c_5)
local c6 = r(c_6)
local c7 = r(c_7)
local c8 = r(c_8)
local c9 = r(c_9)
local c10 = r(c_10)

* Draw survival plots under the proportional hazards model with legend
stcurve , survival at( loghr = (`c1' `c2' `c3' `c4' `c5' `c6' `c7' `c8'  	///
	`c9' `c10' ) )  ylabels(0 (.2)1) 					///
    lcolor(blue red dkgreen cyan green gold magenta orange brown black) 	///
    lwidth(medthick medthick medthick medthick medthick medthick medthick 	///
        medthick medthick medthick) ytitle(Probability of overall survival) 	///
    ylabel(0(.1)1, angle(zero) format(%4.2f)) xtitle(Years since diagnosis) 	///
    xlabel(0(2)14) title("", size(zero)) 					///
    legend( order( 6 "6th     13.5 {&minus} 14.6" 7 "7th     14.7 {&minus} 15.8" ///
        8 "8th     15.9 {&minus} 17.0" 9 "9th     17.1 {&minus} 18.3" 		///
        10 "10th   18.4 {&minus} 19.8" ) subtitle(Decile              Hazard ratio)	///  
    cols(1) position(7) ring(1)  ) graphregion(color(white)) 		/// symxsize(*.5)
    subtitle("Survival curves under proportional hazards model"  		///
        "by decile of cross-validated hazard ratios", ring(0) position(1))name(ph_fit)
graph export overall_survival_ph_fit_legend.pdf, name(ph_fit) replace   

* Draw survival plots under the proportional hazards model without legend
stcurve , survival at( loghr = (`c1' `c2' `c3' `c4' `c5' `c6' `c7' `c8'  	///
	`c9' `c10' ) )  ylabels(0 (.2)1) 					///
    lcolor(blue red dkgreen cyan green gold magenta orange brown black) 	///
    lwidth(medthick medthick medthick medthick medthick medthick medthick 	///
        medthick medthick medthick) ytitle(Probability of overall survival) 	///
    ylabel(0(.1)1, angle(zero) format(%4.2f)) xtitle(Years since diagnosis) 	///
    xlabel(0(2)14) title("", size(zero)) 					///
    legend( off )  graphregion(color(white))	///
    subtitle("Survival curves under proportional hazards model"  		///
        "by decile of cross-validated hazard ratios", ring(0) position(1))name(ph_fit2)
graph export overall_survival_ph_fit_lolegend.pdf, name(ph_fit2) replace     


log close
translate modelFit.c2.log modelFit.c2.pdf, replace
