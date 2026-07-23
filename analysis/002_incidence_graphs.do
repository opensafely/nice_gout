version 16

/*==============================================================================
DO FILE NAME:			Incidence graphs
PROJECT:				OpenSAFELY NICE
DATE: 					18/07/2025
AUTHOR:					M Russell									
DESCRIPTION OF FILE:	Incidence tables and graphs
DATASETS USED:			Incidence and Measures files
OTHER OUTPUT: 			logfiles, printed to folder $Logdir
USER-INSTALLED ADO: 	 
  (place .ado file(s) in analysis folder)						
==============================================================================*/

*Set filepaths
/*
global projectdir "C:/Users/k1754142/OneDrive/PhD Project/OpenSAFELY NICE/nice_gout"
global running_locally = 1 // Running on local machine
*/

global projectdir `c(pwd)'
global running_locally = 0 // Running on OpenSAFELY console

capture mkdir "$projectdir/output/data"
capture mkdir "$projectdir/output/tables"
capture mkdir "$projectdir/output/figures"

*Open a log file
global logdir "$projectdir/logs"
cap log close
log using "$logdir/incidence_graphs.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list and intervention date (passed from yaml)
global arglist diseases intervention_date_1
args $arglist

if $running_locally ==0 {
	foreach var of global arglist {
		local `var' : subinstr local `var' "|" " ", all
		global `var' "``var''"
		di "$`var'"
	}
}

if $running_locally ==1 {
	global diseases "gout"
	global intervention_date_1 "2020-03-01"
}

di "$diseases"
di "$intervention_date_1"

set type double

set scheme plotplainblind

*Import rounded and redacted data for each disease ==============================================*/
local first_disease: word 1 of $diseases
di "`first_disease'"

**Import first file as base dataset
import delimited "$projectdir/output/tables/redacted_counts_`first_disease'.csv", clear
save "$projectdir/output/data/redacted_standardised.dta", replace

**Loop over diseases and years
foreach disease in $diseases {
	if "`disease'" != "`first_disease'"  {
		import delimited "$projectdir/output/tables/redacted_counts_`disease'.csv", clear
		append using "$projectdir/output/data/redacted_standardised.dta"
		save "$projectdir/output/data/redacted_standardised.dta", replace 
		}
}

*Format data =================================================================*/
use "$projectdir/output/data/redacted_standardised.dta", clear

gen mo_year_diagn_s = monthly(mo_year_diagn, "MY") 
format mo_year_diagn_s %tmMon-CCYY
drop mo_year_diagn
rename mo_year_diagn_s mo_year_diagn
order mo_year_diagn, after(measure)
gen year=yofd(dofm(mo_year_diagn))
order year, after(mo_year_diagn)

***Collapse age bands to 20-year bands (this is for 18+ age groups)
bys disease mo_year_diagn measure: egen numerator_18_39 = sum(numerator_18_29 + numerator_30_39)
bys disease mo_year_diagn measure: egen numerator_40_59 = sum(numerator_40_49 + numerator_50_59)
bys disease mo_year_diagn measure: egen numerator_60_79 = sum(numerator_60_69 + numerator_70_79)

bys disease mo_year_diagn measure: egen denominator_18_39 = sum(denominator_18_29 + denominator_30_39)
bys disease mo_year_diagn measure: egen denominator_40_59 = sum(denominator_40_49 + denominator_50_59)
bys disease mo_year_diagn measure: egen denominator_60_79 = sum(denominator_60_69 + denominator_70_79)

gen rate_18_39 = (numerator_18_39/denominator_18_39)*100000
gen rate_40_59 = (numerator_40_59/denominator_40_59)*100000
gen rate_60_79 = (numerator_60_79/denominator_60_79)*100000

**Optional: for demographic subgroups with >50% missing data, convert all to zero
local strat_rates rate_18_39 rate_40_59 rate_60_79 rate_80 rate_18_29 rate_30_39 rate_40_49 rate_50_59 rate_60_69 rate_70_79 ///
	s_rate_male s_rate_female ///
	rate_white rate_mixed rate_black rate_asian rate_other rate_ethunk ///
	rate_imd1 rate_imd2 rate_imd3 rate_imd4 rate_imd5 rate_imdunk ///
	rate_east rate_eastmid rate_london rate_northeast rate_northwest rate_southeast rate_southwest rate_westmid rate_yorkshire rate_regunk

levelsof measure, local(measure_list)

foreach disease in $diseases {
	foreach measure of local measure_list {
		foreach var of local strat_rates {
			
			quietly count if disease == "`disease'" & measure == "`measure'"
			local total = r(N)
			
			quietly count if missing(`var') & disease == "`disease'" & measure == "`measure'"
			local num_missing = r(N)
			
			if `total' > 0 {
				local pct_missing = (`num_missing' / `total') * 100
				
				if `pct_missing' > 50 {
					replace `var' = 0 if disease == "`disease'" & measure == "`measure'"
				}
			}
		}
	}
}

**Generate 3-monthly moving averages
local rate_vars s_rate_all rate_all ///
	s_rate_male s_rate_female ///
	rate_18_39 rate_40_59 rate_60_79 rate_80 rate_18_29 rate_30_39 rate_40_49 rate_50_59 rate_60_69 rate_70_79 ///
	rate_white rate_mixed rate_black rate_asian rate_other rate_ethunk ///
	rate_imd1 rate_imd2 rate_imd3 rate_imd4 rate_imd5 rate_imdunk ///
	rate_east rate_eastmid rate_london rate_northeast rate_northwest rate_southeast rate_southwest rate_westmid rate_yorkshire rate_regunk
	
foreach var of local rate_vars {
	bysort disease measure (mo_year_diagn): ///
		gen `var'_ma = (`var'[_n-1] + `var'[_n] + `var'[_n+1]) / 3
}
save "$projectdir/output/data/redacted_standardised.dta", replace

**Save subset of data for use with forecasting (next step)
preserve
keep if measure=="Incidence"
rename s_rate_all incidence
rename numerator_all numerator
rename denominator_all denominator
keep disease disease_full mo_year_diagn numerator denominator incidence
outsheet * using "$projectdir/output/tables/arima_standardised.csv", comma replace
export delimited using "$projectdir/output/tables/arima_standardised.csv", datafmt replace
restore

*Create graphs of incidence and prevalence rates diagnoses by month, by disease, using rounded/redacted data ========================================================*/

use "$projectdir/output/data/redacted_standardised.dta", clear

levelsof disease, local(disease_list)

foreach disease of local disease_list {
	
	di "`disease'"
	keep if disease=="`disease'"
	replace disease_full = lower(disease_full)
	local disease_full = disease_full[1]
	
	**Label axis titles
	local ytitle "Monthly incidence rate per 100,000 population"
	local ytitleprev "Annual prevalence per 100,000 population"
	local xtitle ""
			
	**Incidence graphs (by month)
	preserve
		keep if measure=="Incidence"
		
		**Set y-axis ranges graphs
		egen rate_max_all = max(s_rate_all)
		egen rate_min_all = min(s_rate_all)
		egen rate_max_sex = max(max(s_rate_male, s_rate_female))
		egen rate_min_sex = min(min(s_rate_male, s_rate_female))
		egen rate_max_age = max(max(rate_18_39, rate_40_59, rate_60_79, rate_80))
		egen rate_min_age = min(min(rate_18_39, rate_40_59, rate_60_79, rate_80))
		egen rate_max_ethn = max(max(rate_white, rate_mixed, rate_black, rate_asian, rate_other, rate_ethunk))
		egen rate_min_ethn = min(min(rate_white, rate_mixed, rate_black, rate_asian, rate_other, rate_ethunk))
		egen rate_max_imd = max(max(rate_imd1, rate_imd2, rate_imd3, rate_imd4, rate_imd5, rate_imdunk))
		egen rate_min_imd = min(min(rate_imd1, rate_imd2, rate_imd3, rate_imd4, rate_imd5, rate_imdunk))
		egen rate_max_region = max(max(rate_east, rate_eastmid, rate_london, rate_northeast, rate_northwest, rate_southeast, rate_southwest, rate_westmid, rate_yorkshire, rate_regunk))
		egen rate_min_region = min(min(rate_east, rate_eastmid, rate_london, rate_northeast, rate_northwest, rate_southeast, rate_southwest, rate_westmid, rate_yorkshire, rate_regunk))

		foreach stem in all sex age ethn imd region {

			quietly summ rate_min_`stem', meanonly
			local rmin = r(min)
			quietly summ rate_max_`stem', meanonly
			local rmax = r(max)

			if missing(`rmin') | missing(`rmax') {
				di as error "Rates not estimable for `stem'; using default y-axis"
				local lower_`stem' = 0
				local upper_`stem' = 1
			}
			else if `rmax' == 0 {
				local lower_`stem' = 0
				local upper_`stem' = 1
			}
			else if `rmax' < 1 {
				local lower_`stem' = round(0.80 * `rmin', 0.01)
				local upper_`stem' = round(1.10 * `rmax', 0.01)
			}
			else if `rmax' < 10 {
				local lower_`stem' = round(0.80 * `rmin', 0.1)
				local upper_`stem' = round(1.10 * `rmax', 0.1)
			}
			else if `rmax' < 100 {
				local lower_`stem' = round(0.80 * `rmin', 1)
				local upper_`stem' = round(1.10 * `rmax', 1)
			}
			else {
				local lower_`stem' = round(0.80 * `rmin', 10)
				local upper_`stem' = round(1.10 * `rmax', 10)
			}

			if `lower_`stem'' < 0.2 local lower_`stem' = 0
			if `upper_`stem'' <= `lower_`stem'' local upper_`stem' = `lower_`stem'' + 1

			if `upper_`stem'' < 2 {
				local format_`stem' = "format(%9.1f)"
			}
			else {
				local format_`stem' = "format(%9.0f)"
			}

			di `lower_`stem''
			di `upper_`stem''

			nicelabels `lower_`stem'' `upper_`stem'', local(ylab_`stem')
			di "`ylab_`stem''"
		}
		
		***Set x-axis format for incidence graphs
		quietly summarize mo_year_diagn, meanonly
		local lb = r(min)
		local ub = r(max)

		local y0 = yofd(dofm(`lb'))
		local y1 = yofd(dofm(`ub'))

		local step = 1 //number of year gaps on axis
		local ystart = cond(mod(`y0',`step'), `y0' + (`step' - mod(`y0',`step')), `y0')
		local yfinish = `y1' + `step'

		local xlabel ""
		forvalues y = `ystart'(`step')`yfinish' {
			local m = ym(`y',1)
			local xlabel `xlabel' `m' "`y'"
		}
		
		***Set vertical intervention line
		local year  = real(substr("$intervention_date_1", 1, 4))
		local month = real(substr("$intervention_date_1", 6, 2))
		local intervention = ym(`year', `month')
		di `intervention'
		local xintlab = `intervention' + 1 
		
		*Adjusted incidence overall (scatter with moving average)
		twoway scatter s_rate_all mo_year_diagn, ytitle("`ytitle'", size(medsmall)) color(emerald%20) msymbol(circle) || line s_rate_all_ma mo_year_diagn, lcolor(emerald) lstyle(solid) ylabel(`ylab_all', `format_all' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(off) xsize(16) ysize(9) name("inc_adj_`disease'", replace) saving("$projectdir/output/figures/inc_`disease'.gph", replace)
			*graph export "$projectdir/output/figures/inc_`disease'.png", replace
			graph export "$projectdir/output/figures/inc_`disease'.svg", replace
	
		*Adjusted incidence by sex (scatter with moving average)
		twoway scatter s_rate_male mo_year_diagn, ytitle("`ytitle'", size(medsmall)) color(eltblue%20) mlcolor(eltblue%20) msymbol(circle) || line s_rate_male_ma mo_year_diagn, lcolor(midblue) lstyle(solid) || scatter s_rate_female mo_year_diagn, color(orange%20) mlcolor(orange%20) msymbol(circle)  || line s_rate_female_ma mo_year_diagn, lcolor(red) lstyle(solid) ylabel(`ylab_sex', `format_sex' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) order(2 "Male" 4 "Female")) xsize(16) ysize(9) name(inc_`disease'_sex, replace) saving("$projectdir/output/figures/inc_`disease'_sex.gph", replace)
			*graph export "$projectdir/output/figures/inc_`disease'_sex.png", replace
			graph export "$projectdir/output/figures/inc_`disease'_sex.svg", replace

		*Adjusted incidence comparison (moving average)
		twoway line rate_all_ma mo_year_diagn, ytitle("`ytitle'", size(medsmall)) lstyle(solid) lcolor(gold)  || line s_rate_all_ma mo_year_diagn, lstyle(solid) lcolor(emerald) ylabel(`ylab_all', `format_all' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) xline(`intervention') title("", size(medium) margin(b=2)) legend(region(fcolor(white%0)) order(1 "Crude" 2 "Adjusted")) xsize(16) ysize(9) name(inc_`disease'_comp, replace) saving("$projectdir/output/figures/inc_`disease'_comp.gph", replace)
			*graph export "$projectdir/output/figures/inc_`disease'_comp.png", replace
			graph export "$projectdir/output/figures/inc_`disease'_comp.svg", replace 
		
		*Unadjusted incidence by 20-year age groups (moving average)
		twoway line rate_18_39_ma mo_year_diagn, lcolor(ltblue) lstyle(solid) ytitle("`ytitle'", size(medsmall)) || line rate_40_59_ma mo_year_diagn, lcolor(ebblue) lstyle(solid) || line rate_60_79_ma mo_year_diagn, lcolor(blue) lstyle(solid) || line rate_80_ma mo_year_diagn, lcolor(navy) lstyle(solid) ylabel(`ylab_age', `format_age' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("Age group", size(small) margin(b=1)) order(1 "18 to 39" 2 "40 to 59" 3 "60 to 79" 4 "80 or above")) xsize(16) ysize(9) name(inc_`disease'_agegroup, replace) saving("$projectdir/output/figures/inc_`disease'_agegroup.gph", replace)
			*graph export "$projectdir/output/figures/inc_`disease'_agegroup.png", replace
			graph export "$projectdir/output/figures/inc_`disease'_agegroup.svg", replace	
		
		*Unadjusted incidence by IMD (moving average)
		twoway line rate_imd1_ma mo_year_diagn, lcolor(ltblue) lstyle(solid) ytitle("`ytitle'", size(medsmall)) || line rate_imd2_ma mo_year_diagn, lcolor(eltblue) lstyle(solid) || line rate_imd3_ma mo_year_diagn, lcolor(ebblue) lstyle(solid) || line rate_imd4_ma mo_year_diagn, lcolor(blue) lstyle(solid) || line rate_imd5_ma mo_year_diagn, lcolor(navy) lstyle(solid) ylabel(`ylab_imd', `format_imd' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("IMD quintile", size(small) margin(b=1)) order(1 "1 Most deprived" 2 "2" 3 "3" 4 "4" 5 "5 Least deprived")) xsize(16) ysize(9) name(inc_`disease'_imd, replace) saving("$projectdir/output/figures/inc_`disease'_imd.gph", replace)
			*graph export "$projectdir/output/figures/inc_`disease'_imd.png", replace
			graph export "$projectdir/output/figures/inc_`disease'_imd.svg", replace
		
		*Unadjusted incidence by ethnicity (moving average)
		twoway line rate_white_ma mo_year_diagn, lcolor(ltblue) lstyle(solid) ytitle("`ytitle'", size(medsmall)) || line rate_mixed_ma mo_year_diagn, lcolor(eltblue) lstyle(solid) || line rate_black_ma mo_year_diagn, lcolor(ebblue) lstyle(solid) || line rate_asian_ma mo_year_diagn, lcolor(blue) lstyle(solid) || line rate_other_ma mo_year_diagn, lcolor(navy) lstyle(solid) ylabel(`ylab_ethn', `format_ethn' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("Ethnicity", size(small) margin(b=1)) order(1 "White" 2 "Mixed" 3 "Black" 4 "Asian" 5 "Chinese or other")) name(inc_`disease'_ethnicity, replace) xsize(16) ysize(9) saving("$projectdir/output/figures/inc_`disease'_ethnicity.gph", replace)
		*graph export "$projectdir/output/figures/inc_`disease'_ethnicity.png", replace
		graph export "$projectdir/output/figures/inc_`disease'_ethnicity.svg", replace
		
		*Unadjusted incidence by region (moving average)
		twoway line rate_east_ma mo_year_diagn, lcolor(emerald) lstyle(solid) ytitle("`ytitle'", size(medsmall)) || line rate_eastmid_ma mo_year_diagn, lcolor(orange) lstyle(solid) || line rate_london_ma mo_year_diagn, lcolor(red) lstyle(solid) || line rate_northeast_ma mo_year_diagn, lcolor(blue) lstyle(solid) || line rate_northwest_ma mo_year_diagn, lcolor(dkgreen) lstyle(solid) || line rate_southeast_ma mo_year_diagn, lcolor(cranberry) lstyle(solid) || line rate_southwest_ma mo_year_diagn, lcolor(navy) lstyle(solid) || line rate_westmid_ma mo_year_diagn, lcolor(maroon) lstyle(solid) || line rate_yorkshire_ma mo_year_diagn, lcolor(teal) lstyle(solid) ylabel(`ylab_region', `format_region' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("Region", size(small) margin(b=1)) order(1 "East" 2 "East Midlands" 3 "London" 4 "North East" 5 "North West" 6 "South East" 7 "South West" 8 "West Midlands" 9 "Yorkshire Humber")) name(inc_`disease'_region, replace) xsize(16) ysize(9) saving("$projectdir/output/figures/inc_`disease'_region.gph", replace)
		*graph export "$projectdir/output/figures/inc_`disease'_region.png", replace
		graph export "$projectdir/output/figures/inc_`disease'_region.svg", replace
		
	restore		
	
	**Prevalence graphs (by year)
	preserve
		keep if measure=="Prevalence"
		
		**Set y-axis ranges graphs
		egen rate_max_all = max(s_rate_all)
		egen rate_min_all = min(s_rate_all)
		egen rate_max_sex = max(max(s_rate_male, s_rate_female))
		egen rate_min_sex = min(min(s_rate_male, s_rate_female))
		egen rate_max_age = max(max(rate_18_39, rate_40_59, rate_60_79, rate_80))
		egen rate_min_age = min(min(rate_18_39, rate_40_59, rate_60_79, rate_80))
		egen rate_max_ethn = max(max(rate_white, rate_mixed, rate_black, rate_asian, rate_other, rate_ethunk))
		egen rate_min_ethn = min(min(rate_white, rate_mixed, rate_black, rate_asian, rate_other, rate_ethunk))
		egen rate_max_imd = max(max(rate_imd1, rate_imd2, rate_imd3, rate_imd4, rate_imd5, rate_imdunk))
		egen rate_min_imd = min(min(rate_imd1, rate_imd2, rate_imd3, rate_imd4, rate_imd5, rate_imdunk))
		egen rate_max_region = max(max(rate_east, rate_eastmid, rate_london, rate_northeast, rate_northwest, rate_southeast, rate_southwest, rate_westmid, rate_yorkshire, rate_regunk))
		egen rate_min_region = min(min(rate_east, rate_eastmid, rate_london, rate_northeast, rate_northwest, rate_southeast, rate_southwest, rate_westmid, rate_yorkshire, rate_regunk))

		foreach stem in all sex age ethn imd region {

			quietly summ rate_min_`stem', meanonly
			local rmin = r(min)
			quietly summ rate_max_`stem', meanonly
			local rmax = r(max)

			if missing(`rmin') | missing(`rmax') {
				di as error "Rates not estimable for `stem'; using default y-axis"
				local lower_`stem' = 0
				local upper_`stem' = 1
			}
			else if `rmax' == 0 {
				local lower_`stem' = 0
				local upper_`stem' = 1
			}
			else if `rmax' < 1 {
				local lower_`stem' = round(0.80 * `rmin', 0.01)
				local upper_`stem' = round(1.10 * `rmax', 0.01)
			}
			else if `rmax' < 10 {
				local lower_`stem' = round(0.80 * `rmin', 0.1)
				local upper_`stem' = round(1.10 * `rmax', 0.1)
			}
			else if `rmax' < 100 {
				local lower_`stem' = round(0.80 * `rmin', 1)
				local upper_`stem' = round(1.10 * `rmax', 1)
			}
			else {
				local lower_`stem' = round(0.80 * `rmin', 10)
				local upper_`stem' = round(1.10 * `rmax', 10)
			}

			if `lower_`stem'' < 0.2 local lower_`stem' = 0
			if `upper_`stem'' <= `lower_`stem'' local upper_`stem' = `lower_`stem'' + 1

			if `upper_`stem'' < 2 {
				local format_`stem' = "format(%9.1f)"
			}
			else {
				local format_`stem' = "format(%9.0f)"
			}

			di `lower_`stem''
			di `upper_`stem''

			nicelabels `lower_`stem'' `upper_`stem'', local(ylab_`stem')
			di "`ylab_`stem''"
		}
		
		***Set x-axis format for incidence graphs
		quietly summarize year, meanonly
		local lb = floor(r(min))
		local ub = ceil(r(max))
		local step = 1
		local xlabel `lb'(`step')`ub'
		
		***Set vertical intervention line
		local intervention = real(substr("$intervention_date_1",1,4))
		display `intervention'
		
		*Adjusted prevalence overall
		twoway connected s_rate_all year, color(emerald%30) msymbol(circle) lcolor(emerald) lstyle(solid) ytitle("`ytitleprev'", size(medsmall)) ylabel(`ylab_all', `format_all' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(1)) xlabel(`xlabel', nogrid labsize(small) axis(1)) xlabel(`intervention' "COVID-19", axis(2) labsize(small) labcolor(navy)) xtitle("", axis(2)) xscale(noline axis(2)) title("", size(medium) margin(b=2)) xline(`intervention') legend(off) xsize(16) ysize(9) name(prev_`disease', replace) saving("$projectdir/output/figures/prev_`disease'.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'.png", replace
			graph export "$projectdir/output/figures/prev_`disease'.svg", replace
			
		*Adjusted prevalence comparison
		twoway connected rate_all year, color(gold%30) msymbol(circle) lstyle(solid) lcolor(gold) ytitle("`ytitleprev'", size(medsmall)) || connected s_rate_all year, color(emerald%30) msymbol(circle) lstyle(solid) lcolor(emerald) ylabel(`ylab_all', `format_all' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) xline(`intervention') title("", size(medium) margin(b=2)) legend(region(fcolor(white%0)) order(1 "Crude" 2 "Adjusted")) xsize(16) ysize(9) name(prev_`disease'_comp, replace) saving("$projectdir/output/figures/prev_`disease'_comp.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'_comp.png", replace
			graph export "$projectdir/output/figures/prev_`disease'_comp.svg", replace
			
		*Adjusted prevalence by sex
		twoway connected s_rate_male year, color(eltblue%30) msymbol(circle) lcolor(midblue) lstyle(solid) || connected s_rate_female year, color(orange%30) msymbol(circle) lcolor(red) lstyle(solid) ytitle("`ytitleprev'", size(medsmall)) ylabel(`ylab_sex', `format_sex' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) order(1 "Male" 2 "Female")) xsize(16) ysize(9) name(prev_`disease'_sex, replace) saving("$projectdir/output/figures/prev_`disease'_sex.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'_sex.png", replace
			graph export "$projectdir/output/figures/prev_`disease'_sex.svg", replace

		*Unadjusted prevalence by 20-year age groups
		twoway connected rate_18_39 year, color(ltblue%30) msymbol(circle) lcolor(ltblue) lstyle(solid) || connected rate_40_59 year, color(ebblue%30) msymbol(circle) lcolor(ebblue) lstyle(solid) || connected rate_60_79 year, color(blue%30) msymbol(circle) lcolor(blue) lstyle(solid) || connected rate_80 year, color(navy%30) msymbol(circle) lcolor(navy) lstyle(solid) ytitle("`ytitleprev'", size(medsmall)) ylabel(`ylab_age', `format_age' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("Age group", size(small) margin(b=1)) order(1 "18 to 39" 2 "40 to 59" 3 "60 to 79" 4 "80 or above")) xsize(16) ysize(9) name(prev_`disease'_agegroup, replace) saving("$projectdir/output/figures/prev_`disease'_agegroup.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'_agegroup.png", replace
			graph export "$projectdir/output/figures/prev_`disease'_agegroup.svg", replace		
			
		*Unadjusted prevalence by IMD
		twoway connected rate_imd1 year, color(ltblue%30) msymbol(circle) lcolor(ltblue) lstyle(solid) || connected rate_imd2 year, color(eltblue%30) msymbol(circle) lcolor(eltblue) lstyle(solid) || connected rate_imd3 year, color(ebblue%30) msymbol(circle) lcolor(ebblue) lstyle(solid) || connected rate_imd4 year, color(blue%30) msymbol(circle) lcolor(blue) lstyle(solid) || connected rate_imd5 year, color(navy%30) msymbol(circle) lcolor(navy) lstyle(solid) ytitle("`ytitleprev'", size(medsmall)) ylabel(`ylab_imd', `format_imd' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("IMD quintile", size(small) margin(b=1)) order(1 "1 Most deprived" 2 "2" 3 "3" 4 "4" 5 "5 Least deprived")) xsize(16) ysize(9) name(prev_`disease'_imd, replace) saving("$projectdir/output/figures/prev_`disease'_imd.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'_imd.png", replace
			graph export "$projectdir/output/figures/prev_`disease'_imd.svg", replace
			
		*Unadjusted prevalence by ethnicity
		twoway connected rate_white year, color(ltblue%30) msymbol(circle) lcolor(ltblue) lstyle(solid) || connected rate_mixed year, color(eltblue%30) msymbol(circle) lcolor(eltblue) lstyle(solid) || connected rate_black year, color(ebblue%30) msymbol(circle) lcolor(ebblue) lstyle(solid) || connected rate_asian year, color(blue%30) msymbol(circle) lcolor(blue) lstyle(solid) || connected rate_other year, color(navy%30) msymbol(circle) lcolor(navy) lstyle(solid) ytitle("`ytitleprev'", size(medsmall)) ylabel(`ylab_ethn', `format_ethn' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("Ethnicity", size(small) margin(b=1)) order(1 "White" 2 "Mixed" 3 "Black" 4 "Asian" 5 "Chinese or other")) xsize(16) ysize(9) name(prev_`disease'_ethnicity, replace) saving("$projectdir/output/figures/prev_`disease'_ethnicity.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'_ethnicity.png", replace
			graph export "$projectdir/output/figures/prev_`disease'_ethnicity.svg", replace
			
		*Unadjusted prevalence by region
		twoway connected rate_east year, color(emerald%30) msymbol(circle) lcolor(emerald) lstyle(solid) || connected rate_eastmid year, color(orange%30) msymbol(circle) lcolor(orange) lstyle(solid) || connected rate_london year, color(red%30) msymbol(circle) lcolor(red) lstyle(solid) || connected rate_northeast year, color(blue%30) msymbol(circle) lcolor(blue) lstyle(solid) || connected rate_northwest year, color(dkgreen%30) msymbol(circle) lcolor(dkgreen) lstyle(solid) || connected rate_southeast year, color(cranberry%30) msymbol(circle) lcolor(cranberry) lstyle(solid) || connected rate_southwest year, color(navy%30) msymbol(circle) lcolor(navy) lstyle(solid) || connected rate_westmid year, color(maroon%30) msymbol(circle) lcolor(maroon) lstyle(solid) || connected rate_yorkshire year, color(teal%30) msymbol(circle) lcolor(teal) lstyle(solid) ytitle("`ytitleprev'", size(medsmall)) ylabel(`ylab_region', `format_region' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel(`intervention' "COVID-19", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline(`intervention') legend(region(fcolor(white%0)) title("Region", size(small) margin(b=1)) order(1 "East" 2 "East Midlands" 3 "London" 4 "North East" 5 "North West" 6 "South East" 7 "South West" 8 "West Midlands" 9 "Yorkshire Humber")) xsize(16) ysize(9) name(prev_`disease'_region, replace) saving("$projectdir/output/figures/prev_`disease'_region.gph", replace)
			*graph export "$projectdir/output/figures/prev_`disease'_region.png", replace
			graph export "$projectdir/output/figures/prev_`disease'_region.svg", replace

	restore			
}

log close	
