version 16

/*==============================================================================
DO FILE NAME:			Produce temporal plots
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Rounded and redacted data tables
USER-INSTALLED ADO: 	 
  (place .ado file(s) in analysis folder)						
==============================================================================*/

*Set filepaths
/*
global projectdir "C:\Users\k1754142\OneDrive\PhD Project\OpenSAFELY NICE\nice_gout"
global running_locally = 1 // Running on local machine
global img png
*/

global projectdir `c(pwd)'
global running_locally = 0 // Running on OpenSAFELY console
global img svg

capture mkdir "$projectdir/output/data"
capture mkdir "$projectdir/output/figures"
capture mkdir "$projectdir/output/tables"

*Open log file
global logdir "$projectdir/logs"
cap log close
log using "$logdir/temporal_plots.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist disease demographic studystart_date studyend_date studyfup_date intervention_date_2
args $arglist

if $running_locally ==0 {
	foreach var of global arglist {
		local `var' : subinstr local `var' "|" " ", all
		global `var' "``var''"
		di "$`var'"
	}
}

if $running_locally ==1 {
	global disease "gout"
	global demographic "agegroup sex ethnicity imd region"
	global studystart_date "2016-07-01"
	global studyend_date "2025-06-30"
	global studyfup_date "2025-12-31"
	global intervention_date_2 "2022-06-01"
}

di "$disease"
di "$studystart_date"
di "$studyend_date"
di "$studyfup_date"
di "$intervention_date_2"

global intervention_date "$intervention_date_2"
di "$intervention_date"

*Start year, end year, intervention date (derived from above)
foreach date in studystart studyend studyfup intervention {
	local year = real(substr("$`date'_date", 1, 4))
	local month = real(substr("$`date'_date", 6, 2))
	local year_month = ym(`year', `month')
	global `date' `year_month'
	global `date'_year `year'
	*global `date'_tm = `year_month'	
	di %tm $`date'
	di %ty $`date'_year
}

set type double

set scheme plotplainblind

*Single line figures for full cohort ==================================*/

**Loop through data tables with different inclusion criteria

**Baseline data (no additional inclusion criteria)
**Individuals with X months minimum duration of follow-up post-diagnosis
**Individuals with X months minimum duration of follow-up post-ULT initiation, assuming ULT initiation within X months of diagnosis
**Individuals who should have been offered ULT at diagnosis or subsequently on the basis of risk factors
**Individuals prescribed febuxostat (with/without MACE)

foreach table in flare_blood febux_mace ultrisk postult postdiagnosis baseline {

	**Import rounded and redacted data tables
	import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
	di "`table'"

	**Extract outcomes of interest from data table
	levelsof outcome_name, local(outcomes)
	di `outcomes'

	**Loop through outcomes of interest for full cohort
	foreach outcome in `outcomes' {

		preserve
		
		**Keep outcome of interest
		keep if outcome_name == "`outcome'"
		di "`outcome'"
		
		**Reshape to long format
		reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
		gen demog_group = substr(demographic, 1, 3)
		gen demog_level = substr(demographic, 5, .)
		replace demog_level = subinstr(demog_level, "_", " ", .)
		rename count_ count
		rename total_ total
		rename prop_ prop
		replace prop = prop*100 //change to %
		drop demographic
		
		**Keep full cohort only
		keep if demog_group == "all"
		
		**Convert date format
		rename month_year month_year_s
		gen month_year = monthly(month_year_s, "MY") 
		format month_year %tmMon-CCYY
		drop month_year_s
		order month_year, after(outcome_desc)
		
		**Generate 3-monthly moving averages for proportions
		sort month_year
		gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

		**Set y-axis format and title
		local format "format(%9.0f)"
		local ytitle "Percentage of patients"

		***Set x-axis format and title
		local xlabel ""
		forvalues y = $studystart_year(1)`= $studyend_year + 1' {
			local m = ym(`y', 1)
			local xlabel `xlabel' `m' "`y'" 
		}
		di as txt `"`xlabel'"'
		local xtitle ""
		
		***Set title
		local outcome_desc = outcome_desc

		***Temporal plot over study period (scatter with moving average)
		twoway scatter prop month_year, ytitle("`ytitle'", size(medsmall)) color(emerald%30) msymbol(circle) || line prop_ma month_year, lcolor(emerald) lstyle(solid) ylabel(, `format' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel($intervention "NICE Guideline", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline($intervention) legend(off) xsize(16) ysize(9) name("`outcome'", replace) saving("$projectdir/output/figures/temporal_plot_`outcome'.gph", replace)
		graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'.$img", replace

		**ITSA graphs
				
		***Extract minimum follow-up duration from variable name
		if regexm("`outcome'", "([0-9]+)m") {
			local fup_months = regexs(1)
		}
		else {
			local fup_months = 0
		}
		di "`fup_months'"
		
		***Define time window			
		local start = $studystart
		local end = ($studyfup - `fup_months')
		keep if inrange(month_year, `start', `end')
		sort month_year
		
		***Skip ITSA if there are any gaps in monthly series or missing data
		quietly count if missing(prop)
		local n_missing = r(N)
		di `n_missing'
		
		gen d = month_year - month_year[_n-1] if _n>1
		quietly count if d > 1 & !missing(d)
		local n_gaps = r(N)
		drop d
		
		if (`n_missing' > 0 | `n_gaps' > 0) {
			di "Skipping ITSA for `outcome': missing proportions or gaps in time series"
		}
		else {
			
			***Do ITSA if no gaps in time series data (Newey Standard Errors with 5 lags)
			sort month_year
			
			****Set time series
			tsset month_year
			
			****Generate time variable from month 0 and calculates months from intervention
			gen t = month_year - `start'
			local t0 = $intervention - `start'

			****Define post-intervention binary indicator and months from intervention
			gen post   = month_year >= $intervention
			gen t_post = (t - `t0')*post

			****Fit segmented regression with NW SEs (with 5 lags) and 
			newey prop c.t i.post c.t_post, lag(5)
			
			****Extract coefficients (annualised) and p-values
			local b_pre = 12*(_b[t])
			local b_step = 12*(_b[1.post])
			local b_chg = 12*(_b[c.t_post])
			local p_pre: display %6.3f 2*ttail(e(df_r), abs(_b[t]      / _se[t]))
			local p_step: display %6.3f 2*ttail(e(df_r), abs(_b[1.post] / _se[1.post]))
			local p_chg: display %6.3f 2*ttail(e(df_r), abs(_b[c.t_post] / _se[c.t_post]))
			
			****Post-intervention trends
			lincom _b[t] + _b[c.t_post]
			local b_post = 12*(r(estimate))
			local p_post: display %6.3f r(p)
			
			****Formatting
			local f_pre: display %9.2f `b_pre'
			local f_step: display %9.2f `b_step'
			local f_chg: display %9.2f `b_chg'
			local f_post: display %9.2f `b_post'
			
			****Generate text box for key ITSA values
			local boxlines ""

			foreach s in ///
				`"Trend before:"' ///
				`"`f_pre'%/yr (p=`p_pre')"' ///
				`"Trend after:"' ///
				`"`f_post'%/yr (p=`p_post')"' ///
				`"Trend change:"' ///
				`"`f_chg'%/yr (p=`p_chg')"' ///
				`"Step change:"' ///
				`"`f_step'% (p=`p_step')"' {

				local boxlines `"`boxlines' `"`s'"'"'
			}
				
			*Where to place box
			quietly summarize month_year if e(sample), meanonly
			local xmax = r(max)
			local xmin = r(min)
			local xrange = `xmax' - `xmin'
			local xbox = `xmax' + 0.05*`xrange'

			quietly summarize prop if e(sample), meanonly
			local ytop = r(max)
			local ybot = r(min)
			local yrange = `ytop' - `ybot'
			local ybox = `ybot' + 0.15*`yrange'

			****Predicted values for plotting
			predict yhat if e(sample)

			****Plot observed and fitted lines
			twoway scatter prop month_year if e(sample), ytitle("`ytitle'", size(medsmall)) color(emerald%30) msymbol(circle) || line yhat month_year if e(sample) & month_year<$intervention, lcolor(emerald) lstyle(solid) || line yhat month_year if e(sample) & month_year>=$intervention, lcolor(emerald) lstyle(solid) ylabel(, `format' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel($intervention "NICE Guideline", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline($intervention) legend(off) xsize(16) ysize(9) ///
			graphregion(margin(r=20)) plotregion(margin(r=15)) ///
			text(`ybox' `xbox' `boxlines', place(e) just(left) box bcolor(white) margin(small) size(small)) ///
			name("`outcome'", replace)
			graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'_itsa.$img", replace
			
		/*  itsa prop if inrange(month_year, `start', `end'), single trperiod($intervention_tm) lag(5) replace posttrend  ///
			figure(xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(1)) xlabel(`xlabel', nogrid labsize(small) axis(1)) xlabel($intervention "NICE Guideline", axis(2) labsize(small) labcolor(navy)) xtitle("", axis(2)) xscale(noline axis(2)) title("", size(small)) lstyle(dashed) subtitle("", size(medsmall)) ytitle("`ytitle'", size(medsmall) margin(small)) ylabel(, nogrid `format' labsize(small)) note("", size(v.small)) legend(off) name("itsa_`outcome'", replace) ///
			) 
			graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'_itsa.$img", replace
		*/
			actest, lag(18)		
		}
		restore
	}
}

*Multi-line figures by demographic characteristics ==================================*/

**Loop through data tables
foreach table in flare_blood febux_mace ultrisk postult postdiagnosis baseline {
	
	***Import rounded and redacted data tables
	import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
	di "`table'"
	
	tempfile table_base
	save `table_base', replace
		
	**Extract outcomes of interest from data table
	levelsof outcome_name, local(outcomes)
	di `outcomes'

	**Loop through outcomes of interest
	foreach outcome in `outcomes' {
		
		use `table_base', clear
	
		keep if outcome_name == "`outcome'"
		di "`outcome'"
		
		**Reshape to long format
		reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
		gen demog_group = substr(demographic, 1, 3)
		gen demog_level = substr(demographic, 5, .)
		replace demog_level = subinstr(demog_level, "_", " ", .)
		replace demog_level = proper(demog_level)
		rename count_ count
		rename total_ total
		rename prop_ prop
		replace prop = prop*100 //change to %
		drop demographic
		
		**Remove not known categories
		drop if regexm(demog_level, "Not Known")
		
		**Save temporary file for that outcome
		tempfile outcome_base
		save `outcome_base', replace
		
		**Loop through demographic variables of interest
		foreach demog_var in $demographic {
			
			use `outcome_base', clear
			
			di "`demog_var'"
			keep if demog_group == substr("`demog_var'", 1, 3)
			
			***Skip if no observations
			count
			if r(N)==0 continue
	
			***Convert date format
			rename month_year month_year_s
			gen month_year = monthly(month_year_s, "MY") 
			format month_year %tmMon-CCYY
			drop month_year_s
			order month_year, after(outcome_desc)
				
			***Generate 3-monthly moving averages for proportions
			bys demog_level (month_year): gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

			***Set y-axis format and title			
			local format "format(%9.0f)"
			local ytitle "Percentage of patients"

			***Set x-axis format and title
			local xlabel ""
			forvalues y = $studystart_year(1)`= $studyend_year + 1' {
				local m = ym(`y', 1)
				local xlabel `xlabel' `m' "`y'"
			}
			di as txt `"`xlabel'"'
			local xtitle ""
			
			***Set title
			local outcome_desc = outcome_desc
			
			***Choose colour palette based on demographic variable (change as needed)
			if "`demog_var'" == "sex" {
				local colours "red midblue"
				local legtitle ""
			}
			else if "`demog_var'" == "agegroup" {
				local colours "ltblue eltblue midblue ebblue blue navy black"
				local legtitle "Age group"
			}
			else if inlist("`demog_var'", "imd") {
				local colours "ltblue eltblue ebblue blue navy"
				local legtitle "IMD quintile"
			}
			else if inlist("`demog_var'", "ethnicity") {
				local colours "ltblue eltblue ebblue blue navy"
				local legtitle "Ethnicity"
			}
			else if inlist("`demog_var'", "region") {
				local colours "emerald orange red blue dkgreen cranberry navy maroon"
				local legtitle "Region"
			}
			else {
				local colours "emerald orange red blue dkgreen cranberry navy maroon teal sienna purple"
				local legtitle ""
			}

			***Store plots and legend labels
			local plots ""
			local legorder ""
			local leglabels ""
			
			***Extract variables of interest from data table
			levelsof demog_level, local(demog_subset)
			
			local i = 0
			foreach subset of local demog_subset {
				di as txt `"`subset'"'
				local ++i
				local colour : word `i' of `colours'
				if "`colour'"=="" local colour "black"

				/*
				****Two plots per subset: scatter and moving average line
				local thisplot scatter prop month_year if demog_level==`"`subset'"', mcolor(`colour'%20) msymbol(circle) || line prop_ma month_year if demog_level==`"`subset'"', lcolor(`colour') lpattern(solid)
				
				if `i' == 1 local plots `thisplot'
				else local plots `plots' || `thisplot'
				di as txt `"`plots'"'

				****Legend: keep only the moving average lines (2,4,6,...)
				local lineidx = 2*`i'
				local legorder `legorder' `lineidx'
				local outcome_disp : subinstr local subset "_" " " , all   //
				local leglabels `leglabels' label(`lineidx' "`outcome_disp'")
				di as txt `"`leglabels'"'
				*/
			
				****Single plot per subset: moving average line
				local thisplot line prop_ma month_year if demog_level==`"`subset'"', lcolor(`colour') lpattern(solid)
				
				if `i' == 1 local plots `thisplot'
				else local plots `plots' || `thisplot'
				di as txt `"`plots'"'

				****Legend: keep only the moving average lines
				local lineidx = `i'
				local legorder `legorder' `lineidx'
				local outcome_disp : subinstr local subset "_" " " , all
				if demog_group == "reg" {
					local outcome_disp = proper("`outcome_disp'")
				}
				else {
				local outcome_disp = lower("`outcome_disp'")
				local outcome_disp = upper(substr("`outcome_disp'",1,1)) + substr("`outcome_disp'",2,.)
				}
				local leglabels `leglabels' label(`lineidx' "`outcome_disp'")
				di as txt `"`leglabels'"'
			}
			
			***Shorten name if long variable
			local gname = strtoname("`table'_`outcome'_`demog_var'")
			local gname = substr("`gname'", 1, 32)

			***Build plots
			twoway `plots' ytitle("`ytitle'", size(medsmall)) ylabel(, `format' nogrid labsize(small)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel($intervention "NICE Guideline", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline($intervention) legend(region(fcolor(white%0)) title("`legtitle'", size(small) margin(b=1)) order(`legorder') `leglabels') xsize(16) ysize(9) name(`gname', replace) saving("$projectdir/output/figures/temporal_plot_`table'_`outcome'_`demog_var'.gph", replace)
			graph export "$projectdir/output/figures/temporal_plot_`table'_`outcome'_`demog_var'.$img", replace
		}	
	}	
}

*Multi-line figures for selected blood tests and comorbidities (grouped binary variables) (full cohort) =============*/

**Outcome groups of interest (change if needed)
local bloods `" "creatinine_within_12m", "hba1c_within_12m", "cholesterol_within_12m" "'
local comorbidities `" "depression_12m", "cva_12m", "chd_12m", "diabetes_12m", "ckd_comb_12m" "'

foreach table in postdiagnosis {
	
	di "`table'"

	**Loop through outcome groups of interest
	foreach group in bloods comorbidities {
		
		di "`group'"
		
		**Import rounded and redacted data tables
		import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
		
		**Keep outcomes of interest
		keep if inlist(outcome_name, ``group'')
		
		**Reshape to long format
		reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
		gen demog_group = substr(demographic, 1, 3)
		gen demog_level = substr(demographic, 5, .)
		replace demog_level = subinstr(demog_level, "_", " ", .)
		rename count_ count
		rename total_ total
		rename prop_ prop
		replace prop = prop*100 //change to %
		drop demographic
		
		**Keep full cohort only
		keep if demog_group == "all"
		
		**Convert date format
		rename month_year month_year_s
		gen month_year = monthly(month_year_s, "MY") 
		format month_year %tmMon-CCYY
		drop month_year_s
		order month_year, after(outcome_desc)
		
		**Generate 3-monthly moving averages for proportions
		sort month_year
		bys outcome_name (month_year): gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

		**Set y-axis format and title
		local format "format(%9.0f)"
		local ytitle "Percentage of patients"
					
		***Set x-axis format and title
		local xlabel ""
		forvalues y = $studystart_year(1)`= $studyend_year + 1' {
			local m = ym(`y', 1)
			local xlabel `xlabel' `m' "`y'"
		}
		local xtitle ""
		
		***Extract variables of interest from data table
		levelsof outcome_name, local(outcomes)
		di as txt `"`outcomes'"'
		
		***Colour palette to cycle through
		local colours "emerald orange blue dkgreen cranberry navy maroon teal sienna purple"

		***Store plots and legend labels
		local plots ""
		local legorder ""
		local leglabels ""

		local i = 0
		foreach outcome of local outcomes {
			local ++i
			local colour : word `i' of `colours'
			if "`colour'"=="" local colour "black"
			
			****Single plot per subset: moving average line
			local thisplot line prop_ma month_year if outcome_name==`"`outcome'"', lcolor(`colour') lpattern(solid)
			
			if `i' == 1 local plots `thisplot'
			else local plots `plots' || `thisplot'
			di as txt `"`plots'"'
			
			****Get matching outcome description
			levelsof outcome_desc if outcome_name=="`outcome'", local(outcome_disp) clean

			****Legend: keep only the moving average lines
			local lineidx = `i'
			local legorder `legorder' `lineidx'
			local leglabels `leglabels' label(`lineidx' `"`outcome_disp'"')
			di as txt `"`leglabels'"'
		}

		***Build plots
		twoway `plots' ytitle("`ytitle'", size(medsmall)) ylabel(0(20)100, `format' nogrid labsize(small)) yscale(range(0(20)100)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel($intervention "NICE Guideline", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline($intervention) legend(order(`legorder') `leglabels') xsize(16) ysize(9) name("`table'_`group'", replace) saving("$projectdir/output/figures/temporal_plot_`table'_`group'.gph", replace)
		graph export "$projectdir/output/figures/temporal_plot_`table'_`group'.$img", replace
	}
}

**Outcome groups of interest (change if needed)
local drugs `" "allopurinol_ongoing_12m", "febuxostat_ongoing_12m" "'
local legtitle "Drug survival"

foreach table in postult {
	
	di "`table'"

	**Loop through outcome groups of interest
	foreach group in drugs {
		
		di "`group'"
		
		**Import rounded and redacted data tables
		import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
		
		**Keep outcomes of interest
		keep if inlist(outcome_name, ``group'')
		
		**Reshape to long format
		reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
		gen demog_group = substr(demographic, 1, 3)
		gen demog_level = substr(demographic, 5, .)
		replace demog_level = subinstr(demog_level, "_", " ", .)
		rename count_ count
		rename total_ total
		rename prop_ prop
		replace prop = prop*100 //change to %
		drop demographic
		
		**Keep full cohort only
		keep if demog_group == "all"
		
		**Convert date format
		rename month_year month_year_s
		gen month_year = monthly(month_year_s, "MY") 
		format month_year %tmMon-CCYY
		drop month_year_s
		order month_year, after(outcome_desc)
		
		**Generate 3-monthly moving averages for proportions
		sort month_year
		bys outcome_name (month_year): gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

		**Set y-axis format and title
		local format "format(%9.0f)"
		local ytitle "Percentage of patients"
					
		***Set x-axis format and title
		local xlabel ""
		forvalues y = $studystart_year(1)`= $studyend_year + 1' {
			local m = ym(`y', 1)
			local xlabel `xlabel' `m' "`y'"
		}
		local xtitle ""
		
		***Extract variables of interest from data table
		levelsof outcome_name, local(outcomes)
		di as txt `"`outcomes'"'
		
		***Colour palette to cycle through
		local colours "emerald orange blue dkgreen cranberry navy maroon teal sienna purple"

		***Store plots and legend labels
		local plots ""
		local legorder ""
		local leglabels ""

		local i = 0
		foreach outcome of local outcomes {
			local ++i
			local colour : word `i' of `colours'
			if "`colour'"=="" local colour "black"
			
			****Single plot per subset: moving average line
			local thisplot line prop_ma month_year if outcome_name==`"`outcome'"', lcolor(`colour') lpattern(solid)
			
			if `i' == 1 local plots `thisplot'
			else local plots `plots' || `thisplot'
			di as txt `"`plots'"'
			
			****Use shortened variable name for legend
			local pos = strpos("`outcome'", "_")
			if `pos' > 0 {
				local outcome_disp = substr("`outcome'", 1, `pos' - 1)
			}
			else {
				local outcome_disp = "`outcome'"
			}
			local outcome_disp = upper(substr("`outcome_disp'",1,1)) + substr("`outcome_disp'",2,.)

			****Legend: keep only the moving average lines
			local lineidx = `i'
			local legorder `legorder' `lineidx'
			local leglabels `leglabels' label(`lineidx' `"`outcome_disp'"')
			di as txt `"`leglabels'"'
		}

		***Build plots
		twoway `plots' ytitle("`ytitle'", size(medsmall)) ylabel(0(20)100, `format' nogrid labsize(small)) yscale(range(0(20)100)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel($intervention "NICE Guideline", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline($intervention) legend(order(`legorder') `leglabels' title("`legtitle'", size(small) margin(b=1))) xsize(16) ysize(9) name("`table'_`group'", replace) saving("$projectdir/output/figures/temporal_plot_`table'_`group'.gph", replace)
		graph export "$projectdir/output/figures/temporal_plot_`table'_`group'.$img", replace
	}
}

*Multi-line figures for categorical variables (full cohort) ==================================

foreach table in ult_drug flare_drug {
	
	***Import rounded and redacted data tables
	import delimited "$projectdir/output/tables/data_table_`table'.csv", clear
	di "`table'"

	***Convert date format
	rename month_year month_year_s
	gen month_year = monthly(month_year_s, "MY") 
	format month_year %tmMon-CCYY
	drop month_year_s
	order month_year, after(outcome_desc)
	
	**Reshape to long format
	reshape long count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demographic) string
	gen demog_group = substr(demographic, 1, 3)
	gen demog_level = substr(demographic, 5, .)
	replace demog_level = subinstr(demog_level, "_", " ", .)
	rename count_ count
	rename total_ total
	rename prop_ prop
	replace prop = prop*100 //change to %
	drop demographic
	
	**Keep full cohort only
	keep if demog_group == "all"
		
	***Generate 3-monthly moving averages for proportions
	bys outcome_name (month_year): gen prop_ma = (prop[_n-1]+prop[_n]+prop[_n+1])/3

	***Set y-axis format and title
	local format "format(%9.0f)"
	local ytitle "Percentage of patients"

	***Set x-axis format and title
	local xlabel ""
	forvalues y = $studystart_year(1)`= $studyend_year + 1' {
		local m = ym(`y', 1)
		local xlabel `xlabel' `m' "`y'"
	}
	local xtitle ""
	
	***Extract variables of interest from data table
	levelsof outcome_name, local(outcomes)
	di as txt `"`outcomes'"'
	
	***Colour palette to cycle through
    local colours "emerald orange blue dkgreen cranberry navy maroon teal sienna purple"

    ***Store plots and legend labels
    local plots ""
    local legorder ""
    local leglabels ""

    local i = 0
    foreach outcome of local outcomes {
        local ++i
        local colour : word `i' of `colours'
        if "`colour'"=="" local colour "black"
		
		****Single plot per subset: moving average line
		local thisplot line prop_ma month_year if outcome_name==`"`outcome'"', lcolor(`colour') lpattern(solid)
		
		if `i' == 1 local plots `thisplot'
		else local plots `plots' || `thisplot'
		di as txt `"`plots'"'

		****Get matching outcome description
		levelsof outcome_desc if outcome_name=="`outcome'", local(outcome_disp) clean

		****Legend: keep only the moving average lines
		local lineidx = `i'
		local legorder `legorder' `lineidx'
		local leglabels `leglabels' label(`lineidx' `"`outcome_disp'"')
		di as txt `"`leglabels'"'
		}

	***Build plots
    twoway `plots' ytitle("`ytitle'", size(medsmall)) ylabel(0(20)100, `format' nogrid labsize(small)) yscale(range(0(20)100)) xaxis(1 2) xtitle("`xtitle'", size(medsmall) margin(medsmall) axis(2)) xlabel(`xlabel', nogrid labsize(small) axis(2)) xlabel($intervention "NICE Guideline", axis(1) labsize(small) labcolor(navy)) xtitle("", axis(1)) xscale(noline axis(1)) title("", size(medium) margin(b=2)) xline($intervention) legend(order(`legorder') `leglabels') xsize(16) ysize(9) name("`table'", replace) saving("$projectdir/output/figures/temporal_plot_`table'.gph", replace)
	graph export "$projectdir/output/figures/temporal_plot_`table'.$img", replace
}

log close