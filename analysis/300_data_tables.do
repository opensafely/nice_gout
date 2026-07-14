version 16

/*==============================================================================
DO FILE NAME:			Generate rounded and redacted data tables using cleaned primary cohort
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Cleaned primary cohort
USER-INSTALLED ADO: 	 
  (place .ado file(s) in analysis folder)						
==============================================================================*/

*Set filepaths
/*
global projectdir "C:\Users\k1754142\OneDrive\PhD Project\OpenSAFELY NICE\nice_gout"
global running_locally = 1 // Running on local machine
*/

global projectdir `c(pwd)'
global running_locally = 0 // Running on OpenSAFELY console

capture mkdir "$projectdir/output/data"
capture mkdir "$projectdir/output/figures"
capture mkdir "$projectdir/output/tables"

*Open log file
global logdir "$projectdir/logs"
cap log close
log using "$logdir/data_tables.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist disease demographic outpatients
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
	global outpatients "rheumatology"
}
di "$disease"
di "$demographic"
di "$outpatients"

set type double

set scheme plotplainblind

*Function to generate a rounded and redacted data table for binary variables of interest, collapsed by month, for full cohort ======================*/
program define rounded_datatable
    syntax varlist(min=1 max=1), timevar(name)
	local var `varlist'
	local time_variable `timevar'
		
	preserve	
		**Store variable label
		local v : variable label `var'
		
		**Collapse dataset by variable of interest (have to be binary variables: yes 1 and no 0)
		collapse (sum) count_un=`var' (count) total_un=`var', by(`time_variable')

		**Round and redact counts
		*gen count_all=count_un
		gen count_all=round(count_un, 5)
		replace count_all = . if count_all<=7
		drop count_un
		*gen total_all=total_un
		gen total_all=round(total_un, 5)
		replace total_all = . if count_all==.
		drop total_un
		gen prop_all = count_all/total_all

		**Save variable name and labels
		gen outcome_name = "`var'"
		gen outcome_desc = "`v'"
		order outcome_name, first
		order outcome_desc, after(outcome_name)
		rename `time_variable' month_year

		**Save temporary dataset
		capture append using "$projectdir/output/data/data_table.dta"
		save "$projectdir/output/data/data_table.dta", replace	
	restore
end

*Function to generate a rounded and redacted data table for binary variables of interest, collapsed by month, comparing demographic variables of interest (demographic variable should have not known labelled uniquely) ======================*/
program define rounded_datatable_demog
    syntax varlist(min=1 max=1), timevar(name) demogvar(name)
	local var `varlist'
	local time_variable `timevar'
	local demog_variable `demogvar'
		
	preserve	
		**Store variable label
		local v : variable label `var'
		
		**Collapse dataset by variable of interest (have to be binary variables: yes 1 and no 0)
		collapse (sum) count_un=`var' (count) total_un=`var', by(`time_variable' `demog_variable')
		
		**Round and redact counts
		*gen count_=count_un
		gen count_ = round(count_un, 5)
		replace count_ = . if count_ <=7
		drop count_un
		*gen total_=total_un
		gen total_=round(total_un, 5)
		replace total_ = . if count_==.
		drop total_un
		gen prop_ = count_/total_
	
		**Save variable name and labels
		gen outcome_name = "`var'"
		rename `time_variable' month_year
		decode `demog_variable', gen(demog_level)
		replace demog_level = subinstr(demog_level, " ", "_", .)
		*drop `demog_variable'
		gen demog_group = substr("`demog_variable'", 1, 3)
		gen demog_lab = demog_group + "_" + demog_level
		replace demog_lab = subinstr(demog_lab, " ", "_", .)
		replace demog_lab = substr(demog_lab, 1, 32)
		drop demog_group demog_level `demog_variable'
		
		**Reshape wide
		reshape wide count_ total_ prop_, i(month_year outcome_name) j(demog_lab) string
		
		**Marge with base dataset
		merge 1:1 month_year outcome_name using "$projectdir/output/data/data_table.dta", nogen update replace
		order outcome_name outcome_desc month_year count_all total_all prop_all, first 
		save "$projectdir/output/data/data_table.dta", replace	
	restore
end

*Baseline data table (no additional inclusion criteria) =================================*/

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set monthly time variable
local time_variable "${disease}_moyear"
di "`time_variable'"

**Loop through outcomes of interest for full cohort (have to be binary variables: yes 1 and no 0)
foreach var of varlist urate_bl_360_repeat had_urate_bl {
	rounded_datatable `var', timevar(`time_variable')
}

**Set demographic variables of interest to compare across
foreach demog_var of varlist $demographic {

	local demog_variable "`demog_var'"
	di "`demog_variable'"

	**Loop through outcomes of interest by demography
	foreach var of varlist urate_bl_360_repeat had_urate_bl {
		rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_baseline.csv", replace

*Events occurring after diagnosis (restricted to those with t months follow-up) =================================*/

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Loop through time periods of interest
foreach t in 12 {
	
	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear

	**Set inclusion criteria - limited to those who had at least t months duration of follow-up post-diagnosis
	keep if has_`t'm_fup==1

	**Set monthly time variable
	local time_variable "${disease}_moyear"
	
	**Store list of additional variables of interest
	*local outpatients_refopa_`t'm ${outpatients}_refopa_`t'm
	local outpatients_refopa_`t'm_risk ${outpatients}_refopa_`t'm_risk
	
	**Outpatient appointment outcomes only (OPA data available from July 2019 onwards)
	preserve

		keep if `time_variable' >= tm(2019m7)

		foreach var of varlist `outpatients_refopa_`t'm_risk' {
			rounded_datatable `var', timevar(`time_variable')
		}

		foreach demog_var of varlist $demographic {

			local demog_variable "`demog_var'"
			di "`demog_variable'"
			
			foreach var of varlist `outpatients_refopa_`t'm_risk' {
				rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
			}
		}

	restore

	**Loop through all other outcomes of interest for full cohort (have to be binary variables: yes 1 and no 0)
	foreach var of varlist chd_12m diabetes_12m cva_12m ckd_comb_12m depression_12m creatinine_within_`t'm hba1c_within_`t'm cholesterol_within_`t'm ult_`t'm {
		rounded_datatable `var', timevar(`time_variable')
	}
	
	**Set demographic variables of interest to compare across
	foreach demog_var of varlist $demographic {

		local demog_variable "`demog_var'"
		di "`demog_variable'"

		**Loop through outcomes of interest by demography
		foreach var of varlist chd_12m diabetes_12m cva_12m ckd_comb_12m depression_12m creatinine_within_`t'm hba1c_within_`t'm cholesterol_within_`t'm ult_`t'm {
			rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
		}
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_postdiagnosis.csv", replace

*Events occurring after ULT initiation (restricted to those with t months follow-up after ULT initiation) =================================*/

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Loop through time periods of interest
foreach t in 12 {
	
	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear

	**Set inclusion criteria - limited to those who had at least t months duration of follow-up post-ULT (no restriction on ULT within 12m)
	keep if has_`t'm_fup_ult==1

	**Set monthly time variable (date of first ULT prescription)
	gen ult_first_date_my = mofd(ult_first_date)
	format ult_first_date_my %tmMon-CCYY
	local time_variable "ult_first_date_my"

	**Loop through outcomes of interest for full cohort (have to be binary variables: yes 1 and no 0)
	foreach var of varlist febuxostat_ongoing_`t'm allopurinol_ongoing_`t'm ult_ongoing_`t'm ult_high repeat_below360_`t'm_ult repeat_after360_`t'm_ult ult_prophylaxis_2 ult_prophylaxis urate_`t'm_ult two_urate_`t'm_ult urate_within_`t'm_ult {
		rounded_datatable `var', timevar(`time_variable')
	}
	
	**Set demographic variables of interest to compare across
	foreach demog_var of varlist $demographic {

		local demog_variable "`demog_var'"
		di "`demog_variable'"

		**Loop through outcomes of interest by demography
		foreach var of varlist febuxostat_ongoing_`t'm allopurinol_ongoing_`t'm ult_ongoing_`t'm ult_high repeat_below360_`t'm_ult repeat_after360_`t'm_ult ult_prophylaxis_2 ult_prophylaxis urate_`t'm_ult two_urate_`t'm_ult urate_within_`t'm_ult {
			rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
		}
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
export delimited using "$projectdir/output/tables/data_table_postult.csv", replace

*Choice of first ULT drug (categorical variable, therefore processed differently) ===========================

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Loop through time periods of interest
foreach t in 12 {

	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear

	**Set monthly time variable (date of first ULT prescription)
	gen ult_first_date_my = mofd(ult_first_date)
	format ult_first_date_my %tmMon-CCYY
	local time_variable "ult_first_date_my"

	**Set inclusion criteria (not limited to those receiving ULT within 12m)

	**Generate binary flare treatments variables and label them
	local varlab : val label ult_first_drug
	levelsof ult_first_drug, local(levels)

	local varlist ""

	foreach drug of local levels {
		local lab : lab `varlab' `drug'
		local name = strtoname("`lab'")
		local newname = substr("`name'", 1, 24) //shortens name if necessary
		gen `newname' = (ult_first_drug == `drug')
		lab var `newname' "`lab'"
		
		local varlist `varlist' `newname'
		di `varlist'
	}

	**Loop through outcomes of interest for full cohort
	foreach var of local varlist {
		rounded_datatable `var', timevar(`time_variable')
	}

	**Set demographic variables of interest to compare across
	foreach demog_var of varlist $demographic {

		local demog_variable "`demog_var'"
		di "`demog_variable'"

		**Loop through outcomes of interest by demography
		foreach var of local varlist {
			rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
		}
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
drop if month_year ==. //drop missing
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_ult_drug.csv", replace

*ULT initiation in individuals who should have been offered ULT at diagnosis or subsequently on the basis of risk factors ===========================*/

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Loop through time periods of interest
foreach t in 12 {
	
	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	**Set monthly time variable (date patient became at risk for ULT)
	gen ult_risk_date_my = mofd(ult_risk_date_dx)
	format ult_risk_date_my %tmMon-CCYY
	local time_variable "ult_risk_date_my"
	
	**Set inclusion criteria - should be limited to those who had at least t months duration of post-risk date
	keep if has_`t'm_fup_risk==1

	**Loop through outcomes of interest for full cohort (have to be binary variables: yes 1 and no 0)
	foreach var of varlist ult_risk_p_`t'm {
		rounded_datatable `var', timevar(`time_variable')
	}
	
	**Set demographic variables of interest to compare across
	foreach demog_var of varlist $demographic {

		local demog_variable "`demog_var'"
		di "`demog_variable'"

		**Loop through outcomes of interest by demography
		foreach var of varlist ult_risk_p_`t'm {
			rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
		}
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
drop if month_year ==. //drop missing
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_ultrisk.csv", replace

*Febuxostat prescribed in people with pre-existing MACE (CHD or CVA) =================================*/

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set inclusion criteria - limited to those who had received febuxostat
keep if febuxostat_first_date!=.

**Set monthly time variable (date of first ULT prescription)
gen febux_first_date_my = mofd(febuxostat_first_date)
format febux_first_date_my %tmMon-CCYY
local time_variable "febux_first_date_my"

**Loop through outcomes of interest for full cohort (have to be binary variables: yes 1 and no 0)
foreach var of varlist febux_mace {
	rounded_datatable `var', timevar(`time_variable')
}

**Set demographic variables of interest to compare across
foreach demog_var of varlist $demographic {

	local demog_variable "`demog_var'"
	di "`demog_variable'"

	**Loop through outcomes of interest by demography
	foreach var of varlist febux_mace {
		rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
export delimited using "$projectdir/output/tables/data_table_febux_mace.csv", replace

*Events in individuals having their first flare after diagnosis (binary variables) ===========================*/

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Loop through time periods of interest
foreach t in 12 {
	
	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	**Set monthly time variable (date of first flare)
	gen first_flare_date_my = mofd(first_flare_overall_date)
	format first_flare_date_my %tmMon-CCYY
	local time_variable "first_flare_date_my"
	
	**Set inclusion criteria - limited to those who had at least 12 months of follow-up after first flare date
	keep if has_`t'm_fup_flare==1 //may not need 12 months

	**Loop through outcomes of interest for full cohort (have to be binary variables: yes 1 and no 0)
	foreach var of varlist post_flare_urate {
		rounded_datatable `var', timevar(`time_variable')
	}
	
	**Set demographic variables of interest to compare across
	foreach demog_var of varlist $demographic {

		local demog_variable "`demog_var'"
		di "`demog_variable'"

		**Loop through outcomes of interest by demography
		foreach var of varlist post_flare_urate {
			rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
		}
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
drop if month_year ==. //drop those not at risk 
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_flare_blood.csv", replace

*Treatment of first flare after diagnosis (categorical variable, therefore processed differently) ===========================

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Loop through time periods of interest
foreach t in 12 {

	**Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear

	**Set monthly time variable (date of first flare)
	gen first_flare_date_my = mofd(first_flare_overall_date)
	format first_flare_date_my %tmMon-CCYY
	local time_variable "first_flare_date_my"

	**Set inclusion criteria - limited to those who had at least 12 months of follow-up after first flare date
	keep if has_`t'm_fup_flare==1

	**Generate binary flare treatments variables and label them
	local varlab : val label first_flare_drug
	levelsof first_flare_drug, local(levels)

	local varlist ""

	foreach drug of local levels {
		local lab : lab `varlab' `drug'
		local name = strtoname("`lab'")
		local newname = substr("`name'", 1, 24) //shortens name if necessary
		gen `newname' = (first_flare_drug == `drug')
		lab var `newname' "`lab'"
		
		local varlist `varlist' `newname'
		di `varlist'
	}

	**Loop through outcomes of interest for full cohort
	foreach var of local varlist {
		rounded_datatable `var', timevar(`time_variable')
	}

	**Set demographic variables of interest to compare across
	foreach demog_var of varlist $demographic {

		local demog_variable "`demog_var'"
		di "`demog_variable'"

		**Loop through outcomes of interest by demography
		foreach var of local varlist {
			rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
		}
	}
}

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
drop if month_year ==. //drop missing
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_flare_drug.csv", replace

/*
*Treatment of any flare after diagnosis (categorical variables) - Nb. would need to have full list of drug use after diagnosis ===========================

**Erase any existing data file
capture erase "$projectdir/output/data/data_table.dta"

**Import cleaned/processed cohort (long format flare data)
use "$projectdir/output/data/flares_long.dta", clear

**Generate binary flare treatments variables and label them
local varlab : val label flare_drug_
levelsof flare_drug_, local(levels)

local varlist ""

foreach drug of local levels {
    local lab : lab `varlab' `drug'
    local name = strtoname("`lab'")
    local newname = substr("`name'", 1, 24) //shortens name if necessary
    gen `newname' = (flare_drug_ == `drug')
    lab var `newname' "`lab'"
	
	local varlist `varlist' `newname'
	di `varlist'
}

**Set monthly time variable
gen flare_date_my = mofd(flare_overall_date_)
format flare_date_my %tmMon-CCYY
local time_variable "flare_date_my"

**Loop through outcomes of interest for full cohort
foreach var of local varlist {
	rounded_datatable `var', timevar(`time_variable')
}

/*
**Set demographic variables of interest to compare across
foreach demog_var of varlist $demographic {

	local demog_variable "`demog_var'"
	di "`demog_variable'"

	**Loop through outcomes of interest by demography
	foreach var of local varlist {
		rounded_datatable_demog `var', timevar(`time_variable') demogvar(`demog_variable')
	}
}
*/

**Export rounded/redacted data table
use "$projectdir/output/data/data_table.dta", clear
drop if month_year ==. //drop missing
sort outcome_name month_year
export delimited using "$projectdir/output/tables/data_table_flares.csv", replace
*/

log close
