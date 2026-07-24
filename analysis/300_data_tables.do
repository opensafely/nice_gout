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
log off

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

program define rounded_datatable_multi
    syntax varlist(min=1), timevar(name) outfile(string)

    preserve
        local time_variable `timevar'

        **Create unrounded numerator and denominator variables
        local collapse_list ""
		local n_outcomes 0

		foreach var of local varlist {
			local ++n_outcomes

			clonevar c`n_outcomes' = `var'
			clonevar n`n_outcomes' = `var'

			local collapse_list `collapse_list' (sum) c`n_outcomes'=c`n_outcomes' (count) n`n_outcomes'=n`n_outcomes'

			**Store original variable name and label
			local outcome_name`n_outcomes' "`var'"
			local outcome_desc`n_outcomes' : variable label `var'
		}

        **One collapse for every outcome
        collapse `collapse_list', by(`time_variable')
        rename `time_variable' month_year

        **Reshape so each outcome becomes a row
        reshape long c n, i(month_year) j(outcome_number)
		rename c count_un
		rename n total_un

        **Round and redact
        gen count_all = round(count_un, 5)
        gen total_all = round(total_un, 5)
        replace count_all = . if count_un <= 7
        replace total_all = . if missing(count_all)
        gen prop_all = count_all / total_all
        drop count_un total_un

        **Recover variable labels
		gen str32 outcome_name = ""
		gen str244 outcome_desc = ""

		forvalues i = 1/`n_outcomes' {
			replace outcome_name = "`outcome_name`i''" ///
				if outcome_number == `i'

			replace outcome_desc = `"`outcome_desc`i''"' ///
				if outcome_number == `i'
		}

		drop outcome_number

        order outcome_name outcome_desc month_year count_all total_all prop_all
        save `"`outfile'"', replace
    restore
end

*Function to generate a rounded and redacted data table for binary variables of interest, collapsed by month, comparing demographic variables of interest (demographic variable should have not known labelled uniquely) ======================*/

program define rounded_datatable_demog_multi
    syntax varlist(min=1), timevar(name) demogvar(name) outfile(string)

    preserve
        local time_variable `timevar'
        local demog_variable `demogvar'

        local collapse_list ""
		local n_outcomes 0

		foreach var of local varlist {
			local ++n_outcomes

			clonevar c`n_outcomes' = `var'
			clonevar n`n_outcomes' = `var'

			local collapse_list `collapse_list' (sum) c`n_outcomes'=c`n_outcomes' (count) n`n_outcomes'=n`n_outcomes'

			local outcome_name`n_outcomes' "`var'"
			local outcome_desc`n_outcomes' : variable label `var'
		}
		
        **One collapse for all outcomes
        collapse `collapse_list', by(`time_variable' `demog_variable')

        rename `time_variable' month_year
        decode `demog_variable', gen(demog_level)
        replace demog_level = subinstr(demog_level, " ", "_", .)
        local demog_group = substr("`demog_variable'", 1, 3)
        gen demog_lab = "`demog_group'_" + demog_level
        replace demog_lab = substr(strtoname(demog_lab), 1, 32)
        drop demog_level `demog_variable'

        **Move outcomes into rows
        reshape long c n, i(month_year demog_lab) j(outcome_number)
		rename c count_un
		rename n total_un
		
        gen count_ = round(count_un, 5)
        gen total_ = round(total_un, 5)
        replace count_ = . if count_un <= 7
        replace total_ = . if missing(count_)
        gen prop_ = count_ / total_
        drop count_un total_un

		gen str32 outcome_name = ""
		gen str244 outcome_desc = ""

		forvalues i = 1/`n_outcomes' {
			replace outcome_name = "`outcome_name`i''" ///
				if outcome_number == `i'

			replace outcome_desc = `"`outcome_desc`i''"' ///
				if outcome_number == `i'
		}

		drop outcome_number

        reshape wide count_ total_ prop_, i(month_year outcome_name outcome_desc) j(demog_lab) string

        save `"`outfile'"', replace
    restore
end

*Baseline data table (no additional inclusion criteria) =================================*/

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set monthly time variable
local time_variable "${disease}_moyear"

**Define all outcomes
local main_outcomes urate_bl_360_repeat had_urate_bl

**Define outcomes requiring demographic breakdowns
local demog_outcomes had_urate_bl

**Generate full-cohort results for all outcomes
tempfile full_results

rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

**Generate selected outcomes by demographic group
local n_demog_files 0

if `"`demog_outcomes'"' != "" {

    foreach demog_var of varlist $demographic {

        local ++n_demog_files
        tempfile demog_result`n_demog_files'

        rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

        local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
    }
}

**Combine demographic result files horizontally
if `n_demog_files' > 0 {

    forvalues i = 1/`n_demog_files' {
        if `i' == 1 {
            use `"`demog_file`i''"', clear
        }
        else {
            merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
        }
    }
    save "$projectdir/output/data/combined_demog.dta", replace
}

**Return to full-cohort results
use "`full_results'", clear

**Merge demographic columns onto selected full-cohort outcomes
if `n_demog_files' > 0 {

    merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog.dta"

    **Using-only observations should not occur
    assert _merge != 2
    drop _merge
}

**Drop missing time periods
drop if missing(month_year)

**Check uniqueness and sort
isid month_year outcome_name
sort outcome_name outcome_desc month_year

**Save final tables
save "$projectdir/output/data/data_table_baseline.dta", replace
export delimited using "$projectdir/output/tables/data_table_baseline.csv", replace

*Events occurring after diagnosis (restricted to those with t months follow-up) =================================*/

foreach t in 12 {

    **Import cleaned/processed cohort
    use "$projectdir/output/data/cohort_processed.dta", clear

    **Set inclusion criteria
    keep if has_`t'm_fup == 1

    **Set monthly time variable
    local time_variable "${disease}_moyear"

    **Define all non-outpatient outcomes
    local main_outcomes chd_12m diabetes_12m cva_12m ckd_comb_12m depression_12m creatinine_within_`t'm hba1c_within_`t'm cholesterol_within_`t'm ult_`t'm
	
	 **Define outpatient outcomes (processed separately due to availability of OPA data)
    local opa_outcomes "${outpatients}_refopa_`t'm_risk"

    **Define outcomes requiring demographic breakdowns
    local demog_outcomes chd_12m diabetes_12m cva_12m ckd_comb_12m depression_12m ult_`t'm

    **Generate outpatient results for full cohort only
    preserve
        keep if `time_variable' >= tm(2019m7)
		
        tempfile opa_results

        rounded_datatable_multi `opa_outcomes', timevar(`time_variable') outfile(`opa_results')
    restore

    **Generate full-cohort results for all other outcomes
    tempfile main_results

    rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`main_results')

    **Generate selected outcomes by demographic group
    local n_demog_files 0

    if `"`demog_outcomes'"' != "" {

        foreach demog_var of varlist $demographic {

            local ++n_demog_files
            tempfile demog_result`n_demog_files'

            rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

            local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
        }
    }

    **Combine demographic result files horizontally
    if `n_demog_files' > 0 {
        forvalues i = 1/`n_demog_files' {
            if `i' == 1 {
                use `"`demog_file`i''"', clear
            }
            else {
                merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
            }
        }
        save "$projectdir/output/data/combined_demog_postdiagnosis.dta", replace
    }

    **Combine full-cohort results vertically
    use "`main_results'", clear
    append using "`opa_results'"

    **Merge demographic columns onto selected outcomes
    if `n_demog_files' > 0 {

        merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_postdiagnosis.dta"

        **Full-cohort-only outcomes should be master-only
        assert _merge != 2
        drop _merge
    }

    **Drop missing time periods
    drop if missing(month_year)

    **Check uniqueness and sort
    isid month_year outcome_name
    sort outcome_name outcome_desc month_year

    **Save final tables
    save "$projectdir/output/data/data_table_postdiagnosis.dta", replace
    export delimited using "$projectdir/output/tables/data_table_postdiagnosis.csv", replace
}

*Events occurring at ULT initiation =================================*/

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set inclusion criteria - limited to patients who initiated ULT
keep if ult_first_date != .

**Set monthly time variable (date of first ULT prescription)
gen ult_first_date_my = mofd(ult_first_date)
format ult_first_date_my %tmMon-CCYY
local time_variable "ult_first_date_my"

**Define all outcomes
local main_outcomes ult_high ult_prophylaxis

**Define outcomes requiring demographic breakdowns
local demog_outcomes ult_high ult_prophylaxis

**Generate full-cohort results
tempfile full_results

rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

**Generate selected outcomes by demographic group
local n_demog_files 0

if `"`demog_outcomes'"' != "" {
    foreach demog_var of varlist $demographic {
		
        local ++n_demog_files
        tempfile demog_result`n_demog_files'

        rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

        local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
    }
}

**Combine demographic result files horizontally
if `n_demog_files' > 0 {
    forvalues i = 1/`n_demog_files' {
        if `i' == 1 {
            use `"`demog_file`i''"', clear
        }
        else {
            merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
        }
    }
    save "$projectdir/output/data/combined_demog_atultinitiation.dta", replace
}

**Return to full-cohort results
use "`full_results'", clear

**Merge demographic columns onto selected outcomes
if `n_demog_files' > 0 {

    merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_atultinitiation.dta"

    assert _merge != 2
    drop _merge
}

**Drop missing time periods
drop if missing(month_year)

**Sort and save final tables
isid month_year outcome_name
sort outcome_name outcome_desc month_year

save "$projectdir/output/data/data_table_atultinitiation.dta", replace
export delimited using "$projectdir/output/tables/data_table_atultinitiation.csv", replace

*Events occurring after ULT initiation (restricted to those with t months follow-up after ULT initiation) =================================*/

foreach t in 12 {

    **Import cleaned/processed cohort
    use "$projectdir/output/data/cohort_processed.dta", clear

    **Set inclusion criteria - limited to those who had at least t months duration of follow-up post-ULT
    keep if has_`t'm_fup_ult == 1

    **Set monthly time variable (date of first ULT prescription)
    gen ult_first_date_my = mofd(ult_first_date)
    format ult_first_date_my %tmMon-CCYY
    local time_variable "ult_first_date_my"

    **Define all outcomes
    local main_outcomes febuxostat_ongoing_`t'm allopurinol_ongoing_`t'm ult_ongoing_`t'm urate_`t'm_ult two_urate_`t'm_ult urate_within_`t'm_ult

    **Define outcomes requiring demographic breakdowns
    local demog_outcomes ult_ongoing_`t'm urate_`t'm_ult two_urate_`t'm_ult urate_within_`t'm_ult

    **Generate full-cohort results
    tempfile full_results

    rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

    **Generate selected outcomes by demographic group
    local n_demog_files 0

    if `"`demog_outcomes'"' != "" {

        foreach demog_var of varlist $demographic {
			
            local ++n_demog_files
            tempfile demog_result`n_demog_files'

            rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

            local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
        }
    }

    **Combine demographic result files horizontally
    if `n_demog_files' > 0 {

        forvalues i = 1/`n_demog_files' {
            if `i' == 1 {
                use `"`demog_file`i''"', clear
            }
            else {
                merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
            }
        }
        save "$projectdir/output/data/combined_demog_postult.dta", replace
    }

    **Return to full-cohort results
    use "`full_results'", clear

    **Merge demographic columns onto selected outcomes
    if `n_demog_files' > 0 {

        merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_postult.dta"

        assert _merge != 2
        drop _merge
    }

    **Drop missing time periods
    drop if missing(month_year)

    **Sort and save final tables
    isid month_year outcome_name
    sort outcome_name outcome_desc month_year

    save "$projectdir/output/data/data_table_postult.dta", replace
    export delimited using "$projectdir/output/tables/data_table_postult.csv", replace
}

*Events occurring after urate target attainment (restricted to those with t months follow-up after target attainment) =================================*/

foreach t in 12 {

    **Import cleaned/processed cohort
    use "$projectdir/output/data/cohort_processed.dta", clear

    **Set inclusion criteria - limited to those who had at least t months duration of follow-up post-target attainment
    keep if has_`t'm_fup_target == 1

    **Check - for dummy data only
    count
    if r(N) == 0 {
        di as text "No observations with has_`t'm_fup_target == 1; skipping `t'-month analysis"
        continue
    }

    **Set monthly time variable (date of target attainment)
    gen target_first_date_my = mofd(urate_below360_ult_date)
    format target_first_date_my %tmMon-CCYY
    local time_variable "target_first_date_my"

    **Define all outcomes
    local main_outcomes repeat_below360_`t'm_ult repeat_after360_`t'm_ult

    **Define outcomes requiring demographic breakdowns
    local demog_outcomes

    **Generate full-cohort results
    tempfile full_results

    rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

    **Generate selected outcomes by demographic group
    local n_demog_files 0

    if `"`demog_outcomes'"' != "" {

        foreach demog_var of varlist $demographic {

            local ++n_demog_files
            tempfile demog_result`n_demog_files'

            rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

            local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
        }
    }

    **Combine demographic result files horizontally
    if `n_demog_files' > 0 {

        forvalues i = 1/`n_demog_files' {
            if `i' == 1 {
                use `"`demog_file`i''"', clear
            }
            else {
                merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
            }
        }
        save "$projectdir/output/data/combined_demog_posttarget.dta", replace
    }

    **Return to full-cohort results
    use "`full_results'", clear

    **Merge demographic columns onto selected outcomes
    if `n_demog_files' > 0 {
		
        merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_posttarget.dta"

        assert _merge != 2
        drop _merge
    }

    **Drop missing time periods
    drop if missing(month_year)

    **Sort and save final tables
    isid month_year outcome_name
    sort outcome_name outcome_desc month_year

    save "$projectdir/output/data/data_table_posttarget.dta", replace
    export delimited using "$projectdir/output/tables/data_table_posttarget.csv", replace
}

*ULT initiation in individuals who should have been offered ULT at diagnosis or subsequently on the basis of risk factors ===========================*/

foreach t in 12 {

    **Import cleaned/processed cohort
    use "$projectdir/output/data/cohort_processed.dta", clear

    **Set inclusion criteria - limited to those who had at least t months duration of follow-up post-risk date
    keep if has_`t'm_fup_risk == 1

    **Set monthly time variable (date patient became at risk for ULT)
    gen ult_risk_date_my = mofd(ult_risk_date_dx)
    format ult_risk_date_my %tmMon-CCYY
    local time_variable "ult_risk_date_my"

    **Define all outcomes
    local main_outcomes ult_risk_p_`t'm

    **Define outcomes requiring demographic breakdowns
    local demog_outcomes ult_risk_p_`t'm

    **Generate full-cohort results
    tempfile full_results

    rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

    **Generate selected outcomes by demographic group
    local n_demog_files 0

    if `"`demog_outcomes'"' != "" {

        foreach demog_var of varlist $demographic {

            local ++n_demog_files
            tempfile demog_result`n_demog_files'

            rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

            local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
        }
    }

    **Combine demographic result files horizontally
    if `n_demog_files' > 0 {

        forvalues i = 1/`n_demog_files' {
            if `i' == 1 {
                use `"`demog_file`i''"', clear
            }
            else {
                merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
            }
        }
        save "$projectdir/output/data/combined_demog_ultrisk.dta", replace
    }

    **Return to full-cohort results
    use "`full_results'", clear

    **Merge demographic columns onto full-cohort rows
    if `n_demog_files' > 0 {

        merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_ultrisk.dta"

        assert _merge != 2
        drop _merge
    }

    **Drop missing time periods
    drop if missing(month_year)

    **Sort and save final tables
    isid month_year outcome_name
    sort outcome_name outcome_desc month_year

    save "$projectdir/output/data/data_table_ultrisk.dta", replace
    export delimited using "$projectdir/output/tables/data_table_ultrisk.csv", replace
}

*Febuxostat prescribed in people with pre-existing MACE (CHD or CVA) =================================*/

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set inclusion criteria - limited to those who had received febuxostat
keep if febuxostat_first_date != .

**Set yearly time variable - yearly rather than monthly due to small numbers
gen febux_first_fy = year(febuxostat_first_date) - (month(febuxostat_first_date) < 7)
local time_variable "febux_first_fy"

**Define all outcomes
local main_outcomes febux_mace

**Define outcomes requiring demographic breakdowns
local demog_outcomes

**Generate full-cohort results
tempfile full_results

rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

**Generate selected outcomes by demographic group
local n_demog_files 0

if `"`demog_outcomes'"' != "" {

    foreach demog_var of varlist $demographic {

        local ++n_demog_files
        tempfile demog_result`n_demog_files'

        rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

        local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
    }
}

**Combine demographic result files horizontally
if `n_demog_files' > 0 {

    forvalues i = 1/`n_demog_files' {
        if `i' == 1 {
            use `"`demog_file`i''"', clear
        }
        else {
            merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
        }
    }
    save "$projectdir/output/data/combined_demog_febux_mace.dta", replace
}

**Return to full-cohort results
use "`full_results'", clear

**Merge demographic columns onto selected outcomes
if `n_demog_files' > 0 {

    merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_febux_mace.dta"

    assert _merge != 2
    drop _merge
}

**Drop missing time periods
drop if missing(month_year)

**Sort and save final tables
isid month_year outcome_name
sort outcome_name outcome_desc month_year

save "$projectdir/output/data/data_table_febux_mace.dta", replace
export delimited using "$projectdir/output/tables/data_table_febux_mace.csv", replace

*Events in individuals having their first flare after diagnosis (binary variables) ===========================*/

foreach t in 12 {

    **Import cleaned/processed cohort
    use "$projectdir/output/data/cohort_processed.dta", clear

    **Set inclusion criteria - limited to those who had at least 12 months of follow-up after first flare date
    keep if has_`t'm_fup_flare == 1

    **Set monthly time variable (date of first flare)
    gen first_flare_date_my = mofd(first_flare_overall_date)
    format first_flare_date_my %tmMon-CCYY
    local time_variable "first_flare_date_my"

    **Define all outcomes
    local main_outcomes post_flare_urate

    **Define outcomes requiring demographic breakdowns
    local demog_outcomes

    **Generate full-cohort results
    tempfile full_results

    rounded_datatable_multi `main_outcomes', timevar(`time_variable') outfile(`full_results')

    **Generate selected outcomes by demographic group
    local n_demog_files 0

    if `"`demog_outcomes'"' != "" {

        foreach demog_var of varlist $demographic {

            local ++n_demog_files
            tempfile demog_result`n_demog_files'

            rounded_datatable_demog_multi `demog_outcomes', timevar(`time_variable') demogvar(`demog_var') outfile(`demog_result`n_demog_files'')

            local demog_file`n_demog_files' `"`demog_result`n_demog_files''"'
        }
    }

    **Combine demographic result files horizontally
    if `n_demog_files' > 0 {

        forvalues i = 1/`n_demog_files' {
            if `i' == 1 {
                use `"`demog_file`i''"', clear
            }
            else {
                merge 1:1 month_year outcome_name outcome_desc using `"`demog_file`i''"', nogen
            }
        }
        save "$projectdir/output/data/combined_demog_flare_blood.dta", replace
    }

    **Return to full-cohort results
    use "`full_results'", clear

    **Merge demographic columns onto selected outcomes
    if `n_demog_files' > 0 {

        merge 1:1 month_year outcome_name outcome_desc using "$projectdir/output/data/combined_demog_flare_blood.dta"

        assert _merge != 2
        drop _merge
    }

    **Drop missing time periods
    drop if missing(month_year)

    **Sort and save final tables
    isid month_year outcome_name
    sort outcome_name outcome_desc month_year

    save "$projectdir/output/data/data_table_flare_blood.dta", replace
    export delimited using "$projectdir/output/tables/data_table_flare_blood.csv", replace
}

*Choice of first ULT drug (categorical variable, therefore processed differently) ===========================

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set inclusion criteria - limited to those receiving ULT at any point
keep if ult_first_date !=.

**Set monthly time variable - date of first ULT prescription
gen ult_first_date_my = mofd(ult_first_date)
format ult_first_date_my %tmMon-CCYY
local time_variable "ult_first_date_my"

**Define outcome
local categorical_var ult_first_drug

**Generate a binary indicator for each category of first ULT drug
local varlab : val label `categorical_var'
levelsof `categorical_var', local(levels)

local outcome_vars

foreach subcat of local levels {

	local lab : label `varlab' `subcat'
	local name = strtoname("`lab'")
	local newname = substr("`name'", 1, 24)

	gen `newname' = (`categorical_var' == `subcat') if !missing(`categorical_var')
	label variable `newname' "`lab'"

	local outcome_vars `outcome_vars' `newname'
}

di "`outcome_vars'"

**Generate rounded and redacted results for all drug categories
tempfile results

rounded_datatable_multi `outcome_vars', timevar(`time_variable') outfile(`results')

**Prepare final table
use "`results'", clear

drop if missing(month_year)

**Optional - Remove outcomes with no non-missing proportions
bys outcome_desc: egen has_data = max(!missing(prop_all))
drop if has_data == 0
drop has_data

**Group all binary drug-category indicators under one outcome
replace outcome_name = "`categorical_var'"

isid month_year outcome_name outcome_desc
sort outcome_name outcome_desc month_year

**Save and export final table
save "$projectdir/output/data/data_table_ult_drug.dta", replace
export delimited using "$projectdir/output/tables/data_table_ult_drug.csv", replace

/*Treatment of first flare after diagnosis (categorical variable, therefore processed differently) ===========================

**Import cleaned/processed cohort
use "$projectdir/output/data/cohort_processed.dta", clear

**Set inclusion criteria - limited to those receiving ULT at any point
keep if first_flare_overall_date !=.

**Set monthly time variable (date of first flare)
gen first_flare_date_my = mofd(first_flare_overall_date)
format first_flare_date_my %tmMon-CCYY
local time_variable "first_flare_date_my"

**Define outcome
local categorical_var first_flare_drug

**Generate binary flare treatments variables and label them
local varlab : val label `categorical_var'
levelsof `categorical_var', local(levels)

local outcome_vars

foreach subcat of local levels {

	local lab : label `varlab' `subcat'
	local name = strtoname("`lab'")
	local newname = substr("`name'", 1, 24)

	gen `newname' = (`categorical_var' == `subcat') if !missing(`categorical_var')
	label variable `newname' "`lab'"

	local outcome_vars `outcome_vars' `newname'
}

di "`outcome_vars'"

**Generate rounded and redacted results for all drug categories
tempfile results

rounded_datatable_multi `outcome_vars', timevar(`time_variable') outfile(`results')

**Prepare final table
use "`results'", clear

drop if missing(month_year)

**Group all binary drug-category indicators under one outcome
replace outcome_name = "`categorical_var'"

isid month_year outcome_name outcome_desc
sort outcome_name outcome_desc month_year

**Save and export final table
save "$projectdir/output/data/data_table_flare_drug.dta", replace
export delimited using "$projectdir/output/tables/data_table_flare_drug.csv", replace
*/

*Treatment of any flare after diagnosis (categorical variables) - Nb. would need to have full list of drug use after diagnosis ===========================

**Import cleaned/processed cohort
use "$projectdir/output/data/flares_long.dta", clear

**Set inclusion criteria - limited to those receiving ULT at any point
keep if flare_overall_date_ !=.

**Set monthly time variable
gen flare_date_my = mofd(flare_overall_date_)
format flare_date_my %tmMon-CCYY
local time_variable "flare_date_my"

local outcome_vars ""
local n = 0

log on

**Loop through categorical variables of interest
foreach categorical_var in flare_drug_ flare_source {

	**Generate binary flare treatments variables and label them
	local varlab : val label `categorical_var'
	levelsof `categorical_var', local(levels)
	
	di as text "Levels: `levels'"

	local clean_name = regexr("`categorical_var'", "_$", "")

	foreach subcat of local levels {
		local ++n
		local lab : label `varlab' `subcat'
		di as text "  Level `subcat' = `lab'"
		local newname = "outcome`n'"
		gen `newname' = (`categorical_var' == `subcat') if !missing(`categorical_var')
		label variable `newname' "`lab'"
		
		**Check binary variable
        tab `newname', missing
        tab `categorical_var' `newname', missing
		
		local outcome_vars `outcome_vars' `newname'
		local outcome_name`n' "`clean_name'"
		local outcome_desc`n' "`lab'"
	}
}
log off

**Generate rounded and redacted results for all drug categories
tempfile results

rounded_datatable_multi `outcome_vars', timevar(`time_variable') outfile(`results')

**Prepare final table
use "`results'", clear

forvalues i = 1/`n' {
	replace outcome_desc = "`outcome_desc`i''" if trim(outcome_name) == "outcome`i'"
	replace outcome_name = "`outcome_name`i''" if trim(outcome_name) == "outcome`i'"
}

drop if missing(month_year)

isid month_year outcome_name outcome_desc
sort outcome_name outcome_desc month_year

**Save and export final table
save "$projectdir/output/data/data_table_flares.dta", replace
export delimited using "$projectdir/output/tables/data_table_flares.csv", replace

log close
