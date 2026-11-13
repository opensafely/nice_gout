version 16

/*==============================================================================
DO FILE NAME:           Clean consultation cohort using dataset definition file
PROJECT:                OpenSAFELY NICE
AUTHOR:                 M Russell
DATASETS USED:          Consultation dataset definition
DATASETS CREATED:       Cleaned patient-level consultation cohort
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
capture mkdir "$projectdir/logs"

*Open log file
global logdir "$projectdir/logs"
cap log close
log using "$logdir/consults_cleaning.log", replace
log off

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Run checks
global checks 0
set more off
set type double

*Set disease and study dates (passed from yaml)
global arglist disease studystart_date studyend_date studyfup_date intervalstart_date
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
    global studystart_date "2016-07-01"
    global studyend_date "2024-06-30"
    global studyfup_date "2026-06-30"
    global intervalstart_date "$studyend_date"
}

*Blood test and medication defined in the consultation dataset
global blood "urate"
global medication "ult"

*Import consultation dataset
import delimited "$projectdir/output/dataset_consults.csv", clear
isid patient_id

*Conversion for dates ====================================================*

**Convert format for variables containing dates that are in string format
ds *date*, has(type string) //check list of variables is appropriate
local string_dates `r(varlist)'

foreach var of local string_dates {
	gen double `var'_num = daily(`var', "YMD")
	quietly count if !missing(`var'_num)
	if r(N) {
		format `var'_num %td
		order `var'_num, after(`var')
		drop `var'
		rename `var'_num `var'
	}
	else {
		drop `var'_num
	}
}

**Convert format for variables containing dates that are in numeric format
ds *date*, has(type numeric) //check list of variables is appropriate
capture ds *date*, has(type numeric)
if !_rc & "`r(varlist)'" != "" {
    foreach var of varlist `r(varlist)' {
        format `var' %td
    }
}

*Inclusion criteria =====================================================*

**Derive the inclusive end of the 12-calendar-month consultation interval
local interval_start = date("$intervalstart_date", "YMD")
local interval_end_month = mofd(`interval_start') + 12
local interval_end = dofm(`interval_end_month') + min(day(`interval_start'), day(dofm(`interval_end_month' + 1) - 1)) - 2

**Check criteria applied in dataset definition
gen ${disease} = !missing(${disease}_consult_date, consult_practice_id) & inrange(${disease}_consult_date, `interval_start', `interval_end') & inrange(age, 18, 110) & inlist(sex, "male", "female") & (missing(date_of_death) | date_of_death > ${disease}_consult_date)

*Generate month and year of consultation, and follow-up ===================*

gen ${disease}_year = year(${disease}_consult_date)
format ${disease}_year %ty
gen ${disease}_mon = month(${disease}_consult_date)
gen ${disease}_moyear = ym(${disease}_year, ${disease}_mon)
format ${disease}_moyear %tmMon-CCYY
lab var ${disease}_moyear "Month/Year of consultation"

**Follow-up flag: administrative data available, alive and registered
gen has_12m_fup = (followup_end_date <= date("$studyfup_date", "YMD")) & (missing(consult_reg_end_date) | consult_reg_end_date > followup_end_date) & (missing(date_of_death) | date_of_death > followup_end_date) if !missing(followup_end_date)

lab var has_12m_fup "At least 12 months of follow-up after consultation"
lab def has_12m_fup 0 "No" 1 "Yes", modify
lab val has_12m_fup has_12m_fup

*Clean and label demographic variables ==================================*

**Age
lab var age "Age at consultation"

**Sex
rename sex sex_s
encode sex_s, gen(sex)
lab def sex 1 "Female" 2 "Male", modify
lab val sex sex
lab var sex "Sex"
drop sex_s
tab sex, missing

*Baseline medication and blood test =====================================*
local Medication = upper("$medication")

**Current medication: last prescription in the six calendar months before consultation (prescription window is applied in the dataset definition)
gen ${medication}_bl = !missing(${medication}_bl_date)
lab var ${medication}_bl "`Medication' prescribed within six months before consultation"
lab def ${medication}_bl 0 "No" 1 "Yes", modify
lab val ${medication}_bl ${medication}_bl

**Baseline blood test: more recent blood test using the window already selected in the dataset definition (24-month window)
replace ${blood}_bl_value = . if !inrange(${blood}_bl_value, 0.05, 2) & !inrange(${blood}_bl_value, 50, 3000)
replace ${blood}_bl_value = . if missing(${blood}_bl_date)
replace ${blood}_bl_date = . if missing(${blood}_bl_value)
replace ${blood}_bl_value = (${blood}_bl_value * 1000) if inrange(${blood}_bl_value, 0.05, 2)
lab var ${blood}_bl_value "Serum ${blood} at baseline, micromol/L"
lab var ${blood}_bl_date "Date of baseline serum ${blood}"

gen had_${blood}_bl = !missing(${blood}_bl_value)
lab var had_${blood}_bl "Serum ${blood} performed at baseline"
lab def had_${blood}_bl 0 "No" 1 "Yes", modify
lab val had_${blood}_bl had_${blood}_bl

**Categorise baseline blood result, using the current pipeline's categories
gen ${blood}_bl_cat = 0 if ${blood}_bl_value >= 360 & !missing(${blood}_bl_value)
replace ${blood}_bl_cat = 1 if ${blood}_bl_value < 360 & !missing(${blood}_bl_value)
replace ${blood}_bl_cat = 9 if missing(${blood}_bl_value)
lab var ${blood}_bl_cat "Baseline serum ${blood} level"
lab def ${blood}_bl_cat 0 ">=360 micromol/L" 1 "<360 micromol/L" 9 "Not known", modify
lab val ${blood}_bl_cat ${blood}_bl_cat

**Agreed denominator: exclude baseline current medication AND documented target
gen controlled_bl = ${medication}_bl == 1 & ${blood}_bl_cat == 1
lab var controlled_bl "Already receiving `Medication' and documented below target at baseline"
gen eligible = controlled_bl == 0
lab var eligible "Eligible for primary target-attainment analysis"
lab def eligible 0 "Already controlled on `Medication'" 1 "Eligible", modify
lab val eligible eligible

*Relevant blood tests ===================================================*

foreach blood in $blood {

    **Find number of sequentially recorded tests, as in the current pipeline
    local max = 0
    capture quietly ds `blood'_value_*
    if !_rc & "`r(varlist)'" != "" {
        foreach v of varlist `r(varlist)' {
            if regexm("`v'", "^`blood'_value_([0-9]+)$") {
                local idx = real(regexs(1))
                local max = max(`max', `idx')
            }
        }
    }
    assert `max' > 0

    **Clean values and convert units using the current urate rules
    forvalues i = 1/`max' {
        replace `blood'_value_`i' = . if !inrange(`blood'_value_`i', 0.05, 2) & !inrange(`blood'_value_`i', 50, 3000)
        replace `blood'_value_`i' = . if missing(`blood'_date_`i')
        replace `blood'_date_`i' = . if missing(`blood'_value_`i')
        replace `blood'_value_`i' = (`blood'_value_`i' * 1000) if inrange(`blood'_value_`i', 0.05, 2)

        **Put the matched prescription number at the end for reshaping
        rename ${medication}_before_`blood'_`i'_date ${medication}_`blood'_date_`i'
    }

    **Temporarily retain only variables required for this blood test
    preserve
    keep patient_id ${disease}_consult_date followup_end_date consult_reg_end_date date_of_death `blood'_value_* `blood'_date_* ${medication}_`blood'_date_*

    **Reshape to long format
    reshape long `blood'_value_ `blood'_date_ ${medication}_`blood'_date_, i(patient_id) j(`blood'_order)

    **Keep a single test from the same day, prioritising non-missing/lower values
    bys patient_id `blood'_date_ (`blood'_value_): gen n = _n
    drop if n > 1
    drop n

    **Keep one empty row for people without a valid test
    bys patient_id (`blood'_date_): gen n = _n
    drop if n > 1 & missing(`blood'_value_)
    drop n

    **Valid follow-up tests: after consultation, within 12 calendar months
    gen `blood'_within_12m = !missing(`blood'_date_, `blood'_value_) & `blood'_date_ > ${disease}_consult_date & `blood'_date_ <= followup_end_date & (missing(consult_reg_end_date) | `blood'_date_ < consult_reg_end_date)
    lab var `blood'_within_12m "Serum `blood' measured within 12 months of consultation"
    lab def `blood'_within_12m 0 "No" 1 "Yes", modify
    lab val `blood'_within_12m `blood'_within_12m

    **Current medication at each test: prescription extracted within six months before it
    gen `blood'_on_${medication} = !missing(${medication}_`blood'_date_, `blood'_date_) & ${medication}_`blood'_date_ < `blood'_date_

    **Primary outcome: at least one result <360 while receiving medication
    gen `blood'_target_12m = `blood'_within_12m == 1 & `blood'_on_${medication} == 1 & `blood'_value_ < 360 & !missing(`blood'_value_)
    lab var `blood'_target_12m "Serum `blood' <360 on `Medication' within 12 months of consultation"
    lab def `blood'_target_12m 0 "Not attained or not tested" 1 "Target attained", modify
    lab val `blood'_target_12m `blood'_target_12m

    **Number of valid test days within the outcome window
    gen `blood'_count_12m = `blood'_within_12m
    lab var `blood'_count_12m "Number of serum `blood' tests within 12 months of consultation"

    **Store labels before collapse, as in the current pipeline
    local summaryvars `blood'_within_12m `blood'_target_12m `blood'_count_12m
    local i = 0
    foreach v of local summaryvars {
        local ++i
        local var`i' "`v'"
        local varlab`i' : variable label `v'
        local vallab`i' : value label `v'
    }
    local nvars = `i'

    **Collapse to one observation per patient
    collapse (max) `blood'_within_12m `blood'_target_12m (sum) `blood'_count_12m, by(patient_id)

    **Restore variable labels and value-label assignments
    forvalues i = 1/`nvars' {
        local v "`var`i''"
        if `"`varlab`i''"' != "" {
            label variable `v' `"`varlab`i''"'
        }
        if "`vallab`i''" != "" {
            label values `v' `vallab`i''
        }
    }

    save "$projectdir/output/data/consults_blood_summary.dta", replace
    restore

    **Merge summaries back into the main wide dataset
    merge 1:1 patient_id using "$projectdir/output/data/consults_blood_summary.dta", assert(match) nogen
}

*Final primary outcome and save =========================================*

**Leave baseline-controlled patients in the file for denominator checking
**Their primary outcome is missing, so they are not counted as eligible failures.
replace ${blood}_target_12m = . if eligible == 0

save "$projectdir/output/data/cohort_consults_processed.dta", replace
log close
