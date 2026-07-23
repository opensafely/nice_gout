version 16

/*==============================================================================
DO FILE NAME:			Produce rounded and redacted summary tables using cleaned primary cohort
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
log using "$logdir/summary_tables.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist disease comorbidities disease_features events admissions bloods medications outpatients
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
	global comorbidities "chd diabetes cva ckd hypertension depression heart_failure liver_disease transplant alcohol"
	global disease_features "tophi chronic_gout"
	global events "flare"
	global admissions "gout"
	global bloods "urate creatinine cholesterol hba1c"
	global medications "ult allopurinol febuxostat benzbromarone probenecid colchicine steroid nsaid diuretic sglt2 ace_arb"
	global outpatients "rheumatology"
}

di "$disease"
di "$comorbidities"
di "$disease_features"
di "$events"
di "$admissions"
di "$bloods"
di "$medications"
di "$outpatients"

set type double

set scheme plotplainblind

*Function to round and redact categorical variables ======================*/
program define rounded_categorical
    syntax varlist(min=1 max=1), outfile(string)
	local outcome `varlist'
	preserve 
		contract `outcome'
		local outcome_desc : variable label `outcome' 
		gen variable = `"`outcome_desc'"'
		decode `outcome', gen(categories)
		replace categories = "Missing" if categories == ""
		gen count = round(_freq, 5)
		egen total = total(count)
		gen percent = round((count/total)*100, 0.01)
		order total, before(percent)
		*egen mincount = min(count)
		replace percent =. if count<=7
		replace total =. if count<=7
		replace count =. if count<=7
		format percent %14.4f
		format count total %14.0f
		list variable categories count total percent
		keep variable categories count total percent
		capture append using `"`outfile'"'
		save `"`outfile'"', replace
    restore
end

*Function to round and redact continuous variables (amend to sum, rather than count, for ordinal variables) ======================*/
program define rounded_continuous
    syntax varlist(min=1 max=1), outfile(string)
	local outcome `varlist'
	preserve 
		local outcome_desc : variable label `outcome'
		collapse (count) count_un=`outcome' total_un=${disease} (mean) mean=`outcome' (sd) stdev=`outcome'
		gen count = round(count_un, 5)
		gen total = round(total_un, 5)
		replace stdev = . if count<=7
		replace mean = . if count<=7
		replace total = . if count<=7
		replace count = . if count<=7
		gen variable = `"`outcome_desc'"'
		order variable, first
		gen categories = "Not applicable"
		order categories, after(variable)
		order count, after(stdev)
		order total, after(count)
		format mean %14.2f
		format stdev %14.2f
		format count %14.0f
		list variable categories mean stdev count total
		keep variable categories mean stdev count total
		capture append using `"`outfile'"'
		save `"`outfile'"', replace
    restore
end

*Baseline summary table ========================*

**Store table name
local cohort "baseline"

**Erase any existing data file
capture erase "$projectdir/output/data/summary_table_`cohort'.dta"

**Load processed dataset
use "$projectdir/output/data/cohort_processed.dta", clear

**Store list of additional variables of interest
local comorbidity_vars_bl
foreach comorbidity in $comorbidities {
    local comorbidity_vars_bl `comorbidity_vars_bl' `comorbidity'_bl
}
di "`comorbidity_vars_bl'"

local disease_feature_vars_bl
foreach feature in $disease_features {
    local disease_feature_vars_bl `disease_feature_vars_bl' `feature'_bl
}
di "`disease_feature_vars_bl'"

local blood_vars_value_bl
foreach blood in $bloods {
    local blood_vars_value_bl `blood_vars_value_bl' `blood'_bl_value
}
di "`blood_vars_value_bl'"

local blood_vars_test_bl
foreach blood in $bloods {
    local blood_vars_test_bl `blood_vars_test_bl' had_`blood'_bl
}
di "`blood_vars_test_bl'"

local outpatients_ref_before ${outpatients}_ref_before
di "`outpatients_ref_before'"

local outpatients_opa_before ${outpatients}_opa_before
di "`outpatients_opa_before'"

**Loop through time periods of interest
foreach t in 12 {

	**Process catergorical outcomes of interest (other than outpatient-based variables)
	foreach outcome of varlist has_`t'm_fup ult_risk_bl sglt2_bl ace_arb_bl diuretic_bl `disease_feature_vars_bl' urate_bl_360_repeat urate_bl_cat diab_bl_cat hba1c_bl_cat ckd_transplant_bl ckd_comb_bl egfr_bl_cat `blood_vars_test_bl' `comorbidity_vars_bl' bmicat smoke region imd ethnicity sex agegroup {
		rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_`cohort'.dta")
	}
}

**Process OPA outcomes only from July 2019 onwards
preserve
    keep if ${disease}_moyear >= tm(2019m7)

    foreach outcome of varlist `outpatients_opa_before' `outpatients_ref_before' {
        rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_`cohort'.dta")
    }
restore

**Process continuous outcomes of interest
foreach outcome of varlist `blood_vars_value_bl' age {
	rounded_continuous `outcome', outfile("$projectdir/output/data/summary_table_`cohort'.dta")
}

**Export to CSV
use "$projectdir/output/data/summary_table_`cohort'.dta", clear
export delimited using "$projectdir/output/tables/summary_table_`cohort'.csv", datafmt replace

*Summary table of disease-specific events occurring within t m of diagnosis ========================*

**Store table name
local cohort "postdiagnosis"

**Loop through time periods of interest
foreach t in 12 {
	
	**Erase any existing data file
	capture erase "$projectdir/output/data/summary_table_`t'm`cohort'.dta"
	
	**Load processed dataset
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	**Set inclusion criteria - limited to those who had at least t months duration of follow-up post-diagnosis
	keep if has_`t'm_fup==1
	
	**Store list of additional variables of interest
	local blood_vars_`t'm
	foreach blood in $bloods {
		local blood_vars_`t'm `blood_vars_`t'm' `blood'_within_`t'm
	}
	di "`blood_vars_`t'm'"
	
	local outpatients_refopa_`t'm ${outpatients}_refopa_`t'm
	di "`outpatients_refopa_`t'm'"

	**Process catergorical outcomes of interest (other than outpatient-based variables)
	foreach outcome of varlist `blood_vars_`t'm' febuxostat_`t'm allopurinol_`t'm ult_first_drug_`t'm has_`t'm_fup_ult ult_cat_`t'm ult_`t'm urate_300_`t'm_cat urate_360_`t'm_cat two_urate_`t'm urate_`t'm {
		rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_`t'm`cohort'.dta")
	}
	
	**Process OPA outcome only from July 2019 onwards
	preserve
		keep if ${disease}_moyear >= tm(2019m7)
		
		rounded_categorical `outpatients_refopa_`t'm', outfile("$projectdir/output/data/summary_table_`t'm`cohort'.dta")
	
	restore

	**Process continuous outcomes of interest
	foreach outcome of varlist urate_count_`t'm lowest_urate_`t'm {
		rounded_continuous `outcome', outfile("$projectdir/output/data/summary_table_`t'm`cohort'.dta")
	}
	
	**Export to CSV
	use "$projectdir/output/data/summary_table_`t'm`cohort'.dta", clear
	export delimited using "$projectdir/output/tables/summary_table_`t'm`cohort'.csv", datafmt replace
}

*Summary table of disease-specific events occurring within t m of ULT initiation ========================*

**Store table name
local cohort "postult"

**Loop through time periods of interest
foreach t in 12 {
	
	**Erase any existing data file
	capture erase "$projectdir/output/data/summary_table_`t'm`cohort'.dta"
	
	**Load processed dataset
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	**Set inclusion criteria - limited to those who had at least t months duration of follow-up post-ULT (no restriction on ULT within 12m)
	keep if has_`t'm_fup_ult==1

	**Process catergorical outcomes of interest //removed ULT prophylaxis for now
	foreach outcome of varlist urate_`t'm_ult_cat two_urate_`t'm_ult urate_within_`t'm_ult has_`t'm_fup_target {
		rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_`t'm`cohort'.dta")
	}

	**Process continuous outcomes of interest
	foreach outcome of varlist urate_count_`t'm_ult lowest_urate_`t'm_ult {
		rounded_continuous `outcome', outfile("$projectdir/output/data/summary_table_`t'm`cohort'.dta")
	}
	
	**Export to CSV
	use "$projectdir/output/data/summary_table_`t'm`cohort'.dta", clear
	export delimited using "$projectdir/output/tables/summary_table_`t'm`cohort'.csv", datafmt replace
}

*Summary table of disease-specific events occurring within t m of urate target attainment ========================*

**Store table name
local cohort "posttarget"

**Loop through time periods of interest
foreach t in 12 {
	
	**Erase any existing data file
	capture erase "$projectdir/output/data/summary_table_`t'm`cohort'.dta"
	
	**Load processed dataset
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	**Set inclusion criteria - limited to those who had at least t months duration of follow-up post-target attainment (no restriction on ULT/target within 12m)
	keep if has_`t'm_fup_target==1
	
	**Check - for dummy data only
	count
	if r(N) == 0 {
		di as text "No observations with has_`t'm_fup_target == 1; skipping `t'-month summary table"
		continue
	}

	**Process catergorical outcomes of interest //different ULT prophylaxis to graphs
	foreach outcome of varlist repeat_below360_`t'm_ult repeat_after360_`t'm_ult {
		rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_`t'm`cohort'.dta")
	}
	
	**Export to CSV - with added check for dummy data
	capture confirm file "$projectdir/output/data/summary_table_`t'm`cohort'.dta"
	if _rc == 0 {
		use "$projectdir/output/data/summary_table_`t'm`cohort'.dta", clear
		export delimited using "$projectdir/output/tables/summary_table_`t'm`cohort'.csv", datafmt replace
	}
	else {
		di as text "No summary table created for `t'-month `cohort'; skipping export."
	}
}

log close
