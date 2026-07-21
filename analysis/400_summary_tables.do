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
		gen count = round(_freq, 5)
		egen total = total(count)
		gen percent = round((count/total)*100, 0.01)
		order total, before(percent)
		egen mincount = min(count)
		replace percent = . if mincount<=7
		replace total = . if mincount<=7
		replace count = . if mincount<=7
		drop mincount
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

	**Process catergorical outcomes of interest //different ULT prophylaxis to graphs
	foreach outcome of varlist urate_`t'm_ult_cat two_urate_`t'm_ult urate_within_`t'm_ult {
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

/*Generate unrounded tables for checking (not outputted) ==============================*

**Baseline table

***Load processed dataset
use "$projectdir/output/data/cohort_processed.dta", clear

***Store list of additional variables of interest (passed from yaml)
local comorbidity_vars_bl
foreach comorbidity in $comorbidities {
    if "`comorbidity_vars_bl'" == "" {
        local comorbidity_vars_bl "`comorbidity'_bl cat %5.1f"
    }
    else {
        local comorbidity_vars_bl "`comorbidity_vars_bl' \ `comorbidity'_bl cat %5.1f"
    }
}
di "`comorbidity_vars_bl'"

local disease_feature_vars_bl
foreach feature in $disease_features {
    if "`disease_feature_vars_bl'" == "" {
        local disease_feature_vars_bl "`feature'_bl cat %5.1f"
    }
    else {
        local disease_feature_vars_bl "`disease_feature_vars_bl' \ `feature'_bl cat %5.1f"
    }
}
di "`disease_feature_vars_bl'"

local blood_vars_bl
foreach blood in $bloods {
    if "`blood_vars_bl'" == "" {
        local blood_vars_bl "had_`blood'_bl cat %5.1f \ `blood'_bl_value contn %5.1f"
    }
    else {
        local blood_vars_bl "`blood_vars_bl' \ had_`blood'_bl cat %5.1f \ `blood'_bl_value contn %5.1f"
    }
}
di "`blood_vars_bl'"

preserve
table1_mc, total(before) onecol nospacelowpercent missing iqrmiddle(",")  ///
	vars(age contn %5.1f \ ///
		 agegroup cat %5.1f \ ///
		 sex cat %5.1f \ ///
		 ethnicity cat %5.1f \ ///
		 imd cat %5.1f \ ///
		 region cat %5.1f \ ///
		 smoke cat %5.1f \ ///
		 bmicat cat %5.1f \ ///
         `comorbidity_vars_bl' \ ///
		 `blood_vars_bl' \ ///
		 egfr_bl_cat cat %5.1f \ ///
		 ckd_comb_bl cat %5.1f \ ///
		 ckd_transplant_bl cat %5.1f \ ///
		 hba1c_bl_cat cat %5.1f \ ///
		 diab_bl_cat cat %5.1f \ ///
		 urate_bl_cat cat %5.1f \ ///
		 urate_bl_360_repeat cat %5.1f \ ///
		 `disease_feature_vars_bl' \ ///
		 diuretic_bl cat %5.1f \ ///
		 sglt2_bl cat %5.1f \ ///
		 ace_arb_bl cat %5.1f \ ///
		 ult_risk_bl cat %5.1f \ ///
		 ${outpatients}_ref_before cat %5.1f \ ///
		 ${outpatients}_opa_before cat %5.1f \ ///
		 has_6m_fup cat %5.1f \ ///
		 has_12m_fup cat %5.1f \ ///
		 )
restore

/*
preserve
table1_mc, total(before) by(diagnosis_year) onecol nospacelowpercent missing iqrmiddle(",")  ///
	vars(age contn %5.1f \ ///
		 agegroup cat %5.1f \ ///	 
		 sex cat %5.1f \ ///
		 ethnicity cat %5.1f \ ///
		 imd cat %5.1f \ ///
		 region cat %5.1f \ ///
		 smoke cat %5.1f \ ///
		 bmicat cat %5.1f \ ///
         `comorbidity_vars_bl' \ ///
		 `blood_vars_bl' \ ///
		 egfr_bl_cat cat %5.1f \ ///
		 ckd_comb_bl cat %5.1f \ ///
		 ckd_transplant_bl cat %5.1f \ ///
		 hba1c_bl_cat cat %5.1f \ ///
		 diab_bl_cat cat %5.1f \ ///
		 urate_bl_cat cat %5.1f \ ///
		 urate_bl_360_repeat cat %5.1f \ ///
		 `disease_feature_vars_bl' \ ///
		 diuretic_bl cat %5.1f \ ///
		 ult_risk_bl cat %5.1f \ ///
		 ${outpatients}_ref_before cat %5.1f \ ///
		 ${outpatients}_opa_before cat %5.1f \ ///
		 has_6m_fup cat %5.1f \ ///
		 has_12m_fup cat %5.1f \ ///
		 )
restore
*/

**Disease-specific events occurring within 6/12m of diagnosis

***Loop through time periods of interest
foreach t in 12 {
	
	***Load processed dataset
	use "$projectdir/output/data/cohort_processed.dta", clear

	***Set inclusion criteria - limited to those who had at least t months duration of follow-up post-diagnosis
	keep if has_`t'm_fup==1
	
	***Store list of additional variables of interest (passed from yaml)
	local blood_vars_`t'm
	foreach blood in $bloods {
		if "`blood_vars_`t'm'" == "" {
			local blood_vars_`t'm "`blood'_within_`t'm cat %5.1f"
		}
		else {
			local blood_vars_`t'm  "`blood_vars_`t'm' \ `blood'_within_`t'm cat %5.1f"
		}
	}
	di "`blood_vars_`t'm'"
	
	preserve
	table1_mc, total(before) onecol nospacelowpercent missing iqrmiddle(",")  ///
		vars(urate_`t'm cat %5.1f \ ///
			 two_urate_`t'm cat %5.1f \ ///
			 urate_count_`t'm contn %5.1f \ ///
			 urate_360_`t'm_cat cat %5.1f \ ///
			 urate_300_`t'm_cat cat %5.1f \ ///
			 lowest_urate_`t'm contn %5.1f \ ///
			 ult_`t'm cat %5.1f \ ///
			 ult_cat_`t'm cat %5.1f \ ///
			 has_`t'm_fup_ult cat %5.1f \ ///
			 ult_first_drug_`t'm cat %5.1f \ ///
			 allopurinol_`t'm cat %5.1f \ ///
			 febuxostat_`t'm cat %5.1f \ ///
			 ${outpatients}_refopa_`t'm cat %5.1f \ ///
			 `blood_vars_`t'm' cat %5.1f \ ///
			 )
	restore	
/*	
	preserve
	keep if has_`t'm_fup==1
	table1_mc, total(before) by(diagnosis_year) onecol nospacelowpercent missing iqrmiddle(",")  ///
		vars(urate_`t'm cat %5.1f \ ///
			 two_urate_`t'm cat %5.1f \ ///
			 urate_count_`t'm contn %5.1f \ ///
			 urate_360_`t'm_cat cat %5.1f \ ///
			 urate_300_`t'm_cat cat %5.1f \ ///
			 lowest_urate_`t'm contn %5.1f \ ///
			 ult_`t'm cat %5.1f \ ///
			 ult_cat_`t'm cat %5.1f \ ///
			 has_`t'm_fup_ult cat %5.1f \ ///
			 ult_first_drug_`t'm cat %5.1f \ ///
			 allopurinol_`t'm cat %5.1f \ ///
			 febuxostat_`t'm cat %5.1f \ ///
			 ${outpatients}_refopa_`t'm cat %5.1f \ ///
			 `blood_vars_`t'm' cat %5.1f \ ///
			 )
	restore	
*/
}
	
**Disease-specific events occurring within 6/12m of ULT initiation

***Loop through time periods of interest
foreach t in 12 {
	
	***Load processed dataset
	use "$projectdir/output/data/cohort_processed.dta", clear

	***Set inclusion criteria - limited to those who initiated ULT and who had at least t months duration of follow-up post-diagnosis
	keep if ult_`t'm & has_`t'm_fup_ult==1
	
	preserve
	table1_mc, total(before) onecol nospacelowpercent missing iqrmiddle(",")  ///
		vars(urate_within_`t'm_ult cat %5.1f \ ///
			 two_urate_`t'm_ult cat %5.1f \ ///
			 urate_count_`t'm_ult contn %5.1f \ ///
			 urate_`t'm_ult_cat cat %5.1f \ ///
			 lowest_urate_`t'm_ult contn %5.1f \ ///
			 ult_scripts_`t'm contn %5.1f \ ///
			 ult_prophylaxis_`t'm cat %5.1f \ ///
			 repeat_after360_`t'm_ult cat %5.1f \ ///
			 repeat_below360_`t'm_ult cat %5.1f \ ///
			 )
	restore	
/*	
	preserve
	keep if ult_`t'm & has_`t'm_fup_ult==1
	table1_mc, total(before) by(diagnosis_year) onecol nospacelowpercent missing iqrmiddle(",")  ///
		vars(urate_within_`t'm_ult cat %5.1f \ ///
			 two_urate_`t'm_ult cat %5.1f \ ///
			 urate_count_`t'm_ult contn %5.1f \ ///
			 urate_`t'm_ult_cat cat %5.1f \ ///
			 lowest_urate_`t'm_ult contn %5.1f \ ///
			 ult_scripts_`t'm contn %5.1f \ ///
			 ult_prophylaxis cat %5.1f \ ///
			 repeat_after360_ult cat %5.1f \ ///
			 repeat_below360_ult cat %5.1f \ /// 
			 )
	restore	
*/
}
*/
log close
