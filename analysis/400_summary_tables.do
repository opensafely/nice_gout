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
global arglist disease studystart_date studyend_date studyfup_date comorbidities disease_features events admissions bloods medications outpatients
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
	global studyend_date "2025-06-30"
	global studyfup_date "2026-06-30"
	global comorbidities "chd diabetes cva ckd hypertension depression heart_failure liver_disease transplant alcohol hernia"
	global disease_features "tophi chronic_gout"
	global events "flare"
	global admissions "gout"
	global bloods "urate creatinine cholesterol hba1c"
	global medications "ult allopurinol febuxostat colchicine steroid nsaid diuretic sglt2 ace_arb"
	global outpatients "rheumatology"
}

di "$disease"
di "$studystart_date"
di "$studyend_date"
di "$studyfup_date"
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
    syntax varlist(min=1 max=1), outfile(string) [group(string)]
	if "`group'" == "" local group "Not applicable"
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
		gen str80 exposure_group = "`group'"
		order exposure_group, first
		list exposure_group variable categories count total percent
		keep exposure_group variable categories count total percent
		capture append using `"`outfile'"'
		save `"`outfile'"', replace
    restore
end

*Function to round and redact continuous variables (amend to sum, rather than count, for ordinal variables) ======================*/
program define rounded_continuous
    syntax varlist(min=1 max=1), outfile(string) [group(string)]
	if "`group'" == "" local group "Not applicable"
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
		gen str80 exposure_group = "`group'"
		order exposure_group, first
		list exposure_group variable categories mean stdev count total
		keep exposure_group variable categories mean stdev count total
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
	foreach outcome of varlist has_`t'm_fup_target has_`t'm_fup_ult ult_first_drug febuxostat_ever allopurinol_ever ult_ever has_`t'm_fup ult_risk_bl sglt2_bl ace_arb_bl diuretic_bl `disease_feature_vars_bl' urate_bl_360_repeat urate_bl_cat diab_bl_cat hba1c_bl_cat ckd_transplant_bl ckd_comb_bl egfr_bl_cat `blood_vars_test_bl' `comorbidity_vars_bl' bmicat smoke region imd ethnicity sex agegroup {
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

	**Process catergorical outcomes of interest
	foreach outcome of varlist febuxostat_ongoing_`t'm allopurinol_ongoing_`t'm ult_ongoing_`t'm urate_`t'm_ult_recode urate_`t'm_ult_cat two_urate_`t'm_ult urate_within_`t'm_ult has_`t'm_fup_target {
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

	**Process catergorical outcomes of interest
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

*Summary table of events occurring at the time of ULT initiation ========================*

**Store table name
local cohort "atultinitiation"

**Erase any existing data file
capture erase "$projectdir/output/data/summary_table_`cohort'.dta"

**Load processed dataset
use "$projectdir/output/data/cohort_processed.dta", clear

**Set inclusion criteria - limited to those who initiated ULT at any point after diagnosis, up to study end point
keep if ult_first_date !=. & (ult_first_date <= date("$studyend_date", "YMD"))

**Check - for dummy data only
count
if r(N) == 0 {
	di as text "No observations; skipping summary table"
	continue
}

**Process catergorical outcomes of interest
foreach outcome of varlist has_12m_fup_ult ult_high ult_prophylaxis_3m ult_prophylaxis ult_first_drug febuxostat_ever allopurinol_ever ult_ever {
	rounded_categorical `outcome', outfile("$projectdir/output/data/summary_table_`cohort'.dta")
}

**Export to CSV - with added check for dummy data
capture confirm file "$projectdir/output/data/summary_table_`cohort'.dta"
if _rc == 0 {
	use "$projectdir/output/data/summary_table_`cohort'.dta", clear
	export delimited using "$projectdir/output/tables/summary_table_`cohort'.csv", datafmt replace
}
else {
	di as text "No summary table created; skipping export."
}

*Summary table of events for landmark survival analyses ========================*

**Store table name
local cohort "landmark"

**Erase any existing data file
capture erase "$projectdir/output/data/summary_table_`cohort'.dta"

**Load processed dataset
use "$projectdir/output/data/cohort_processed.dta", clear

***Remove later
lab var age_land "Age, years"

**Set inclusion criteria, as per landmark criteria

***Define landmark date
local landmark_date ult_landmark //date of first ULT drug + 12 months

***Censor criteria
gen study_end = date("$studyfup_date", "YMD")
format study_end %td
local study_end_date study_end //end of study follow-up period
local death_date date_of_death //date of death
local dereg_date reg_end_date //end of practice registration
egen censor_date = rowmin(`study_end_date' `death_date' `dereg_date') //first of the above dates
format censor_date %td
label var censor_date "Censoring date"

***Outcome status at baseline/landmark variables
local outcome_free_baseline ckd_free_ult //CKD, defined using single eGFR <60 or CKD code at or before ULT initiation date
local outcome_free_landmark ckd_free_landmark //CKD, defined using single eGFR <60 or CKD code at or before ULT initiation date + 12 months

**Apply landmark eligibility criteria
keep if !missing(`landmark_date') & !missing(censor_date) & (censor_date > `landmark_date') //landmark date present and before censor date
keep if `outcome_free_baseline' ==1 //outcome not present before cohort entry
keep if `outcome_free_landmark' ==1 //outcome not present before landmark

**Apply exposure criteria
local primary_exposure urate_12m_ult_cat //separate category coded if urate not checked (1/0/9)
keep if !missing(`primary_exposure')

**Generate follow-up variable after landmark
gen landmark_12m_fup = (censor_date >= ult_landmark + 365) if !missing(ult_landmark, censor_date)
label define landmark_fup_lab 0 "Less than 365 days after landmark" 1 "At least 365 days after landmark", replace
label values landmark_12m_fup landmark_fup_lab
label variable landmark_12m_fup "Available follow-up after landmark"

**Categorical adjustment variables of interest
local categorical_vars ///
	landmark_12m_fup ace_arb_land sglt2_land diuretic_land alcohol_land hypertension_land cva_land chd_land heart_failure_land diabetes_land smoke bmicat ethnicity imd sex               

**Continuous adjustment variables of interest
local continuous_vars ///
     urate_before_ult_value egfr_before_ult_value age_land
	
local baseline_cohort "$projectdir/output/data/baseline_landmark_cohort.dta"
local baseline_results "$projectdir/output/data/summary_table_`cohort'.dta"

quietly save "`baseline_cohort'", replace
capture erase "`baseline_results'"

quietly levelsof `primary_exposure', local(exposure_levels)
local exposure_label : value label `primary_exposure'

**Loop through observed exposure groups
foreach exposure_level of local exposure_levels {

    quietly use "`baseline_cohort'", clear
    keep if `primary_exposure' == `exposure_level'

    **Default to numeric value if no label is available
    local group_label "`exposure_level'"

    if "`exposure_label'" != "" {
        local group_label : label `exposure_label' `exposure_level'
    }

    **Categorical variables
    foreach var of local categorical_vars {
        rounded_categorical `var', outfile("`baseline_results'") group("`group_label'")
    }

    **Continuous variables
    foreach var of local continuous_vars {
        rounded_continuous `var', outfile("`baseline_results'") group("`group_label'")
    }
	
	**Monitoring summaries: require a full year after landmark
    keep if landmark_12m_fup == 1

    quietly count
    if r(N) > 0 {

        foreach var in creat_any_12m_land creat_two_12m_land {
            rounded_categorical `var', outfile("`baseline_results'") group("`group_label'")
        }

        rounded_continuous creat_n_12m_land, outfile("`baseline_results'") group("`group_label'")
    }
}

**Export to CSV - with check for dummy data
capture confirm file "`baseline_results'"
if _rc == 0 {
    use "`baseline_results'", clear
    order exposure_group variable categories
    export delimited using "$projectdir/output/tables/summary_table_`cohort'.csv", datafmt replace
}
else {
    di as text "No summary table created; skipping export."
}

*Landmark table for flowchart of eligibility

**Store table name
local cohort "landmark_inclusion"

**Erase existing output
local results "$projectdir/output/data/summary_table_`cohort'.dta"
capture erase "`results'"

**Load processed dataset
use "$projectdir/output/data/cohort_processed.dta", clear

**Include everyone who started ULT
keep if !missing(ult_first_date)

**Alive at ULT plus 12 months
gen alive_landmark = (missing(date_of_death) | date_of_death > ult_landmark) if !missing(ult_landmark)

**Deregistered on or before ULT plus 12 months
gen deregistered_landmark = (!missing(reg_end_date) & reg_end_date <= ult_landmark) if !missing(ult_landmark)

**CKD at ULT initiation
gen ckd_baseline = (ckd_free_ult == 0) if !missing(ckd_free_ult)

**CKD free at landmark
gen ckd_free_at_landmark = ckd_free_landmark if !missing(ult_landmark)

****Meets all inclusion criteria
gen eligible_landmark = !missing(ult_landmark) & (ult_landmark < date("$studyfup_date", "YMD")) & has_12m_fup_ult == 1 & alive_landmark == 1 & deregistered_landmark == 0 & ckd_baseline == 0 & ckd_free_at_landmark == 1

**Labels
label define inclusion_yesno 0 "No" 1 "Yes", replace

foreach var in alive_landmark deregistered_landmark ckd_baseline ckd_free_at_landmark eligible_landmark {
    label values `var' inclusion_yesno
}

label var alive_landmark "Alive at ULT plus 12 months"
label var deregistered_landmark "Deregistered by ULT plus 12 months"
label var ckd_baseline "CKD at ULT initiation"
label var ckd_free_at_landmark "CKD free at ULT plus 12 months"
label var eligible_landmark "Meets all landmark inclusion criteria"

**Output each characteristic among all ULT initiators
foreach var in has_12m_fup_ult alive_landmark deregistered_landmark ckd_baseline ckd_free_at_landmark eligible_landmark {
    rounded_categorical `var', outfile("`results'") group("All ULT initiators")
}

**Export to CSV
capture confirm file "`results'"

if _rc == 0 {
    use "`results'", clear
    order exposure_group variable categories

    export delimited using "$projectdir/output/tables/summary_table_`cohort'.csv", datafmt replace
}

log close
