version 16

/*==============================================================================
DO FILE NAME:			Clean primary cohort using dataset definition file
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Primary dataset defintion
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
log using "$logdir/cohort_cleaning.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list, characteristics of interest, and study dates (passed from yaml)
global arglist disease studystart_date studyend_date studyfup_date nice_date demographic comorbidities disease_features events admissions bloods medications outpatients
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
	global studyfup_date "2025-12-31"
	global nice_date "2022-06-01"
	global demographic "agegroup sex ethnicity imd region"
	global comorbidities "chd diabetes cva ckd hypertension depression heart_failure liver_disease transplant alcohol"
	global disease_features "tophi chronic_gout"
	global events "flare"
	global admissions "gout"
	global bloods "urate creatinine cholesterol hba1c"
	global medications "ult allopurinol allopurinol_high febuxostat febuxostat_high benzbromarone probenecid colchicine steroid nsaid diuretic sglt2 ace_arb"
	global outpatients "rheumatology"
}

di "$disease"
di "$studystart_date"
di "$studyend_date"
di "$studyfup_date"
di "$nice_date"
di "$demographic"
di "$comorbidities"
di "$disease_features"
di "$events"
di "$admissions"
di "$bloods"
di "$medications"
di "$outpatients"

*Start year, end year and number of study years (derived from above)
global base_year = year(date("$studystart_date", "YMD"))
global end_year = year(date("$studyend_date", "YMD"))
global max_year = $end_year - $base_year
di "$base_year"
di "$end_year"
di "$max_year"

*Define max number of sequentially recorded prescriptions for an example drug (same for each drug in dataset definition)
local max_prescription = 0
local first_med: word 1 of $medications
display "`first_med'"
capture quietly ds `first_med'_date_*
if !_rc & "`r(varlist)'" != "" {
	foreach v of varlist `r(varlist)' {
		if regexm("`v'","^`first_med'_date_([0-9]+)$") {
			local idx = real(regexs(1))
			local max_prescription = max(`max_prescription', `idx')
		}
	}
}
global max_prescription = `max_prescription'
di "$max_prescription"

set type double

set scheme plotplainblind

*Import primary dataset
import delimited "$projectdir/output/dataset_primary.csv", clear

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

*Inclusion criteria (amend as necessary) ===============================================*/

di "${disease}_inc_case"
codebook ${disease}_inc_case

**Check criteria applied in dataset definition
gen ${disease} = 1 if ${disease}_inc_case=="T" & ((${disease}_inc_date >= date("$studystart_date", "YMD")) & (${disease}_inc_date <= date("$studyend_date", "YMD"))) & (${disease}_age >=18 & ${disease}_age <= 110) & (sex!="") & ${disease}_pre_reg=="T" & ${disease}_alive_inc=="T"
recode ${disease} .=0
tab ${disease}, missing //should all be 1
codebook ${disease}_inc_date //check

**Apply additional disease-specific inclusion criteria (amend as necessary)

**Drop cases where first ULT prescription was >30 days before first recorded diagnosis code (could also apply in dataset definition, but retained here to assess how many are dropped)
codebook ult_first_date //check
gen ult_first = 1 if ult_first_date!=.
tab ult_first if ult_first_date!=. & ult_first_date<${disease}_inc_date
tab ult_first if ult_first_date!=. & (ult_first_date+30)<${disease}_inc_date //>30 days before first recorded diagnosis code
tab ult_first if ult_first_date!=. & (ult_first_date+60)<${disease}_inc_date //>60 days before first recorded diagnosis code
drop if ${disease}_inc_date!=. & ult_first_date!=. & (ult_first_date+30)<${disease}_inc_date //drop if first ULT script more than 30 days before first recorded diagnosis code 
codebook ${disease}_inc_date //check

*Generate month and year of diagnosis, and define duration of follow-up post-diagnosis ===========================================================*/

gen ${disease}_year = year(${disease}_inc_date)
format ${disease}_year %ty
gen ${disease}_mon = month(${disease}_inc_date)
gen ${disease}_moyear = ym(${disease}_year, ${disease}_mon)
format ${disease}_moyear %tmMon-CCYY
lab var ${disease}_moyear "Month/Year of diagnosis"

**Separate into 12-month time windows from July (with relation to diagnosis date)
gen diagnosis_year = (floor((${disease}_inc_date - date("$studystart_date", "YMD")) / 365.25) + 1) if inrange(${disease}_inc_date, date("$studystart_date", "YMD"), date("$studyend_date", "YMD"))
lab var diagnosis_year "Year of diagnosis"
forvalues i = 1/$max_year {
    local start = $base_year + `i' - 1
    local end = `start' + 1
    label define diagnosis_year_lbl `i' "July `start'-June `end'", add
}
lab val diagnosis_year diagnosis_year_lbl
tab diagnosis_year, missing

**Proportion of patients with at least 6 or 12 months of GP registration after diagnosis
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	gen has_`t'm_fup=1 if (reg_end_date!=. & (reg_end_date >= (${disease}_inc_date + `days')) & ((${disease}_inc_date + `days') <= (date("$studyfup_date", "YMD")))) | (reg_end_date==. & ((${disease}_inc_date + `days') <= (date("$studyfup_date", "YMD"))))
	recode has_`t'm_fup .=0
	lab var has_`t'm_fup "At least `t' months of follow-up after diagnosis"
	lab def has_`t'm_fup 0 "No" 1 "Yes"
	lab val has_`t'm_fup has_`t'm_fup
	tab has_`t'm_fup
	tab ${disease}_moyear has_`t'm_fup
}

**Generate variable for date of diagnosis before vs. after NICE guideline publication
gen post_nice_diag = 1 if ${disease}_inc_date >= date("$nice_date","YMD")
recode post_nice_diag .=0
lab var post_nice_diag "Diagnosed after NICE guideline"
lab def post_nice_diag 0 "No" 1 "Yes"
lab val post_nice_diag post_nice_diag
tab post_nice_diag

*Clean and label demographic and comorbidity variables ====================================*/

**Age
rename ${disease}_age age
lab var age "Age at diagnosis"
tabstat age, stat(n mean sd p50 p25 p75)

gen age_decile = age/10
lab var age_decile "Age decile at diagnosis"

***Define 10-year age bands
recode age 18/29.9999 = 1 /// 
		   30/39.9999 = 2 ///
           40/49.9999 = 3 ///
		   50/59.9999 = 4 ///
	       60/69.9999 = 5 ///
		   70/79.9999 = 6 ///
		   80/max = 7, gen(agegroup_broad) 

label define agegroup_broad	1 "18 to 29" ///
							2 "30 to 39" ///
							3 "40 to 49" ///
							4 "50 to 59" ///
							5 "60 to 69" ///
							6 "70 to 79" ///
							7 "80 or above", modify
						
label values agegroup_broad agegroup_broad
lab var agegroup_broad "Age group"
order agegroup_broad, after(age)
tab agegroup_broad, missing

***Define broader age bands
recode age 18/39.9999 = 1 /// 
		   40/59.9999 = 2 ///
           60/79.9999 = 3 ///
		   80/max = 4, gen(agegroup) 

label define agegroup 	1 "18 to 39" ///
						2 "40 to 59" ///
						3 "60 to 79" ///
						4 "80 or above", modify
						
label values agegroup agegroup
lab var agegroup "Age group"
order agegroup, after(age)
tab agegroup, missing

**Sex
rename sex sex_s
encode sex_s, gen(sex)
lab def sex 1 "Female" 2 "Male", modify
lab val sex sex
lab var sex "Sex"
drop sex_s
tab sex, missing

**Ethnicity
gen ethnicity_n = 1 if ethnicity == "White"
replace ethnicity_n = 2 if ethnicity == "Asian or Asian British"
replace ethnicity_n = 3 if ethnicity == "Black or Black British"
replace ethnicity_n = 4 if ethnicity == "Mixed"
replace ethnicity_n = 5 if ethnicity == "Chinese or Other Ethnic Groups"
replace ethnicity_n = 9 if ethnicity == "Unknown"

label define ethnicity_n	1 "White"  						///
							2 "Asian"		///
							3 "Black"  	///
							4 "Mixed"						///
							5 "Chinese or Other" ///
							9 "Not known", modify
							
label val ethnicity_n ethnicity_n
lab var ethnicity_n "Ethnicity"
drop ethnicity
rename ethnicity_n ethnicity 
tab ethnicity, missing

**Practice region (at time of primary diagnosis)
replace region="Not known" if region==""
replace region="Yorkshire Humber" if region=="Yorkshire and The Humber"
encode region, gen(nuts_region)
replace nuts_region = 9 if region=="Not known"
drop region
rename nuts_region region
lab var region "Region"
tab region, missing

**Index of multiple deprivation (at time of primary diagnosis)
gen imd = 1 if imd_quintile == "1 (most deprived)"
replace imd = 2 if imd_quintile == "2"
replace imd = 3 if imd_quintile == "3"
replace imd = 4 if imd_quintile == "4"
replace imd = 5 if imd_quintile == "5 (least deprived)"
replace imd = 9 if imd_quintile == "Unknown"

label define imd 1 "1 most deprived" 2 "2" 3 "3" 4 "4" 5 "5 least deprived" 9 "Not known", modify
label val imd imd 
lab var imd "Index of multiple deprivation"
drop imd_quintile
tab imd, missing

**Index of multiple deprivation (lastest address) - sense check; remove later
gen imd_latest = 1 if imd_quintile_latest == "1 (most deprived)"
replace imd_latest = 2 if imd_quintile_latest == "2"
replace imd_latest = 3 if imd_quintile_latest == "3"
replace imd_latest = 4 if imd_quintile_latest == "4"
replace imd_latest = 5 if imd_quintile_latest == "5 (least deprived)"
replace imd_latest = 9 if imd_quintile_latest == "Unknown"

label define imd_latest 1 "1 most deprived" 2 "2" 3 "3" 4 "4" 5 "5 least deprived" 9 "Not known", modify
label val imd_latest imd_latest 
lab var imd_latest "Index of multiple deprivation"
drop imd_quintile_latest
tab imd_latest, missing

**Body Mass Index
***Recode values that are more likely to be erroneous
replace bmi_value = . if !inrange(bmi_value, 10, 80)

***Restrict to BMI recorded within 10 years of primary diagnosis date and aged > 16 years old
gen bmi_time = (${disease}_inc_date - bmi_date)/365.25
gen bmi_age = age - bmi_time
replace bmi_value = . if bmi_age < 16 
replace bmi_value = . if bmi_time > 10 & bmi_time !=. 
replace bmi_value = . if bmi_date == . 
replace bmi_date = . if bmi_value == . 
replace bmi_time = . if bmi_value == . 
replace bmi_age = . if bmi_value == . 

***Create BMI categories
gen bmicat = .
recode bmicat . = 1 if bmi_value < 18.5
recode bmicat . = 2 if bmi_value < 25
recode bmicat . = 3 if bmi_value < 30
recode bmicat . = 4 if bmi_value < 35
recode bmicat . = 5 if bmi_value < 40
recode bmicat . = 6 if bmi_value >= 40 & bmi_value!=.
replace bmicat = 9 if bmi_value == .

label define bmicat 1 "Underweight (<18.5)" 	///
					2 "Normal (18.5-24.9)"		///
					3 "Overweight (25-29.9)"	///
					4 "Obese I (30-34.9)"		///
					5 "Obese II (35-39.9)"		///
					6 "Obese III (40+)"			///
					9 "Not known"
					
label values bmicat bmicat
lab var bmicat "BMI"
order bmicat, after (bmi_value)
drop bmi_age bmi_time
tab bmicat, missing

**Smoking status
gen smoke = 1 if smoking_status == "N"
replace smoke = 2 if smoking_status == "E"
replace smoke = 3 if smoking_status == "S"
replace smoke = 9 if smoking_status == "M"
replace smoke = 9 if smoking_status == "" 
label define smoke 1 "Never" 2 "Former" 3 "Current" 9 "Not known"
label values smoke smoke
lab var smoke "Smoking status"
drop smoking_status ever_smoked most_recent_smoking_code
tab smoke, missing

***Create non-missing 3-category variable for current smoking (assumes missing smoking is never smoking)
recode smoke 9 = 1, gen(smoke_nomiss)
order smoke_nomiss, after(smoke)
label values smoke_nomiss smoke
lab var smoke_nomiss "Smoking status"
tab smoke_nomiss, missing

**Clinical comorbidities at baseline and after diagnosis (using recorded codes only for now)
foreach comorbidity in $comorbidities {
    local lbl : subinstr local comorbidity "_" " ", all
	local lbl = strproper("`lbl'")
	di "`lbl'"
	gen `comorbidity'_bl = 1 if (`comorbidity'_date <= ${disease}_inc_date) & `comorbidity'_date!=.
	recode `comorbidity'_bl .=0
	lab define `comorbidity'_bl 0 "No" 1 "Yes", modify
	lab val `comorbidity'_bl `comorbidity'_bl
	lab var `comorbidity'_bl "`lbl'"
	order `comorbidity'_bl, after(`comorbidity'_date)
	tab `comorbidity'_bl, missing

	gen `comorbidity'_new = 1 if (`comorbidity'_date > ${disease}_inc_date) & `comorbidity'_date!=.
	recode `comorbidity'_new .=0
	lab define `comorbidity'_new 0 "No" 1 "Yes", modify //ignores baseline disease
	lab var `comorbidity'_new "`lbl'"
	lab val `comorbidity'_new `comorbidity'_new
	order `comorbidity'_new, after(`comorbidity'_bl)
	tab `comorbidity'_new, missing
	
	gen `comorbidity'_12m = 1 if (`comorbidity'_date <= ${disease}_inc_date) | ((`comorbidity'_date > ${disease}_inc_date) & (`comorbidity'_date <= (${disease}_inc_date + 365))) & `comorbidity'_date!=.
	recode `comorbidity'_12m .=0
	lab define `comorbidity'_12m 0 "No" 1 "Yes", modify
	lab var `comorbidity'_12m "`lbl'"
	lab val `comorbidity'_12m `comorbidity'_12m
	order `comorbidity'_12m, after(`comorbidity'_new)
	tab `comorbidity'_12m, missing
}

***Re-label variables with acronyms (amend manually)
lab var ckd_bl "CKD"
lab var ckd_new "CKD"
lab var ckd_12m "CKD"
lab var diabetes_bl "T2DM"
lab var diabetes_new "T2DM"
lab var diabetes_12m "T2DM"
lab var chd_bl "CHD"
lab var chd_new "CHD"
lab var chd_12m "CHD"
lab var cva_bl "Stroke/TIA"
lab var cva_new "Stroke/TIA"
lab var cva_12m "Stroke/TIA"
lab var liver_disease_bl "Chronic liver disease"
lab var liver_disease_new "Chronic liver disease"
lab var liver_disease_12m "Chronic liver disease"
lab var transplant_bl "Solid organ transplant"
lab var transplant_new "Solid organ transplant"
lab var transplant_12m "Solid organ transplant"
lab var alcohol_bl "Excess alcohol"
lab var alcohol_new "Excess alcohol"
lab var alcohol_12m "Excess alcohol"

**Disease-specific features at baseline and after diagnosis (passed from yaml) =================================*/

foreach feature in $disease_features {
    local lbl : subinstr local feature "_" " ", all
	local lbl = strproper("`lbl'")
	di "`lbl'"
	gen `feature'_bl = 1 if (`feature'_date <= ${disease}_inc_date) & `feature'_date!=.
	recode `feature'_bl .=0
	lab define `feature'_bl 0 "No" 1 "Yes", modify
	lab val `feature'_bl `feature'_bl
	lab var `feature'_bl "`lbl'"
	order `feature'_bl, after(`feature'_date)
	tab `feature'_bl, missing
	
	gen `feature'_new = 1 if (`feature'_date > ${disease}_inc_date) & `feature'_date!=.
	recode `feature'_new .=0
	lab define `feature'_new 0 "No" 1 "Yes", modify //ignores baseline disease
	lab var `feature'_new "`lbl'"
	lab val `feature'_new `feature'_new
	order `feature'_new, after(`feature'_bl)
	tab `feature'_new, missing
}

save "$projectdir/output/data/cohort_generic.dta", replace

*Disease-specific medications (amend as necessary) =================================*/

use "$projectdir/output/data/cohort_generic.dta", clear

**Specify drug of interest: diuretics =============*/
foreach drug in diuretic sglt2 ace_arb {

	local Drug = strproper("`drug'") //first letter capitalised for labelling
	di "`drug'"

	***Drug use at baseline (prescribed within 6 months prior to index diagnosis)
	codebook `drug'_bl_date //check then remove
	
	gen `drug'_bl = 1 if (${disease}_inc_date - `drug'_bl_date) <= 183 & `drug'_bl_date != . //not necessary
	recode `drug'_bl .=0
	lab define `drug'_bl 0 "No" 1 "Yes", modify
	lab val `drug'_bl `drug'_bl
	lab var `drug'_bl "`Drug' at baseline"
}

**Specify drug of interest: urate-lowering therapy (ULT) =============*/
local drug "ult"

local Drug = upper("`drug'") //all caps for labelling
di "`drug'"

***Separate first prescription into 12-month time windows from July
gen `drug'_year = (floor((`drug'_first_date - date("$studystart_date", "YMD")) / 365.25) + 1) if inrange(`drug'_first_date, date("$studystart_date", "YMD"), date("$studyend_date", "YMD"))
lab var `drug'_year "Year of first `Drug' prescription"
forvalues i = 1/$max_year {
    local start = $base_year + `i' - 1
    local end = `start' + 1
    label define `drug'_year_lbl `i' "July `start'-June `end'", add
}
lab val `drug'_year `drug'_year_lbl
tab `drug'_year, missing

***What was the first drug prescribed within a class (amend list as necessary)
gen `drug'_first_d="" if `drug'_first_date!=.
foreach var of varlist allopurinol_first_date febuxostat_first_date benzbromarone_first_date probenecid_first_date {
	replace `drug'_first_d="`var'" if `drug'_first_date==`var' & `drug'_first_date!=. & `var'!=.
	}
gen `drug'_first_drug_s = strproper(substr(`drug'_first_d, 1, strpos(`drug'_first_d, "_") - 1)) if `drug'_first_d!=""  
encode `drug'_first_drug_s, gen(`drug'_first_drug)
lab var `drug'_first_drug "First prescribed `Drug' drug"
tab `drug'_first_drug, missing
drop `drug'_first_d

***Prescriptions of drug class within a timeframe after diagnosis
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	****First drug prescribed within 6/12m after diagnosis
	gen `drug'_first_drug_`t'm_s = `drug'_first_drug_s if (`drug'_first_date <= (${disease}_inc_date + `days')) & `drug'_first_date!=.
	encode `drug'_first_drug_`t'm_s, gen(`drug'_first_drug_`t'm)
	lab var `drug'_first_drug_`t'm "First prescribed `Drug' drug within `t' months of diagnosis"
	tab `drug'_first_drug_`t'm, missing
	drop `drug'_first_drug_`t'm_s

	****Proportion of patients with at least 6/12 months of GP registration after first prescription - can subsequently limit this to those who initiated within 6/12m of diagnosis 
	gen has_`t'm_fup_`drug'=1 if (reg_end_date!=. & (reg_end_date >= (`drug'_first_date + `days')) & ((`drug'_first_date + `days') <= (date("$studyfup_date", "YMD")))) | (reg_end_date==. & ((`drug'_first_date + `days') <= (date("$studyfup_date", "YMD"))))
	recode has_`t'm_fup_`drug' .=0 //includes those who didn't initiate drug
	lab var has_`t'm_fup_`drug' "At least `t' months of follow-up after `Drug' initiation"
	lab def has_`t'm_fup_`drug' 0 "No" 1 "Yes"
	lab val has_`t'm_fup_`drug' has_`t'm_fup_`drug'
	tab has_`t'm_fup_`drug'
	tab ${disease}_moyear has_`t'm_fup_`drug'
}

**Generate variable for date of first drug with relation to NICE guideline publication
gen post_nice_`drug' = 1 if `drug'_first_date >= date("$nice_date","YMD") & `drug'_first_date!=.
replace post_nice_`drug' = 0 if `drug'_first_date < date("$nice_date","YMD") & `drug'_first_date!=.
lab var post_nice_`drug' "`Drug' initiated after NICE guideline"
lab def post_nice_`drug' 0 "No" 1 "Yes"
lab val post_nice_`drug' post_nice_`drug'
tab post_nice_`drug'

***Was the first prescription for a drug a high dose (amend list as necessary)
gen `drug'_high = 0 if `drug'_first_date!=.

foreach var of varlist allopurinol_high_first_date febuxostat_high_first_date {
	replace `drug'_high=1 if `drug'_first_date==`var' & `drug'_first_date!=. & `var'!=.
}
	lab var `drug'_high "First prescribed `Drug' was high dose"
	lab define `drug'_high 0 "No" 1 "Yes"
	lab val `drug'_high `drug'_high 
	tab `drug'_high, missing

***Prescriptions for individual drugs within class within a timeframe after diagnosis (amend list as necessary)
foreach med in `drug' allopurinol febuxostat benzbromarone probenecid  {
	if "`med'" == "`drug'" {
		local druglabel = "`Drug'" 
	}
	else {
		local druglabel = strproper("`med'")
	}
	
	****Ever prescribed drug before study end date
	gen `med'_ever = 1 if `med'_first_date!=. 
	recode `med'_ever .=0
	lab var `med'_ever "Ever prescribed `druglabel'"
	lab define `med'_ever 0 "No" 1 "Yes"
	lab val `med'_ever `med'_ever

	***Was prophylaxis prescribed on same day at first ULT drug - Nb. restricted to first initiation of each ULT drug
	gen `med'_prophylaxis = .
	forval i = 1/$max_prescription {
		display $max_prescription
		display `i'
	    replace `med'_prophylaxis = 1 if (`med'_first_date == colchicine_date_`i') & colchicine_date_`i'!=. & `med'_first_date!=.
		replace `med'_prophylaxis = 1 if (`med'_first_date == nsaid_date_`i') & nsaid_date_`i'!=. & `med'_first_date!=.
		replace `med'_prophylaxis = 1 if (`med'_first_date == steroid_date_`i') & steroid_date_`i'!=. & `med'_first_date!=.
	}
	recode `med'_prophylaxis .=0 if `med'_first_date !=.
	lab var `med'_prophylaxis "Prophylaxis prescribed at same time as ULT initiation"
	lab define `med'_prophylaxis 0 "No" 1 "Yes"
	lab val `med'_prophylaxis `med'_prophylaxis

	****Time from diagnosis to first prescription
	gen time_to_`med' = (`med'_first_date - ${disease}_inc_date) if `med'_first_date!=.
	
	****Generate variable for those who had prescription within 6/12m of diagnosis
	foreach t in 6 12 {
		local days = int((`t'/12)*365.25)
		di `days'
	 
		gen `med'_`t'm = 1 if time_to_`med'<=`days' & time_to_`med'!=.
		recode `med'_`t'm .=0
		lab var `med'_`t'm "`druglabel' prescribed within `t' months of diagnosis"
		lab define `med'_`t'm 0 "No" 1 "Yes", modify
		lab val `med'_`t'm `med'_`t'm
		tab `med'_`t'm, missing
		tab `med'_`t'm if has_`t'm_fup==1, missing //for those with at least 6/12m of available follow-up after diagnosis
		
		****Was prophylaxis prescribed on same day at first ULT drug
		gen `med'_prophylaxis_`t'm = `med'_prophylaxis if `med'_`t'm == 1
		lab var `med'_prophylaxis_`t'm "Prophylaxis prescribed at same time as ULT initiation"
		lab define `med'_prophylaxis_`t'm 0 "No" 1 "Yes"
		lab val `med'_prophylaxis_`t'm `med'_prophylaxis_`t'm
			
		****Number of prescriptions issued in 6/12m after first script (Nb. doses will be double counted if different strengths issued)
		lab var `med'_count_`t'm "Number of prescriptions for `druglabel' within `t' months of diagnosis"
		tabstat `med'_count_`t'm, stats (n mean sd p50 p25 p75)
		
		****Number of prescriptions issued in 6/12m of first prescription (to assess continued prescribing) - Nb. only looks at prescriptions occurring after diagnosis
		reshape long `med'_date_, i(patient_id) j(`med'_order)
		gen `med'_within_`t'm = 1 if ((`med'_date_>=`med'_first_date) & (`med'_date_<=(`med'_first_date + `days'))) & `med'_date_!=.
		bys patient_id (`med'_within_`t'm): gen n=_n if `med'_within_`t'm!=.
		by patient_id: egen `med'_scripts_`t'm = max(n)
		drop n `med'_within_`t'm
		
		****Was the patient still prescribed drug at 6/12 months (+/- 2 month window) after start date
		gen `med'_ongoing_`t'm = 1 if ((`med'_date_>=(`med'_first_date + (`days' - 60))) & (`med'_date_<=(`med'_first_date + (`days' + 60)))) & `med'_date_!=.	
		sort patient_id `med'_ongoing_`t'm
		by patient_id: replace `med'_ongoing_`t'm = `med'_ongoing_`t'm[_n-1] if missing(`med'_ongoing_`t'm)
		
		reshape wide `med'_date_, i(patient_id) j(`med'_order)
		lab var `med'_scripts_`t'm "Number of prescriptions for `druglabel' within `t' months of first prescription"
		tabstat `med'_scripts_`t'm, stats (n mean sd p50 p25 p75)
		recode `med'_ongoing_`t'm .=0 if `med'_first_date!=.
		lab var `med'_ongoing_`t'm "Still prescribed `druglabel' `t' months after initiation"
		lab def `med'_ongoing_`t'm 0 "No" 1 "Yes", modify
		lab val `med'_ongoing_`t'm `med'_ongoing_`t'm
		tab `med'_ongoing_`t'm, missing
	}
}

***Generate categorical variables for drug initiation for box plots (amend list as necessary)
foreach med in `drug' {

	if "`med'" == "`drug'" {
		local druglabel = "`Drug'" 
	}
	else {
		local druglabel = strproper("`med'")
	}
	
	foreach t in 6 12 {
		local days = int((`t'/12)*365.25)
		di `days'
		local half_days = `days'/2
		local half = `t'/2
	 
		****Categorical variable, overall
		gen `med'_cat_`t'm=1 if time_to_`med'<=`half_days' & time_to_`med'!=. 
		replace `med'_cat_`t'm=2 if time_to_`med'>`half_days' & time_to_`med'<=`days' & time_to_`med'!=.
		replace `med'_cat_`t'm=3 if time_to_`med'>`days' | time_to_`med'==.
		lab define `med'_cat_`t'm 1 "Within `half' months" 2 "`half'-`t' months" 3 "No prescription within `t' months", modify
		lab val `med'_cat_`t'm `med'_cat_`t'm
		lab var `med'_cat_`t'm "`druglabel' prescribed within `t' months of diagnosis"
		tab `med'_cat_`t'm, missing
		
		****Categorical variable by year (for box plots)
		forvalues i = 1/$max_year {
			local start = $base_year + `i' - 1
			di "`start'"
			local end = `start' + 1
			gen `med'_cat_`t'm_`start'=`med'_cat_`t'm if diagnosis_year==`i'
			replace `med'_cat_`t'm_`start' = 9 if `med'_cat_`t'm_`start'==.
			lab define `med'_cat_`t'm_`start' 1 "Within 3 months" 2 "`half'-`t' months" 3 "No prescription within `t' months" 9 "Not diagnosed in this year", modify
			lab val `med'_cat_`t'm_`start' `med'_cat_`t'm_`start'
			lab var `med'_cat_`t'm_`start' "`druglabel' initiation within `t' months of diagnosis, July `start'-June `end'"
			tab `med'_cat_`t'm_`start', missing
		}	
	}
}

***For those prescribed febuxostat, did they have pre-existing CHD or CVA?
gen febux_mace = 0 if febuxostat_first_date!=.
replace febux_mace = 1 if febuxostat_first_date!=. & ((chd_date!=. & (chd_date <= febuxostat_first_date)) | (cva_date!=. & (cva_date <= febuxostat_first_date)))
lab var febux_mace "Febuxostat prescribed in those with MACE"
lab define febux_mace 0 "No" 1 "Yes"
lab val febux_mace febux_mace

save "$projectdir/output/data/cohort_meds.dta", replace

*Relevant blood tests (passed from yaml, but amend thresholds as necessary) ===================================================*/

use "$projectdir/output/data/cohort_meds.dta", clear

*local blood "urate"
foreach blood in $bloods {
	di "`blood'"
	
	**Set thresholds for plausible low and high values (amend as necessary)
	if "`blood'" == "creatinine" {
		local low 20
		local high 3000
	}
	else if "`blood'" == "hba1c" {
		local low 10
		local high 200
	}
	else if "`blood'" == "urate" {
		local low 0.05
		local high 3000
	}
	else if "`blood'" == "cholesterol" {
		local low 0.5
		local high 20
	}
	else {
		local low 0
		local high 100000
	}
	di "`low'"
	di "`high'"

	***Find number of tests
	local max = 0
	capture quietly ds `blood'_value_*
	if !_rc & "`r(varlist)'" != "" {
		foreach v of varlist `r(varlist)' {
			if regexm("`v'","^`blood'_value_([0-9]+)$") {
				local idx = real(regexs(1))
				local max = max(`max', `idx')
			}
		}
	}
	di "`max'"

	***Set implausible values to missing (as defined above)
	if `max'!=0 {
		forval i = 1/`max' 	{
			codebook `blood'_value_`i' //check
			tabstat `blood'_value_`i', stats(n mean sd p50 p25 p75)
			summ `blood'_value_`i' if inrange(`blood'_value_`i', `low', `high')
			replace `blood'_value_`i' = . if !inrange(`blood'_value_`i', `low', `high')
			replace `blood'_value_`i' = . if `blood'_date_`i' == . 
			replace `blood'_date_`i' = . if `blood'_value_`i' == . 
			codebook `blood'_value_`i' //check
			tabstat `blood'_value_`i', stats(n mean sd p50 p25 p75)
		}
	}

	/*
	***Nb for urate, depending on code, urate is in range 0.05 - 2 mmol/L or 50 - 2000 micromol/L; need also to consider mg/dL
			summ urate_value_`i' if inrange(urate_value_`i', 0.05, 2) // for mmol/L
			summ urate_value_`i' if inrange(urate_value_`i', 50, 2000) // for micromol/L
			summ urate_value_`i' if urate_value_`i'==. | urate_value_`i'==0 // missing or zero
			summ urate_value_`i' if ((!inrange(urate_value_`i', 0.05, 2)) & (!inrange(urate_value_`i', 50, 2000)) & urate_value_`i'!=. & urate_value_`i'!=0) // not missing or zero or in above ranges	
			replace urate_value_`i' = . if ((!inrange(urate_value_`i', 0.05, 2)) & (!inrange(urate_value_`i', 50, 2000))) // keep values that are consistent with mmol/L or micromol/L
			replace urate_value_`i' = (urate_value_`i'*1000) if inrange(urate_value_`i', 0.05, 2) //x1000 for values that are consistent with mmol/L or micromol/L
	*/

	***Reshape to long format
	reshape long `blood'_value_ `blood'_date_, i(patient_id) j(`blood'_order)

	***Keep only single blood test from same day, priotising ones not missing
	bys patient_id `blood'_date_ (`blood'_value_): gen n=_n
	drop if n>1 
	drop n

	***For individuals with at least one blood test, remove subsequent missing values
	bys patient_id (`blood'_date_): gen n=_n
	drop if n>1 & `blood'_value_==.
	drop n

	***Time from diagnosis date to blood test
	gen time_to_`blood' = `blood'_date_ - ${disease}_inc_date if `blood'_date_!=. 
	
	***Check date of blood test with respect to first ULT prescription (amend as necessary)
	gen `blood'_after_ult=1 if `blood'_date_>ult_first_date & `blood'_date_!=. & ult_first_date!=.
	replace `blood'_after_ult=0 if `blood'_date_<=ult_first_date & `blood'_date_!=. & ult_first_date!=.

	***Code baseline blood tests as value closest to diagnosis date, within 2 months before or after diagnosis date, assuming it was not after starting ULT
	gen abs_time_to_`blood' = abs(time_to_`blood') if time_to_`blood'!=. & time_to_`blood'<=60 & time_to_`blood'>=-60 & `blood'_after_ult!=1 & `blood'_value_!=.
	bys patient_id (abs_time_to_`blood'): gen n=_n 
	gen `blood'_bl_value=`blood'_value_ if n==1 & abs_time_to_`blood'!=. 
	lab var `blood'_bl_value "Serum `blood' at baseline"
	gen had_`blood'_bl = 1 if `blood'_bl_value!=.
	lab var had_`blood'_bl "Serum `blood' performed at baseline"
	lab define had_`blood'_bl 0 "No" 1 "Yes", modify
	lab val had_`blood'_bl had_`blood'_bl
	gen `blood'_bl_date=`blood'_date_ if n==1 & abs_time_to_`blood'!=.
	format `blood'_bl_date %td 
	drop n
	by patient_id: replace `blood'_bl_value = `blood'_bl_value[_n-1] if missing(`blood'_bl_value)
	by patient_id: replace `blood'_bl_date = `blood'_bl_date[_n-1] if missing(`blood'_bl_date)
	by patient_id: replace had_`blood'_bl = had_`blood'_bl[_n-1] if missing(had_`blood'_bl)
	recode had_`blood'_bl .=0

	***Code blood tests within 6/12 months after diagnosis
	foreach t in 6 12 {
		local days = int((`t'/12)*365.25)
		di `days'
		gen `blood'_within_`t'm = 1 if time_to_`blood' >=0 & time_to_`blood' <=`days'
		local Blood = upper(substr("`blood'",1,1)) + substr("`blood'",2,.)
		lab var `blood'_within_`t'm "`Blood'"
		lab define `blood'_within_`t'm 0 "No" 1 "Yes"
		lab val `blood'_within_`t'm `blood'_within_`t'm
		sort patient_id `blood'_within_`t'm
		by patient_id: replace `blood'_within_`t'm = `blood'_within_`t'm[_n-1] if missing(`blood'_within_`t'm)
		recode `blood'_within_`t'm .=0
	}
	
	***Reshape to wide format
	drop `blood'_after_ult time_to_`blood' abs_time_to_`blood'
	reshape wide `blood'_value_ `blood'_date_, i(patient_id) j(`blood'_order)
	order `blood'_bl_value, after(`blood'_bl_date) 
}

**Generate eGFR from serum creatinine (using CKD-EPI formula with no ethnicity) ============================*/ 
reshape long creatinine_value_ creatinine_date_, i(patient_id) j(creatinine_order)
codebook creatinine_value if creatinine_order==1 //check
codebook creatinine_date if creatinine_order==1 //check
 
gen SCr_adj = creatinine_value/88.4

gen min = .
replace min = SCr_adj/0.7 if sex==1
replace min = SCr_adj/0.9 if sex==2
replace min = min^-0.329  if sex==1
replace min = min^-0.411  if sex==2
replace min = 1 if min<1

gen max = .
replace max=SCr_adj/0.7 if sex==1
replace max=SCr_adj/0.9 if sex==2
replace max=max^-1.209
replace max=1 if max>1

gen egfr_value_=min*max*141
replace egfr_value_=egfr_value_*(0.993^age)
replace egfr_value_=egfr_value_*1.018 if sex==1
label var egfr_value_ "eGFR"
drop min max SCr_adj

gen egfr_date_ = creatinine_date_
format egfr_date_ %td
codebook egfr_value_ if creatinine_order==1 //check

***Date of first eGFR <60 (from 24 months before diagnosis to end of study follow-up) - consider amending 24-month limit for survival analyses
gen egfr_ckd = 1 if egfr_value_< 60 & !missing(egfr_value_)
bys patient_id egfr_ckd (egfr_date_): gen n=_n if egfr_ckd!=.
gen first_egfr_ckd_date = egfr_date_ if n==1
format first_egfr_ckd_date %td
drop n
sort patient_id first_egfr_ckd_date
by patient_id: replace first_egfr_ckd_date = first_egfr_ckd_date[_n-1] if missing(first_egfr_ckd_date)
lab var first_egfr_ckd_date "Date of first eGFR <60"

***Date of second eGFR <60 at least 90 days after first eGFR <60
gen subsequent_egfr_ckd = 1 if egfr_ckd == 1 & !missing(first_egfr_ckd_date) & (egfr_date_ >= (first_egfr_ckd_date + 90))
bys patient_id subsequent_egfr_ckd (egfr_date_): gen n=_n if subsequent_egfr_ckd==1
gen second_egfr_ckd_date = egfr_date_ if n==1
format second_egfr_ckd_date %td
drop n
sort patient_id second_egfr_ckd_date
by patient_id: replace second_egfr_ckd_date = second_egfr_ckd_date[_n-1] if missing(second_egfr_ckd_date)
lab var second_egfr_ckd_date "Date of second eGFR <60, at least 90 days after first"

***Baseline eGFR value with respect to diagnosis date
gen egfr_bl_date = egfr_date_ if creatinine_bl_date == egfr_date_ & creatinine_bl_date!=. & egfr_value_!=. & egfr_date_!=.
format egfr_bl_date %td
gen egfr_bl_value = egfr_value_ if creatinine_bl_date == egfr_date_ & creatinine_bl_date!=. & egfr_value_!=. & egfr_date_!=.
sort patient_id egfr_bl_date
by patient_id: replace egfr_bl_date = egfr_bl_date[_n-1] if missing(egfr_bl_date)
sort patient_id egfr_bl_value
by patient_id: replace egfr_bl_value = egfr_bl_value[_n-1] if missing(egfr_bl_value)
lab var egfr_bl_value "Baseline eGFR at diagnosis"
lab var egfr_bl_date "Date of eGFR at diagnosis"

***Baseline eGFR value with respect to ULT initiation date (up to 12 months before)
gen time_egfr_before_ult = (ult_first_date - egfr_date_) if egfr_date_!=. & ult_first_date!=.
gen egfr_before_ult = 1 if (time_egfr_before_ult<=365) & (time_egfr_before_ult>0) & time_egfr_before_ult!=. & egfr_value_!=. //blood tests on or before initiating ULT
bys patient_id egfr_before_ult (time_egfr_before_ult): gen n=_n if egfr_before_ult==1
gen egfr_before_ult_value = egfr_value_ if n==1
lab var egfr_before_ult_value "Baseline eGFR before initiating ULT"
sort patient_id egfr_before_ult_value
by patient_id: replace egfr_before_ult_value = egfr_before_ult_value[_n-1] if missing(egfr_before_ult_value)
tabstat egfr_before_ult_value, stats(n mean sd p50 p25 p75)
drop n egfr_before_ult time_egfr_before_ult

drop egfr_ckd subsequent_egfr_ckd

reshape wide creatinine_value_ creatinine_date_ egfr_value_ egfr_date_, i(patient_id) j(creatinine_order)

***Categorise baseline eGFR into CKD stages
gen egfr_bl_cat = .
recode egfr_bl_cat . = 3 if egfr_bl_value < 30
recode egfr_bl_cat . = 2 if egfr_bl_value < 60
recode egfr_bl_cat . = 1 if egfr_bl_value < .
replace egfr_bl_cat = 9 if egfr_bl_value >= .

label define egfr_bl_cat 	1 ">=60" 		///
							2 "30-59"		///
							3 "<30"			///
							9 "Not known"
					
label val egfr_bl_cat egfr_bl_cat
lab var egfr_bl_cat "eGFR at baseline"
tab egfr_bl_cat, missing

***Categorise baseline eGFR into more granular CKD stages
gen egfr_bl_finecat = .
recode egfr_bl_finecat . = 6 if egfr_bl_value < 15
recode egfr_bl_finecat . = 5 if egfr_bl_value < 30
recode egfr_bl_finecat . = 4 if egfr_bl_value < 45
recode egfr_bl_finecat . = 3 if egfr_bl_value < 60
recode egfr_bl_finecat . = 2 if egfr_bl_value < 90
recode egfr_bl_finecat . = 1 if egfr_bl_value < .
replace egfr_bl_finecat = 9 if egfr_bl_value >= .

label define egfr_bl_finecat 	1 ">=90" 		///
								2 "60-89"		///
								3 "45-59"		///
								4 "30-44"		///
								5 "15-29"		///
								6 "<15"			///
								9 "Not known"
					
label val egfr_bl_finecat egfr_bl_finecat
lab var egfr_bl_finecat "eGFR at baseline"
tab egfr_bl_finecat, missing

***Generate baseline CKD code that combines CKD coding (stages 3-5) + eGFR (stages 3-5)
gen ckd_comb_bl = 0
replace ckd_comb_bl = 1 if egfr_bl_value != . & egfr_bl_value < 60
replace ckd_comb_bl = 1 if ckd_bl == 1
label define ckd_comb_bl 0 "No" 1 "Yes"
label val ckd_comb_bl ckd_comb_bl
label var ckd_comb_bl "CKD"
tab ckd_comb_bl, missing

***Generate baseline CKD code that combines CKD coding (stages 3-5) + eGFR (stages 3b-5) - bespoke, for referral metric
gen ckd_comb_bl_3b = 0
replace ckd_comb_bl_3b = 1 if egfr_bl_value != . & egfr_bl_value < 45
replace ckd_comb_bl_3b = 1 if ckd_bl == 1
label define ckd_comb_bl_3b 0 "No" 1 "Yes"
label val ckd_comb_bl_3b ckd_comb_bl_3b
label var ckd_comb_bl_3b "Chronic kidney disease (stages 3b-5)"
tab ckd_comb_bl_3b, missing

**Generate date of first CKD (code or eGFR), including those after diagnosis
gen ckd_comb = 0
replace ckd_comb = 1 if first_egfr_ckd_date!=. //any CKD blood test from 2 years before diagnosis to study end date
replace ckd_comb = 1 if ckd_bl == 1 | ckd_new == 1 //any CKD code ever
label define ckd_comb 0 "No" 1 "Yes"
label val ckd_comb ckd_comb
label var ckd_comb "CKD"
tab ckd_comb, missing

gen temp1_ckd_date = first_egfr_ckd_date if first_egfr_ckd_date!=.
gen temp2_ckd_date = ckd_date if ckd_date!=.
gen first_ckd_comb_date = min(temp1_ckd_date,temp2_ckd_date) 
format first_ckd_comb_date %td 
drop temp1_ckd_date temp2_ckd_date

gen ckd_comb_12m = 0
replace ckd_comb_12m = 1 if first_ckd_comb_date!=. & (first_ckd_comb_date <= (${disease}_inc_date + 365))
label define ckd_comb_12m 0 "No" 1 "Yes"
label val ckd_comb_12m ckd_comb_12m
label var ckd_comb_12m "CKD"
tab ckd_comb_12m, missing

**Categorise HbA1c at baseline ============================*/
codebook hba1c_bl_value //check reasonable (i.e. mmol, not %) - would be screened out by the above
codebook hba1c_bl_date //check

gen hba1c_bl_cat = 0 if hba1c_bl_value < 58
replace hba1c_bl_cat = 1 if hba1c_bl_value >= 58 & hba1c_bl_val !=.
replace hba1c_bl_cat = 9 if hba1c_bl_cat ==. 
label define hba1c_bl_cat 0 "HbA1c <58mmol/mol" 1 "HbA1c >=58mmol/mol" 9 "Not known"
label val hba1c_bl_cat hba1c_bl_cat
lab var hba1c_bl_cat "HbA1c at baseline"
tab hba1c_bl_cat, missing

***Create combined diabetes code that combines diabetes coding + HBA1c
gen diab_bl_cat = 1 if diabetes_bl==0
replace diab_bl_cat = 2 if diabetes_bl==1 & hba1c_bl_cat==0
replace diab_bl_cat = 3 if diabetes_bl==1 & hba1c_bl_cat==1
replace diab_bl_cat = 4 if diabetes_bl==1 & hba1c_bl_cat==9

label define diab_bl_cat 1 "No diabetes" 			///
						2 "Diabetes with HbA1c <58mmol/mol"		///
						3 "Diabetes with HbA1c >58mmol/mol" 	///
						4 "Diabetes with no recorded HbA1c"
label values diab_bl_cat diab_bl_cat
lab var diab_bl_cat "Type 2 diabetes mellitus with HbA1c categorisation"
tab diab_bl_cat, missing

**For urate levels, perform further cleaning (remove if not necessary) ================================*/

**Categorise baseline urate levels
gen urate_bl_cat=0 if urate_bl_value>=360 & urate_bl_value!=.
replace urate_bl_cat=1 if urate_bl_value<360 & urate_bl_value!=.
replace urate_bl_cat=9 if urate_bl_value==.
lab var urate_bl_cat "Baseline serum urate level"
lab define urate_bl_cat 0 ">=360 micromol/L" 1 "<360 micromol/L" 9 "Not known", modify
lab val urate_bl_cat urate_bl_cat
tab urate_bl_cat, missing

**If baseline urate level <360, check whether repeat test was performed within 1 month of diagnosis (or within 1 month of test if this was up to 1 month after diagnosis)
gen urate_bl_360_repeat = 0 if urate_bl_cat==2

***Find max number of urate levels
local max = 0
capture quietly ds urate_date_*
if !_rc & "`r(varlist)'" != "" {
	foreach v of varlist `r(varlist)' {
		if regexm("`v'","^urate_date_([0-9]+)$") {
			local idx = real(regexs(1))
			local max = max(`max', `idx')
		}
	}
}
di "`max'"

***Find matching urate levels within 3 months after baseline urate (if baseline <360)
forval i = 1/`max'	{
	replace urate_bl_360_repeat = 1 if urate_bl_cat==2 & ((urate_date_`i' > urate_bl_date) & (urate_date_`i' <= (urate_bl_date + 90)) & urate_date_`i'!=.)
}
lab var urate_bl_360_repeat "Repeat urate level within three months if baseline urate <360"
lab def urate_bl_360_repeat 0 "No" 1 "Yes"
lab val urate_bl_360_repeat urate_bl_360_repeat

***Reshape to long format
reshape long urate_value_ urate_date_, i(patient_id) j(urate_order)

***Define first serum urate <360/300 micromol/L after diagnosis 
gen time_to_urate = (urate_date_ - ${disease}_inc_date) if urate_date_!=. //time to urate level
gen urate_after_diag = 1 if time_to_urate>0 & time_to_urate!=. & urate_value_!=.

foreach target in 300 360 {
	gen urate_below`target' = 1 if urate_after_diag==1 & urate_value_<`target'
	bys patient_id urate_below`target' (urate_date_): gen n=_n
	gen urate_below`target'_date = urate_date_ if n==1 & urate_below`target'==1
	format urate_below`target'_date %td
	lab var urate_below`target'_date "Date of first serum urate <`target' micromol/L after diagnosis"
	sort patient_id urate_below`target'_date
	by patient_id: replace urate_below`target'_date = urate_below`target'_date[_n-1] if missing(urate_below`target'_date)
	gen urate_below`target'_value = urate_value_ if n==1 & urate_below`target'==1
	lab var urate_below`target'_value "Value of first serum urate <`target' micromol/L after diagnosis"
	sort patient_id urate_below`target'_value
	by patient_id: replace urate_below`target'_value = urate_below`target'_value[_n-1] if missing(urate_below`target'_value)
	tabstat urate_below`target'_value, stats(n mean sd p50 p25 p75)
	drop n urate_below`target'
}

sort patient_id urate_after_diag
by patient_id: replace urate_after_diag = urate_after_diag[_n-1] if missing(urate_after_diag)
recode urate_after_diag .=0
lab var urate_after_diag "Urate test performed after diagnosis"
lab def urate_after_diag 0 "No" 1 "Yes"
lab val urate_after_diag urate_after_diag
tab urate_after_diag, missing

***Define proportion of patients who attained serum urate <360/300 micromol/L within 6/12 months of diagnosis, irrespective of ULT use - can subsequently limit this to those who had 6m/12m of follow-up post-diagnosis
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'

	gen urate_`t'm = 1 if (time_to_urate>0 & time_to_urate<=`days') & urate_value_!=. //test done within 6/12 months of diagnosis
	bys patient_id (urate_`t'm): gen n=_n if urate_`t'm!=.
	by patient_id: egen urate_count_`t'm = max(n) //number of tests within 6/12 months
	recode urate_count_`t'm .=0
	lab var urate_count_`t'm "Number of urate levels within `t' months of diagnosis"
	gen two_urate_`t'm = 1 if urate_count_`t'm >=2 & urate_count_`t'm!=. //two or more urate tests performed within 6/12 months of diagnosis
	recode two_urate_`t'm .=0 //includes those who didn't receive ULT
	lab var two_urate_`t'm "Two or more serum urate tests performed within `t' months of diagnosis"
	lab def two_urate_`t'm 0 "No" 1 "Yes", modify
	lab val two_urate_`t'm two_urate_`t'm
	drop n
	sort patient_id urate_`t'm
	by patient_id: replace urate_`t'm = urate_`t'm[_n-1] if missing(urate_`t'm)
	recode urate_`t'm .=0
	lab var urate_`t'm "Serum urate test performed within `t' months of diagnosis"
	lab def urate_`t'm 0 "No" 1 "Yes", modify
	lab val urate_`t'm urate_`t'm
	gen urate_value_`t'm = urate_value_ if (time_to_urate>0 & time_to_urate<=`days') & urate_value_!=. //test values within 6/12 months of diagnosis
	bys patient_id (urate_value_`t'm): gen n=_n if urate_value_`t'm!=.
	gen lowest_urate_`t'm = urate_value_`t'm if n==1 //lowest urate value within 6/12 months of diagnosis
	lab var lowest_urate_`t'm "Lowest serum urate level within `t' months of diagnosis"
	sort patient_id (lowest_urate_`t'm)
	by patient_id: replace lowest_urate_`t'm = lowest_urate_`t'm[_n-1] if missing(lowest_urate_`t'm)
	drop n urate_value_`t'm
	
	foreach target in 300 360 { 
		***Categorical variable (coded missing)
		gen urate_`target'_`t'm_cat = 1 if lowest_urate_`t'm<`target' & lowest_urate_`t'm!=.
		replace urate_`target'_`t'm_cat = 0 if lowest_urate_`t'm>=`target' & lowest_urate_`t'm!=.
		replace urate_`target'_`t'm_cat = 9 if lowest_urate_`t'm==.
		lab var urate_`target'_`t'm_cat  "Serum urate <`target' micromol/L within `t' months of diagnosis"
		lab def urate_`target'_`t'm_cat 0 ">=`target' micromol/L" 1 "<`target' micromol/L" 9 "Not known", modify
		lab val urate_`target'_`t'm_cat urate_`target'_`t'm_cat
		tab urate_`target'_`t'm_cat, missing
		
		***Binary variable (uncoded missing)
		gen urate_`target'_`t'm = 1 if lowest_urate_`t'm<`target' & lowest_urate_`t'm!=.
		replace urate_`target'_`t'm = 0 if lowest_urate_`t'm>=`target' & lowest_urate_`t'm!=.
		lab var urate_`target'_`t'm "Serum urate <`target' micromol/L within `t' months of diagnosis"
		lab def urate_`target'_`t'm 0 ">=`target' micromol/L" 1 "<`target' micromol/L", modify
		lab val urate_`target'_`t'm urate_`target'_`t'm
		tab urate_`target'_`t'm, missing
	}
	
	tab urate_`t'm, missing
	tabstat urate_count_`t'm, stats(n mean sd p50 p25 p75)
	tab two_urate_`t'm, missing
	tabstat lowest_urate_`t'm, stats(n mean sd p50 p25 p75)
}
drop time_to_urate

***Define first serum urate <360/300 micromol/L after initiating ULT 
gen time_ult_to_urate = (urate_date_ - ult_first_date) if urate_date_!=. & ult_first_date!=.
gen urate_after_ult = 1 if time_ult_to_urate>0 & time_ult_to_urate!=. & urate_value_!=. //blood tests after initiating ULT

foreach target in 300 360 { 
	gen urate_below`target'_ult = 1 if urate_after_ult==1 & urate_value_<`target'
	bys patient_id urate_below`target'_ult (urate_date_): gen n=_n
	gen urate_below`target'_ult_date = urate_date_ if n==1 & urate_below`target'_ult==1
	format urate_below`target'_ult_date %td
	lab var urate_below`target'_ult_date "Date of first serum urate <`target' micromol/L after initiating ULT"
	sort patient_id urate_below`target'_ult_date 
	by patient_id: replace urate_below`target'_ult_date = urate_below`target'_ult_date[_n-1] if missing(urate_below`target'_ult_date)
	gen urate_below`target'_ult_value = urate_value_ if n==1 & urate_below`target'_ult==1
	lab var urate_below`target'_ult_value "Value of first serum urate <`target' micromol/L after initiating ULT"
	sort patient_id urate_below`target'_ult_value
	by patient_id: replace urate_below`target'_ult_value = urate_below`target'_ult_value[_n-1] if missing(urate_below`target'_ult_value)
	tabstat urate_below`target'_ult_value, stats(n mean sd p50 p25 p75)
	drop n urate_below`target'_ult
}	

sort patient_id urate_after_ult
by patient_id: replace urate_after_ult = urate_after_ult[_n-1] if missing(urate_after_ult)
recode urate_after_ult .=0
lab var urate_after_ult "Urate test performed after initiating ULT"
lab def urate_after_ult 0 "No" 1 "Yes"
lab val urate_after_ult urate_after_ult
tab urate_after_ult, missing

foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	***For those who attained serum urate target (at any point), was urate repeated after 12 months (between 6 and 18 months afterwards) - can then limit this to a timeframe after ULT initiation
	gen repeat_after360_`t'm_ult = 1 if urate_below360_ult_date!=. & urate_date_!=. & (urate_date_ > (urate_below360_ult_date + 183)) & (urate_date_ <= (urate_below360_ult_date + 548))
	gen repeat_below360_`t'm_ult = 1 if repeat_after360_`t'm_ult == 1 & urate_value_<360
	sort patient_id repeat_after360_`t'm_ult
	by patient_id: replace repeat_after360_`t'm_ult = repeat_after360_`t'm_ult[_n-1] if missing(repeat_after360_`t'm_ult)
	lab var repeat_after360_`t'm_ult "Repeat serum urate test performed after achieving target <360 micromol/L following ULT initiation"
	recode repeat_after360_`t'm_ult .=0 if urate_below360_ult_date!=.
	lab def repeat_after360_`t'm_ult 0 "No" 1 "Yes"
	lab val repeat_after360_`t'm_ult repeat_after360_`t'm_ult
	tab repeat_after360_`t'm_ult, missing 
	sort patient_id repeat_below360_`t'm_ult
	by patient_id: replace repeat_below360_`t'm_ult = repeat_below360_`t'm_ult[_n-1] if missing(repeat_below360_`t'm_ult)
	lab var repeat_below360_`t'm_ult "Repeat serum urate level remains <360 micromol/L following ULT initiation" 
	recode repeat_below360_`t'm_ult .=0 if repeat_after360_`t'm_ult!=.
	lab def repeat_below360_`t'm_ult 0 "No" 1 "Yes"
	lab val repeat_below360_`t'm_ult repeat_below360_`t'm_ult
	tab repeat_below360_`t'm_ult, missing 

	***Define proportion of patients who attained serum urate <360 micromol/L within 6/12 months of initiating ULT - can subsequently limit this to those who initiated ULT within 6/12m of diagnosis and those with 6/12m of follow-up post-ULT
	gen urate_within_`t'm_ult = 1 if (time_ult_to_urate>0 & time_ult_to_urate<=`days') & urate_value_!=. //test done within 6/12 months of ULT initiation
	bys patient_id (urate_within_`t'm_ult): gen n=_n if urate_within_`t'm_ult!=.
	by patient_id: egen urate_count_`t'm_ult = max(n) //number of tests within 6/12 months of ULT initiation
	recode urate_count_`t'm_ult .=0
	lab var urate_count_`t'm_ult "Number of serum urate tests within `t' months of ULT initiation"
	gen two_urate_`t'm_ult = 1 if urate_count_`t'm_ult >=2 & urate_count_`t'm_ult!=. //two or more urate tests performed within 6/12 months of ULT initiation
	recode two_urate_`t'm_ult .=0 //includes those who didn't receive ULT
	lab var two_urate_`t'm_ult "Two or more serum urate tests performed within `t' months of ULT initiation"
	lab def two_urate_`t'm_ult 0 "No" 1 "Yes", modify
	lab val two_urate_`t'm_ult two_urate_`t'm_ult
	drop n
	sort patient_id urate_within_`t'm_ult
	by patient_id: replace urate_within_`t'm_ult= urate_within_`t'm_ult[_n-1] if missing(urate_within_`t'm_ult)
	recode urate_within_`t'm_ult .=0 //includes those who didn't receive ULT
	lab var urate_within_`t'm_ult "Serum urate test performed within `t' months of ULT initiation"
	lab def urate_within_`t'm_ult 0 "No" 1 "Yes", modify
	lab val urate_within_`t'm_ult urate_within_`t'm_ult
	gen urate_value_`t'm_ult = urate_value_ if (time_ult_to_urate>0 & time_ult_to_urate<=`days') & urate_value_!=. //test values within 6/12 months of ULT initiation
	bys patient_id (urate_value_`t'm_ult): gen n=_n if urate_value_`t'm_ult!=.
	gen lowest_urate_`t'm_ult = urate_value_`t'm_ult if n==1 //lowest urate value within 6/12 months of ULT initiation
	lab var lowest_urate_`t'm_ult "Lowest urate level within `t' months of ULT initiation"
	sort patient_id (lowest_urate_`t'm_ult)
	by patient_id: replace lowest_urate_`t'm_ult = lowest_urate_`t'm_ult[_n-1] if missing(lowest_urate_`t'm_ult)
	drop n urate_value_`t'm_ult
	
	***Binary variable, overall (uncoded missing)
	gen urate_`t'm_ult = 1 if lowest_urate_`t'm_ult<360 & lowest_urate_`t'm_ult!=.
	replace urate_`t'm_ult = 0 if lowest_urate_`t'm_ult>=360 & lowest_urate_`t'm_ult!=.
	lab var urate_`t'm_ult  "Serum urate <`target' micromol/L within `t' months of ULT initiation"
	lab def urate_`t'm_ult 0 ">=360 micromol/L" 1 "<360 micromol/L", modify
	lab val urate_`t'm_ult urate_`t'm_ult
	
	***Categorical variable, overall (coded missing)
	gen urate_`t'm_ult_cat = 1 if lowest_urate_`t'm_ult<360 & lowest_urate_`t'm_ult!=.
	replace urate_`t'm_ult_cat = 0 if lowest_urate_`t'm_ult>=360 & lowest_urate_`t'm_ult!=.
	replace urate_`t'm_ult_cat = 9 if lowest_urate_`t'm_ult==. //includes those who didn't receive ULT
	lab var urate_`t'm_ult_cat  "Serum urate <`target' micromol/L within `t' months of ULT initiation"
	lab def urate_`t'm_ult_cat 0 ">=360 micromol/L" 1 "<360 micromol/L" 9 "Not known", modify
	lab val urate_`t'm_ult_cat urate_`t'm_ult_cat
	
	***Binary variable, overall (missing recoded as not attained)
	gen urate_`t'm_ult_recode = urate_`t'm_ult
	replace urate_`t'm_ult_recode = 0 if urate_`t'm_ult ==. //includes those who didn't receive ULT
	lab var urate_`t'm_ult_recode  "Serum urate <`target' micromol/L within `t' months of ULT initiation"
	lab def urate_`t'm_ult_recode 0 ">=360 micromol/L or not known" 1 "<360 micromol/L", modify
	lab val urate_`t'm_ult_recode urate_`t'm_ult_recode
		
	tab urate_within_`t'm_ult, missing
	tabstat urate_count_`t'm_ult, stats(n mean sd p50 p25 p75)
	tab two_urate_`t'm_ult, missing
	tab urate_`t'm_ult_cat, missing
	tabstat lowest_urate_`t'm_ult, stats(n mean sd p50 p25 p75)
	
	***Categorical variable by year of ULT initiation (for box plots)
	forvalues i = 1/$max_year {
		local start = $base_year + `i' - 1
		di "`start'"
		local end = `start' + 1
		gen urate_`t'm_ult_cat_`start'=urate_`t'm_ult_cat if ult_year==`i'
		replace urate_`t'm_ult_cat_`start' = 10 if urate_`t'm_ult_cat_`start'==.
		lab define urate_`t'm_ult_cat_`start' 0 ">=360 micromol/L" 1 "<360 micromol/L" 9 "Not known" 10 "Initial ULT prescription not in this year", modify
		lab val urate_`t'm_ult_cat_`start' urate_`t'm_ult_cat_`start'
		lab var urate_`t'm_ult_cat_`start' "Serum urate <`target' micromol/L within `t' months of ULT initiation, July `start'-June `end'"
		tab urate_`t'm_ult_cat_`start', missing
	}
}

drop time_ult_to_urate

**Baseline urate level in the 12 months on or before ULT initiation
gen time_urate_before_ult = (ult_first_date - urate_date_) if urate_date_!=. & ult_first_date!=.
gen urate_before_ult = 1 if (time_urate_before_ult<=365) & (time_urate_before_ult>0) & time_urate_before_ult!=. & urate_value_!=. //blood tests on or before initiating ULT
bys patient_id urate_before_ult (time_urate_before_ult): gen n=_n if urate_before_ult==1
gen urate_before_ult_value = urate_value_ if n==1
lab var urate_before_ult_value "Baseline serum urate level before initiating ULT"
sort patient_id urate_before_ult_value
by patient_id: replace urate_before_ult_value = urate_before_ult_value[_n-1] if missing(urate_before_ult_value)
tabstat urate_before_ult_value, stats(n mean sd p50 p25 p75)
drop n time_urate_before_ult

reshape wide urate_value_ urate_date_, i(patient_id) j(urate_order)

save "$projectdir/output/data/cohort_bloods.dta", replace

*For clinical events (flares, admissions, ED attendances), perform bespoke cleaning (amend as necessary) ==============================================================================*/

use "$projectdir/output/data/cohort_bloods.dta", clear

**Criteria for defining flares adapted from https://jamanetwork.com/journals/jama/fullarticle/2794763: 1) presence of a non-index diagnostic code for gout flare (specified by dedicated flare codelist); 2) non-index admission with primary gout diagnostic code; 3) non-index ED attendance with primary gout diagnostic code; 4) any non-index gout diagnostic code AND prescription for a flare treatment on same day as that code. Exclude events that occur within 14 days of one another (handled in dataset definition for 1, 2 and 3)

**Store list of admissions with primary gout diagnostic codes
preserve
reshape long gout_adm_date_, i(patient_id) j(admission_order)
gen flare_overall_date=gout_adm_date_
format %td flare_overall_date
keep patient_id flare_overall_date
save "$projectdir/output/data/adm_dates_long.dta", replace
restore

**Store list of ED attendances with primary gout diagnostic codes
preserve
reshape long gout_ed_date_, i(patient_id) j(emerg_order)
gen flare_overall_date=gout_ed_date_
format %td flare_overall_date
keep patient_id flare_overall_date
save "$projectdir/output/data/emerg_dates_long.dta", replace
restore

**Store list of gout flare code dates with primary gout diagnostic codes
preserve
reshape long flare_date_ , i(patient_id) j(flare_order)
gen flare_overall_date=flare_date_
format %td flare_overall_date
keep patient_id flare_overall_date
save "$projectdir/output/data/flare_dates_long.dta", replace
restore

**Store list of gout consultations with prescriptions for flare medications on the same date
reshape long gout_cons_date_ , i(patient_id) j(consult_order)
gen code_and_tx_date_=.
format %td code_and_tx_date_

***Find matching consult and treatment dates
forval i = 1/$max_prescription {
	replace code_and_tx_date_ = gout_cons_date_ if gout_cons_date_ == colchicine_date_`i' & colchicine_date_`i'!=. & gout_cons_date_!=. | gout_cons_date_ == nsaid_date_`i' & nsaid_date_`i'!=. & gout_cons_date_!=. | gout_cons_date_ == steroid_date_`i' & steroid_date_`i'!=. & gout_cons_date_!=.
}
replace code_and_tx_date_=. if (code_and_tx_date_ < (${disease}_inc_date+14)) & code_and_tx_date_!=. //remove events within 14 days after diagnosis 

***Remove events that occur within 14 days of one another (repeat this until no further events within 14 days of one another)
sort patient_id code_and_tx_date_
local changed 1
while `changed' {
    bys patient_id (code_and_tx_date_): gen n=_n
    quietly count if n>1 & ((code_and_tx_date_-14)<(code_and_tx_date_[_n-1])) & code_and_tx_date_!=. & code_and_tx_date_[_n-1]!=.
    local changed = r(N)
    replace code_and_tx_date_=. if n>1 & ((code_and_tx_date_-14)<(code_and_tx_date_[_n-1])) & code_and_tx_date_!=. & code_and_tx_date_[_n-1]!=.
    drop n
	***Check log to see how many are being removed in each iteration
}
gen flare_overall_date=code_and_tx_date_
format %td flare_overall_date
keep patient_id flare_overall_date
save "$projectdir/output/data/code_and_tx_dates_long.dta", replace

**Append admissions, ED attendances and flare codes
append using "$projectdir/output/data/emerg_dates_long.dta"
append using "$projectdir/output/data/adm_dates_long.dta"
append using "$projectdir/output/data/flare_dates_long.dta"

**Remove events that occur within 14 days of one another (repeat this until no further events within 14 days of one another)
sort patient_id flare_overall_date
local changed 1
while `changed' {
    bys patient_id (flare_overall_date): gen n=_n
    quietly count if n>1 & ((flare_overall_date-14)<(flare_overall_date[_n-1])) & flare_overall_date!=. & flare_overall_date[_n-1]!=.
    local changed = r(N)
    replace flare_overall_date=. if n>1 & ((flare_overall_date-14)<(flare_overall_date[_n-1])) & flare_overall_date!=. & flare_overall_date[_n-1]!=.
    drop n
	***Check log to see how many are being removed in each iteration
}

**Generate overall flare counts and first post-diagnosis flare date
bys patient_id (flare_overall_date): gen n=_n if flare_overall_date!=.
by patient_id: egen flare_overall_count=max(n) //count of flares/admissions/ED after diagnosis
lab var flare_overall_count "Number of gout flares after diagnosis"
recode flare_overall_count .=0
drop n
bys patient_id (flare_overall_date): gen n=_n
drop if n>1 & flare_overall_date==. //drop missing values after the first row
gen first_flare_overall_date = flare_overall_date if n==1
format first_flare_overall_date %td
lab var first_flare_overall_date "Date of first flare after diagnosis"
by patient_id: replace first_flare_overall_date= first_flare_overall_date[_n-1] if missing(first_flare_overall_date)
rename n flare_overall_order
rename flare_overall_date flare_overall_date_
save "$projectdir/output/data/flares.dta", replace

**Keep long format and merge original data to obtain flare treatment and blood data
use "$projectdir/output/data/flares.dta", clear
merge m:1 patient_id using "$projectdir/output/data/cohort_bloods.dta", keep(match) nogen

**Flare dates that received colchicine vs. NSAIDs vs. steroids on same day
gen flare_drug_ = 0 if flare_overall_date_!=.

forval i = 1/$max_prescription {
	replace flare_drug_ = 1 if flare_overall_date_!=. & flare_overall_date_ == colchicine_date_`i' & colchicine_date_`i'!=. 
	replace flare_drug_ = 2 if flare_overall_date_!=. & flare_overall_date_ == nsaid_date_`i' & nsaid_date_`i'!=. 
	replace flare_drug_ = 3 if flare_overall_date_!=. & flare_overall_date_ == steroid_date_`i' & steroid_date_`i'!=. 
}

lab var flare_drug_ "Drug used for treatment of gout flare"
lab define flare_drug_ 0 "No drug" 1 "Colchicine" 2 "NSAID" 3 "Corticosteroid"
lab val flare_drug_ flare_drug_

**Save a long format dta for downstream flare analyses
preserve
keep patient_id flare_overall_date_ flare_drug_ $demographic
save "$projectdir/output/data/flares_long.dta", replace
restore

**Revert to wide format
reshape wide flare_overall_date_ flare_drug_, i(patient_id) j(flare_overall_order)

**Check whether at least 12 months of follow-up after first flare
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	gen has_`t'm_fup_flare=1 if (reg_end_date!=. & (reg_end_date >= (first_flare_overall_date + `days')) & ((first_flare_overall_date + `days') <= (date("$studyfup_date", "YMD")))) | (reg_end_date==. & ((first_flare_overall_date + `days') <= (date("$studyfup_date", "YMD"))))
	recode has_`t'm_fup_flare .=0 //includes those who didn't have flare
	lab var has_`t'm_fup_flare "At least `t' months of follow-up after first flare"
	lab def has_`t'm_fup_flare 0 "No" 1 "Yes"
	lab val has_`t'm_fup_flare has_`t'm_fup_flare
	tab has_`t'm_fup_flare
}

**Categorise patients who have one of more flares at any time after diagnosis
tabstat flare_overall_count, stats (n mean p50 p25 p75)
gen any_flares = 1 if flare_overall_count>=1 & flare_overall_count!=.
recode any_flares .=0
lab var any_flares "At least one additional flare after diagnosis"
lab def any_flares 0 "No" 1 "Yes", modify
lab val any_flares any_flares
tab any_flares

**Categorise patients who have one of more flares within t months after diagnosis
foreach t in 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	gen any_flares_`t'm = 1 if first_flare_overall_date!=. & (first_flare_overall_date < (${disease}_inc_date + `days'))
	recode any_flares_`t'm .=0
	lab var any_flares_`t'm  "At least one additional flare within `t' months after diagnosis"
	lab def any_flares_`t'm  0 "No" 1 "Yes", modify
	lab val any_flares_`t'm  any_flares_`t'm 
	tab any_flares_`t'm
}

**Treatment of first flare after diagnosis (i.e. received colchicine vs. NSAIDs vs. steroids on same day)
gen first_flare_drug = 0 if first_flare_overall_date!=.

forval i = 1/$max_prescription {
	replace first_flare_drug = 1 if first_flare_overall_date!=. & first_flare_overall_date == colchicine_date_`i' & colchicine_date_`i'!=. 
	replace first_flare_drug = 2 if first_flare_overall_date!=. & first_flare_overall_date == nsaid_date_`i' & nsaid_date_`i'!=. 
	replace first_flare_drug = 3 if first_flare_overall_date!=. & first_flare_overall_date == steroid_date_`i' & steroid_date_`i'!=. 
}

lab var first_flare_drug "Drug used for treatment of gout flare"
lab define first_flare_drug 0 "No drug" 1 "Colchicine" 2 "NSAID" 3 "Corticosteroid", modify
lab val first_flare_drug first_flare_drug

**Whether urate test was performed within 3 months of first non-index flare
generate post_flare_urate = 0 if first_flare_overall_date!=.

***Find max number of urate levels
local max = 0
capture quietly ds urate_date_*
if !_rc & "`r(varlist)'" != "" {
	foreach v of varlist `r(varlist)' {
		if regexm("`v'","^urate_date_([0-9]+)$") {
			local idx = real(regexs(1))
			local max = max(`max', `idx')
		}
	}
}
di "`max'"

***Find matching urate levels within 3 months after baseline urate (if baseline <360)
forval i = 1/`max'	{
	replace post_flare_urate = 1 if first_flare_overall_date!=. & ((urate_date_`i' > first_flare_overall_date) & (urate_date_`i' <= (first_flare_overall_date + 90)) & urate_date_`i'!=.)
}
lab var post_flare_urate "Urate level within three months of first non-index flare"
lab def post_flare_urate 0 "No" 1 "Yes"
lab val post_flare_urate post_flare_urate

*Work out when patients would be classified as "offer" ULT===========

**Multiple or troublesome flares = first flare more than 14 days after diagnosis = first_flare_overall_date
**CKD stages 3 to 5 (glomerular filtration rate [GFR] categories G3 to G5) = first appearance of CKD code or single eGFR <60 = first_ckd_comb_date 
**Diuretic therapy = use of diuretic within 6m before diagnosis or first use after diagnosis = diuretic_bl_date or diuretic_date_1
**Tophi = first appearance of tophaceous gout code = tophi_date
**Chronic gouty arthritis = first appearance of chronic gouty arthritis code = chronic_gout_date

**Binary variable - at any point during study period or before
gen ult_risk_ever = 1 if first_flare_overall_date!=. | first_ckd_comb_date!=. | diuretic_bl_date!=. | diuretic_date_1!=. | tophi_date!=. | chronic_gout_date!=.
recode ult_risk_ever .=0	
lab var ult_risk_ever "Presence of ULT risk factors"
lab def ult_risk_ever 0 "No" 1 "Yes", modify
lab val ult_risk_ever ult_risk_ever

**Date they first became at risk
gen ult_risk_date = min(first_flare_overall_date,first_ckd_comb_date,diuretic_bl_date,diuretic_date_1,tophi_date,chronic_gout_date) 
format %td ult_risk_date
lab var ult_risk_date "Onset of ULT risk factors"
gen ult_risk_date_dx = ult_risk_date
replace ult_risk_date_dx = ${disease}_inc_date if (ult_risk_date < ${disease}_inc_date) //recode risk date as primary diagnosis date if occurred before then
format %td ult_risk_date_dx
lab var ult_risk_date_dx "Onset of ULT risk factors"

**Check whether at least 12 months of follow-up after at risk date
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	gen has_`t'm_fup_risk=1 if (reg_end_date!=. & (reg_end_date >= (ult_risk_date_dx + `days')) & ((ult_risk_date_dx + `days') <= (date("$studyfup_date", "YMD")))) | (reg_end_date==. & ((ult_risk_date_dx + `days') <= (date("$studyfup_date", "YMD"))))
	recode has_`t'm_fup_risk .=0 //includes those who were not at risk
	lab var has_`t'm_fup_risk "At least `t' months of follow-up after becoming at risk for ULT"
	lab def has_`t'm_fup_risk 0 "No" 1 "Yes"
	lab val has_`t'm_fup_risk has_`t'm_fup_risk
	tab has_`t'm_fup_risk
}

**ULT risk factors at baseline
gen ult_risk_bl = 1 if ult_risk_date<=${disease}_inc_date & ult_risk_date!=.
recode ult_risk_bl 0=.
lab var ult_risk_bl "Presence of ULT risk factors at diagnosis"
lab def ult_risk_bl 0 "No" 1 "Yes", modify
lab val ult_risk_bl ult_risk_bl

**Check whether ULT prescribed when became high risk (or at diagnosis if before gout diagnosis) - a prescription 6m before or up to 6/12m afterwards
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days' 
	reshape long ult_date_, i(patient_id) j(ult_order)
	gen ult_risk_p_`t'm = 1 if ult_risk_date_dx!=. & ult_date_!=. & (ult_date_ >= (ult_risk_date_dx - 183)) & (ult_date_ <= (ult_risk_date_dx + `days')) 
	sort patient_id ult_risk_p_`t'm
	by patient_id: replace ult_risk_p_`t'm= ult_risk_p_`t'm[_n-1] if missing(ult_risk_p_`t'm)
	reshape wide ult_date_, i(patient_id) j(ult_order)
	recode ult_risk_p_`t'm .=0 if ult_risk_date_dx!=. 
	lab var ult_risk_p_`t'm "Prescribed ULT within `t' months of onset of ULT risk factors"
	lab def ult_risk_p_`t'm 0 "No" 1 "Yes"
	lab val ult_risk_p_`t'm ult_risk_p_`t'm
	tab ult_risk_p_`t'm, missing
}

save "$projectdir/output/data/cohort_events.dta", replace

*Admissions and referrals (amend as necessary) ================================================================*/

use "$projectdir/output/data/cohort_events.dta", clear

**Any recorded appointment with specialty (passed from YAML; can loop if necessary) in the 12 months before diagnosis (restricted by dataset definition)
codebook ${outpatients}_opa_date //check
gen ${outpatients}_opa_before = 1 if (${outpatients}_opa_date < ${disease}_inc_date) & ${outpatients}_opa_date!=.
lab var ${outpatients}_opa_before "Outpatient appointment in ${outpatients} within 12 months before diagnosis date"
lab define ${outpatients}_opa_before 0 "No" 1 "Yes"
recode ${outpatients}_opa_before .=0
lab val ${outpatients}_opa_before ${outpatients}_opa_before
tab ${outpatients}_opa_before, missing

**Any recorded appointment with specialty on or after diagnosis, if no appointment before diagnosis
gen ${outpatients}_opa_after = 1 if (${outpatients}_opa_date >= ${disease}_inc_date) & ${outpatients}_opa_date!=.
lab var ${outpatients}_opa_after "Outpatient appointment in ${outpatients} after diagnosis date"
lab define ${outpatients}_opa_after 0 "No" 1 "Yes"
recode ${outpatients}_opa_after .=0 // includes those who had first appointment before diagnosis
lab val ${outpatients}_opa_after ${outpatients}_opa_after
tab ${outpatients}_opa_after, missing
gen ${outpatients}_opa_after_date = ${outpatients}_opa_date if (${outpatients}_opa_date >= ${disease}_inc_date) & ${outpatients}_opa_date!=.
lab var ${outpatients}_opa_after_date "Date of first outpatient appointment in ${outpatients} after diagnosis date"
format ${outpatients}_opa_after_date %td

**Any recorded referral to specialty (passed from YAML; can loop if necessary) in the 12 months before diagnosis
codebook ${outpatients}_ref_date //check
gen ${outpatients}_ref_before = 1 if (${outpatients}_ref_date < ${disease}_inc_date) & ${outpatients}_ref_date!=.
lab var ${outpatients}_ref_before "Referral to ${outpatients} within 12 months before diagnosis date"
lab define ${outpatients}_ref_before 0 "No" 1 "Yes"
recode ${outpatients}_ref_before .=0
lab val ${outpatients}_ref_before ${outpatients}_ref_before
tab ${outpatients}_ref_before, missing

**Any recorded referral to specialty on or after diagnosis, if no referral before diagnosis
gen ${outpatients}_ref_after = 1 if (${outpatients}_ref_date >= ${disease}_inc_date) & ${outpatients}_ref_date!=.
lab var ${outpatients}_ref_after "Referral to ${outpatients} after diagnosis date"
lab define ${outpatients}_ref_after 0 "No" 1 "Yes"
recode ${outpatients}_ref_after .=0 // includes those who had first referral before diagnosis
lab val ${outpatients}_ref_after ${outpatients}_ref_after
tab ${outpatients}_ref_after, missing
gen ${outpatients}_ref_after_date = ${outpatients}_ref_date if (${outpatients}_ref_date >= ${disease}_inc_date) & ${outpatients}_ref_date!=.
lab var ${outpatients}_ref_after_date "Date of first referral to ${outpatients} after diagnosis date"
format ${outpatients}_ref_after_date %td

**Any recorded referral to specialty, appointment in specialty, or either/or, within 12 months before or after diagnosis
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	gen ${outpatients}_opa_`t'm = 1 if (${outpatients}_opa_date >= (${disease}_inc_date - `days')) & (${outpatients}_opa_date <= (${disease}_inc_date + `days')) & ${outpatients}_opa_date!=.
	lab var ${outpatients}_opa_`t'm "Outpatient appointment in ${outpatients} within `t' months before or after diagnosis date"
	lab define ${outpatients}_opa_`t'm 0 "No" 1 "Yes"
	recode ${outpatients}_opa_`t'm .=0
	lab val ${outpatients}_opa_`t'm ${outpatients}_opa_`t'm
	tab ${outpatients}_opa_`t'm, missing

	gen ${outpatients}_ref_`t'm = 1 if (${outpatients}_ref_date >= (${disease}_inc_date - `days')) & (${outpatients}_ref_date <= (${disease}_inc_date + `days')) & ${outpatients}_ref_date!=.
	lab var ${outpatients}_ref_`t'm "Referral to ${outpatients} within `t' months before or after diagnosis date"
	lab define ${outpatients}_ref_`t'm 0 "No" 1 "Yes"
	recode ${outpatients}_ref_`t'm .=0
	lab val ${outpatients}_ref_`t'm ${outpatients}_ref_`t'm
	tab ${outpatients}_ref_`t'm, missing
	
	gen ${outpatients}_refopa_`t'm = 1 if ${outpatients}_opa_`t'm==1 | ${outpatients}_ref_`t'm==1
	lab var ${outpatients}_refopa_`t'm "Referral and/or appointment with ${outpatients} within `t' months before or after diagnosis date"
	lab define ${outpatients}_refopa_`t'm 0 "No" 1 "Yes"
	recode ${outpatients}_refopa_`t'm .=0
	lab val ${outpatients}_refopa_`t'm ${outpatients}_refopa_`t'm
	tab ${outpatients}_refopa_`t'm, missing
}

**Date of first referral/outpatient appointment within 12 months of diagnosis
gen ${outpatients}_refopa_date = ${outpatients}_opa_date if ${outpatients}_opa_12m==1 & ${outpatients}_ref_12m==0
replace ${outpatients}_refopa_date = ${outpatients}_ref_date if ${outpatients}_ref_12m==1 & ${outpatients}_opa_12m==0
replace ${outpatients}_refopa_date = min(${outpatients}_ref_date, ${outpatients}_opa_date) if ${outpatients}_ref_12m==1 & ${outpatients}_opa_12m==1
format ${outpatients}_refopa_date %td

**Consider referring a person with gout to a rheumatology service if they have CKD stages 3b to 5 or they have had an organ transplant
gen ckd_transplant_bl = 1 if ckd_comb_bl_3b==1 | transplant_bl==1 //at baseline
recode ckd_transplant_bl .=0
lab var ckd_transplant_bl "Chronic kidney disease (stages 3b-5) or solid organ transplant"
lab def ckd_transplant_bl 0 "No" 1 "Yes"
lab val ckd_transplant_bl ckd_transplant_bl
tab ckd_transplant_bl, missing

**Any recorded referral to specialty, appointment in specialty, or either/or, within 12 months before or after diagnosis if the patient had CKD stages 3b to 5 or they have had an organ transplant at diagnosis
foreach t in 6 12 {
	local days = int((`t'/12)*365.25)
	di `days'
	
	gen ${outpatients}_refopa_`t'm_risk = 1 if ${outpatients}_refopa_`t'm==1 & ckd_transplant_bl==1
	replace ${outpatients}_refopa_`t'm_risk = 0 if ${outpatients}_refopa_`t'm==0 & ckd_transplant_bl==1
	lab var ${outpatients}_refopa_`t'm_risk "Referral and/or appointment with ${outpatients} within `t' months before or after diagnosis date"
	lab def ${outpatients}_refopa_`t'm_risk 0 "No" 1 "Yes", modify
	lab val ${outpatients}_refopa_`t'm_risk ${outpatients}_refopa_`t'm_risk
	tab ${outpatients}_refopa_`t'm_risk, missing
}

save "$projectdir/output/data/cohort_processed_outpatients.dta", replace

*Further cleaning for landmark survival analyes=====================

use "$projectdir/output/data/cohort_processed_outpatients.dta", clear

**CKD status with relation to ULT initiation date
gen ckd_pre_ult = 0
replace ckd_pre_ult = 1 if (first_ckd_comb_date <= ult_first_date) & first_ckd_comb_date !=. & ult_first_date !=.
replace ckd_pre_ult =. if ult_first_date ==.
label def ckd_pre_ult 0 "No" 1 "Yes"
label val ckd_pre_ult ckd_pre_ult
label var ckd_pre_ult "Evidence of CKD before ULT initiation"

gen ckd_free_ult = 0
replace ckd_free_ult = 1 if ((first_ckd_comb_date > ult_first_date) & first_ckd_comb_date !=. & ult_first_date !=.) | first_ckd_comb_date ==.
replace ckd_free_ult =. if ult_first_date ==.
label def ckd_free_ult 0 "No" 1 "Yes"
label val ckd_free_ult ckd_free_ult
label var ckd_free_ult "No evidence of CKD at ULT initiation"

**Landmark date: 12 months after ULT initiation
gen ult_landmark = (ult_first_date + 365) if ult_first_date !=.
format ult_landmark %td
label var ult_landmark "ULT initiation + 12 months"

**CKD status with relation to landmark date
gen ckd_pre_landmark = 0
replace ckd_pre_landmark = 1 if (first_ckd_comb_date <= ult_landmark) & first_ckd_comb_date !=. & ult_landmark !=.
replace ckd_pre_landmark =. if ult_landmark ==.
label def ckd_pre_landmark 0 "No" 1 "Yes"
label val ckd_pre_landmark ckd_pre_landmark
label var ckd_pre_landmark "Evidence of CKD before ULT landmark"

gen ckd_free_landmark = 0
replace ckd_free_landmark = 1 if ((first_ckd_comb_date > ult_landmark) & first_ckd_comb_date !=. & ult_landmark !=.) | first_ckd_comb_date ==.
replace ckd_free_landmark =. if ult_landmark ==.
label def ckd_free_landmark 0 "No" 1 "Yes"
label val ckd_free_landmark ckd_free_landmark
label var ckd_free_landmark "No evidence of CKD at ULT landmark"

**Covariates with relation to landmark date

***Age at landmark
gen time_diag_land = (ult_landmark - ${disease}_inc_date) if ult_landmark!=.
gen time_diag_land_yr = time_diag_land/365.25 if ult_landmark!=.
gen age_land = (age + time_diag_land_yr) if ult_landmark!=.
gen age_land_decile = (age_land/10) if ult_landmark!=.
lab var age_land_decile "Age at landmark, decile"
drop time_diag_land  time_diag_land_yr

***Comorbidities at landmark
foreach comorbidity in $comorbidities {
    local lbl : subinstr local comorbidity "_" " ", all
	local lbl = strproper("`lbl'")
	di "`lbl'"
	
	gen `comorbidity'_land = 0
	replace `comorbidity'_land = 1 if (`comorbidity'_date <= ult_landmark) & `comorbidity'_date!=. & ult_landmark!=.
	replace `comorbidity'_land =. if ult_landmark ==.
	lab define `comorbidity'_land 0 "No" 1 "Yes", modify
	lab var `comorbidity'_land "`lbl' at or before landmark"
	lab val `comorbidity'_land `comorbidity'_land
	tab `comorbidity'_land, missing
}

lab var ckd_land "CKD at or before landmark"
lab var diabetes_land "T2DM at or before landmark"
lab var chd_land "CHD at or before landmark"
lab var cva_land "Stroke/TIA at or before landmark"
lab var liver_disease_land "Chronic liver disease at or before landmark"
lab var transplant_land "Solid organ transplant at or before landmark"
lab var alcohol_land "Excess alcohol at or before landmark"

***Drugs at landmark (prescriptions within 6 months before)
foreach drug in diuretic sglt2 ace_arb {
	
	local Drug = strproper("`drug'") //first letter capitalised for labelling
	di "`drug'"

	gen `drug'_land = 0 if !missing(ult_landmark)

	forval i = 1/$max_prescription {
		replace `drug'_land = 1 if !missing(`drug'_date_`i') & !missing(ult_landmark) & (`drug'_date_`i' >= (ult_landmark - 183)) & (`drug'_date_`i' <= ult_landmark)
	}
	lab var `drug'_land "`Drug' prescription within 6m before landmark"
	lab def `drug'_land 0 "No" 1 "Yes", modify
	lab val `drug'_land `drug'_land
}

save "$projectdir/output/data/cohort_processed_prepractice.dta", replace

*Generate practice-level summary counts=================================

**Import practice-level measures
import delimited "$projectdir/output/measures/measures_practice_$disease.csv", clear
rename numerator practice_${disease}_n
rename denominator practice_list_n
lab var practice_list_n "Practice list size"
rename ratio practice_${disease}_ratio
lab var practice_${disease}_ratio "Disease to list size ratio"
order practice_id, first
keep practice*
save "$projectdir/output/data/measures_practice.dta", replace

**Merge with cleaned dataset
use "$projectdir/output/data/cohort_processed_prepractice.dta", clear
merge m:1 practice_id using "$projectdir/output/data/measures_practice.dta"
drop if _merge==2
drop _merge
tabstat practice_${disease}_n, stat(n mean sd p50 p25 p75)
tabstat practice_list_n, stat(n mean sd p50 p25 p75)
tabstat practice_${disease}_ratio, stat(n mean sd p50 p25 p75)

save "$projectdir/output/data/cohort_processed.dta", replace

log close
