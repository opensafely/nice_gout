version 16

/*==============================================================================
DO FILE NAME:			Incidence cleaning
PROJECT:				OpenSAFELY NICE 
DATE: 					29/09/2025
AUTHOR:					M Russell									
DATASETS USED:			Incidence and Measures files
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
capture mkdir "$projectdir/output/tables"
capture mkdir "$projectdir/output/figures"

*Open a log file
global logdir "$projectdir/logs"
cap log close
log using "$logdir/incidence_cleaning.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list and study dates (passed from yaml)
global arglist diseases studystart_date studyend_date
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
	global studystart_date "2016-07-01"
	global studyend_date "2025-06-30"
}

di "$diseases"
di "$studystart_date"
di "$studyend_date"

set type double

set scheme plotplainblind

*Import measures files ================================================================*/

**Derive years from study dates
local start_year = real(substr("$studystart_date", 1, 4))
local end_year = real(substr("$studyend_date", 1, 4)) - 1
local years
forvalues y = `start_year'/`end_year' {
    local years "`years' `y'"
}
local years: list retokenize years
di "`years'"

local first_year: word 1 of `years'

local first_disease: word 1 of $diseases
di "`first_disease'"

**Import first measures file as base dataset
import delimited "$projectdir/output/measures/measures_dataset_`first_disease'_`first_year'.csv", clear
save "$projectdir/output/data/measures_appended.dta", replace

**Loop over diseases and years
foreach disease in $diseases {
	foreach year in `years' {
		if (("`disease'" != "`first_disease'") | ("`year'" != "`first_year'"))  {
		import delimited "$projectdir/output/measures/measures_dataset_`disease'_`year'.csv", clear
		append using "$projectdir/output/data/measures_appended.dta"
		save "$projectdir/output/data/measures_appended.dta", replace 
		}
	}
}

sort measure interval_start sex age
save "$projectdir/output/data/measures_appended.dta", replace 

*Generate rates ===========================================================*/

use "$projectdir/output/data/measures_appended.dta", clear 

**Format dates
rename interval_start interval_start_s
gen interval_start = date(interval_start_s, "YMD") 
format interval_start %td
drop interval_start_s interval_end

**Month/Year of interval
gen year_diag=year(interval_start)
format year_diag %ty
gen month_diag=month(interval_start)
gen mo_year_diagn=ym(year_diag, month_diag)
format mo_year_diagn %tmMon-CCYY
lab var mo_year_diagn "Month/Year of Diagnosis"

**Code incidence and prevalence measures
gen measure_inc = 1 if substr(measure,-10,.) == "_incidence"
recode measure_inc .=0
gen measure_prev = 1 if substr(measure,-11,.) == "_prevalence"
recode measure_prev .=0

**Code IMD and ethnicity measures
gen measure_inc_imd = 1 if substr(measure,-7,.) == "inc_imd"
recode measure_inc_imd .=0
gen measure_inc_ethnicity = 1 if substr(measure,-8,.) == "inc_ethn"
recode measure_inc_ethnicity .=0
gen measure_prev_imd = 1 if substr(measure,-8,.) == "prev_imd"
recode measure_prev_imd .=0
gen measure_prev_ethnicity = 1 if substr(measure,-9,.) == "prev_ethn"
recode measure_prev_ethnicity .=0

**Code region measures
gen measure_inc_region = 1 if substr(measure,-10,.) == "inc_region"
recode measure_inc_region .=0
gen measure_prev_region = 1 if substr(measure,-11,.) == "prev_region"
recode measure_prev_region .=0

**Code any incidence and prevalence measures
gen measure_inc_any = 1 if measure_inc ==1 | measure_inc_imd==1 | measure_inc_ethnicity==1 | measure_inc_region==1
recode measure_inc_any .=0
gen measure_prev_any = 1 if measure_prev ==1 | measure_prev_imd==1 | measure_prev_ethnicity==1 | measure_prev_region==1
recode measure_prev_any .=0

**Label diseases 
gen diseases_ = substr(measure, 1, strlen(measure) - 10) if measure_inc==1
replace diseases_ = substr(measure, 1, strlen(measure) - 11) if measure_prev==1
replace diseases_ = substr(measure, 1, strlen(measure) - 8) if measure_inc_imd==1
replace diseases_ = substr(measure, 1, strlen(measure) - 9) if measure_inc_ethnicity==1
replace diseases_ = substr(measure, 1, strlen(measure) - 11) if measure_inc_region==1
replace diseases_ = substr(measure, 1, strlen(measure) - 9) if measure_prev_imd==1
replace diseases_ = substr(measure, 1, strlen(measure) - 10) if measure_prev_ethnicity==1
replace diseases_ = substr(measure, 1, strlen(measure) - 12) if measure_prev_region==1
rename diseases_ disease
gen disease_full = strproper(subinstr(disease, "_", " ",.))
order disease_full, after(disease)

*Generate incidence and prevalence by months across ages, sexes, IMD and ethnicity ===================================*/

**For overall cohort
sort disease mo_year_diagn measure
bys disease mo_year_diagn measure: egen numerator_all = sum(numerator)
bys disease mo_year_diagn measure: egen denominator_all = sum(denominator)

***Redact and round counts
replace numerator_all =. if numerator_all<=7 | denominator_all<=7
replace denominator_all =. if numerator_all<=7 | numerator_all==. | denominator_all<=7
replace numerator_all = round(numerator_all, 5)
replace denominator_all = round(denominator_all, 5)

***Generate incidence rates per 100,000 population
gen rate_all = (numerator_all/denominator_all) if (numerator_all!=. & denominator_all!=.)
replace rate_all =. if (numerator_all==. | denominator_all==.)
gen rate_all_100000 = rate_all*100000

**For males
bys disease mo_year_diagn measure: egen numerator_male = sum(numerator) if sex=="male"
bys disease mo_year_diagn measure: egen denominator_male = sum(denominator) if sex=="male"

***Redact and round counts
replace numerator_male =. if numerator_male<=7 | denominator_male<=7
replace denominator_male =. if numerator_male<=7 | numerator_male==. | denominator_male<=7
replace numerator_male = round(numerator_male, 5)
replace denominator_male = round(denominator_male, 5)

***Generate incidence rates per 100,000 population
gen rate_male = (numerator_male/denominator_male) if (numerator_male!=. & denominator_male!=.)
replace rate_male =. if (numerator_male==. | denominator_male==.)
gen rate_male_100000 = rate_male*100000

sort disease mo_year_diagn measure_prev_any measure_inc_any rate_male_100000 
by disease mo_year_diagn measure_prev_any measure_inc_any (rate_male_100000): replace rate_male_100000 = rate_male_100000[_n-1] if missing(rate_male_100000)
sort disease mo_year_diagn measure_prev_any measure_inc_any numerator_male 
by disease mo_year_diagn measure_prev_any measure_inc_any (numerator_male): replace numerator_male = numerator_male[_n-1] if missing(numerator_male)
sort disease mo_year_diagn measure_prev_any measure_inc_any denominator_male 
by disease mo_year_diagn measure_prev_any measure_inc_any (denominator_male): replace denominator_male = denominator_male[_n-1] if missing(denominator_male)

**For females
bys disease mo_year_diagn measure: egen numerator_female = sum(numerator) if sex=="female"
bys disease mo_year_diagn measure: egen denominator_female = sum(denominator) if sex=="female"

***Redact and round counts
replace numerator_female =. if numerator_female<=7 | denominator_female<=7
replace denominator_female =. if numerator_female<=7 | numerator_female==. | denominator_female<=7
replace numerator_female = round(numerator_female, 5)
replace denominator_female = round(denominator_female, 5)

***Generate incidence rates per 100,000 population
gen rate_female = (numerator_female/denominator_female) if (numerator_female!=. & denominator_female!=.)
replace rate_female =. if (numerator_female==. | denominator_female==.)
gen rate_female_100000 = rate_female*100000

sort disease mo_year_diagn measure_prev_any measure_inc_any rate_female_100000 
by disease mo_year_diagn measure_prev_any measure_inc_any (rate_female_100000): replace rate_female_100000 = rate_female_100000[_n-1] if missing(rate_female_100000)
sort disease mo_year_diagn measure_prev_any measure_inc_any numerator_female 
by disease mo_year_diagn measure_prev_any measure_inc_any (numerator_female): replace numerator_female = numerator_female[_n-1] if missing(numerator_female)
sort disease mo_year_diagn measure_prev_any measure_inc_any denominator_female 
by disease mo_year_diagn measure_prev_any measure_inc_any (denominator_female): replace denominator_female = denominator_female[_n-1] if missing(denominator_female)

**For age groups (18+); if including 0-18, need to amend below and dataset definitions
replace age = "80" if age == "age_greater_equal_80"

foreach var in 18_29 30_39 40_49 50_59 60_69 70_79 80 {
	replace age = "`var'" if age=="age_`var'"
	bys disease mo_year_diagn measure: egen numerator_`var' = sum(numerator) if age=="`var'"
	bys disease mo_year_diagn measure: egen denominator_`var' = sum(denominator) if age=="`var'"

	***Redact and round counts
	replace numerator_`var' =. if numerator_`var'<=7 | denominator_`var'<=7
	replace denominator_`var' =. if numerator_`var'<=7 | numerator_`var'==. | denominator_`var'<=7
	replace numerator_`var' = round(numerator_`var', 5)
	replace denominator_`var' = round(denominator_`var', 5)

	***Generate incidence rates per 100,000 population
	gen rate_`var' = (numerator_`var'/denominator_`var') if (numerator_`var'!=. & denominator_`var'!=.)
	replace rate_`var' =. if (numerator_`var'==. | denominator_`var'==.)
	gen rate_`var'_100000 = rate_`var'*100000

	sort disease mo_year_diagn measure_prev_any measure_inc_any rate_`var'_100000 
	by disease mo_year_diagn measure_prev_any measure_inc_any (rate_`var'_100000): replace rate_`var'_100000 = rate_`var'_100000[_n-1] if missing(rate_`var'_100000)
	sort disease mo_year_diagn measure_prev_any measure_inc_any numerator_`var'
	by disease mo_year_diagn measure_prev_any measure_inc_any (numerator_`var'): replace numerator_`var' = numerator_`var'[_n-1] if missing(numerator_`var')
	sort disease mo_year_diagn measure_prev_any measure_inc_any denominator_`var' 
	by disease mo_year_diagn measure_prev_any measure_inc_any (denominator_`var'): replace denominator_`var' = denominator_`var'[_n-1] if missing(denominator_`var')
}

**For ethnicity
bys disease mo_year_diagn measure: egen numerator_white = sum(numerator) if ethnicity=="White"
bys disease mo_year_diagn measure: egen denominator_white = sum(denominator) if ethnicity=="White"

bys disease mo_year_diagn measure: egen numerator_mixed = sum(numerator) if ethnicity=="Mixed"
bys disease mo_year_diagn measure: egen denominator_mixed = sum(denominator) if ethnicity=="Mixed"

bys disease mo_year_diagn measure: egen numerator_black = sum(numerator) if ethnicity=="Black or Black British"
bys disease mo_year_diagn measure: egen denominator_black = sum(denominator) if ethnicity=="Black or Black British"

bys disease mo_year_diagn measure: egen numerator_asian = sum(numerator) if ethnicity=="Asian or Asian British"
bys disease mo_year_diagn measure: egen denominator_asian = sum(denominator) if ethnicity=="Asian or Asian British"

bys disease mo_year_diagn measure: egen numerator_other = sum(numerator) if ethnicity=="Chinese or Other Ethnic Groups"
bys disease mo_year_diagn measure: egen denominator_other = sum(denominator) if ethnicity=="Chinese or Other Ethnic Groups"

bys disease mo_year_diagn measure: egen numerator_ethunk = sum(numerator) if ethnicity=="Unknown"
bys disease mo_year_diagn measure: egen denominator_ethunk = sum(denominator) if ethnicity=="Unknown"

***Redact and round counts
foreach var in white mixed black asian other ethunk {
	replace numerator_`var' =. if numerator_`var'<=7 | denominator_`var'<=7
	replace denominator_`var' =. if numerator_`var'<=7 | numerator_`var'==. | denominator_`var'<=7
	replace numerator_`var' = round(numerator_`var', 5)
	replace denominator_`var' = round(denominator_`var', 5)
	
	***Generate incidence rates per 100,000 population
	gen rate_`var' = (numerator_`var'/denominator_`var') if (numerator_`var'!=. & denominator_`var'!=.)
	replace rate_`var' =. if (numerator_`var'==. | denominator_`var'==.)
	gen rate_`var'_100000 = rate_`var'*100000

	sort disease mo_year_diagn measure_prev_any measure_inc_any rate_`var'_100000 
	by disease mo_year_diagn measure_prev_any measure_inc_any (rate_`var'_100000): replace rate_`var'_100000 = rate_`var'_100000[_n-1] if missing(rate_`var'_100000)
	sort disease mo_year_diagn measure_prev_any measure_inc_any numerator_`var'
	by disease mo_year_diagn measure_prev_any measure_inc_any (numerator_`var'): replace numerator_`var' = numerator_`var'[_n-1] if missing(numerator_`var')
	sort disease mo_year_diagn measure_prev_any measure_inc_any denominator_`var' 
	by disease mo_year_diagn measure_prev_any measure_inc_any (denominator_`var'): replace denominator_`var' = denominator_`var'[_n-1] if missing(denominator_`var')
}

*For IMD
bys disease mo_year_diagn measure: egen numerator_imd1 = sum(numerator) if imd=="1 (most deprived)"
bys disease mo_year_diagn measure: egen denominator_imd1 = sum(denominator) if imd=="1 (most deprived)"

bys disease mo_year_diagn measure: egen numerator_imd2 = sum(numerator) if imd=="2"
bys disease mo_year_diagn measure: egen denominator_imd2 = sum(denominator) if imd=="2"

bys disease mo_year_diagn measure: egen numerator_imd3 = sum(numerator) if imd=="3"
bys disease mo_year_diagn measure: egen denominator_imd3 = sum(denominator) if imd=="3"

bys disease mo_year_diagn measure: egen numerator_imd4 = sum(numerator) if imd=="4"
bys disease mo_year_diagn measure: egen denominator_imd4 = sum(denominator) if imd=="4"

bys disease mo_year_diagn measure: egen numerator_imd5 = sum(numerator) if imd=="5 (least deprived)"
bys disease mo_year_diagn measure: egen denominator_imd5 = sum(denominator) if imd=="5 (least deprived)"

bys disease mo_year_diagn measure: egen numerator_imdunk = sum(numerator) if imd=="Unknown"
bys disease mo_year_diagn measure: egen denominator_imdunk = sum(denominator) if imd=="Unknown"

***Redact and round counts
foreach var in imd1 imd2 imd3 imd4 imd5 imdunk {
	replace numerator_`var' =. if numerator_`var'<=7 | denominator_`var'<=7
	replace denominator_`var' =. if numerator_`var'<=7 | numerator_`var'==. | denominator_`var'<=7
	replace numerator_`var' = round(numerator_`var', 5)
	replace denominator_`var' = round(denominator_`var', 5)
	
	***Generate incidence rates per 100,000 population
	gen rate_`var' = (numerator_`var'/denominator_`var') if (numerator_`var'!=. & denominator_`var'!=.)
	replace rate_`var' =. if (numerator_`var'==. | denominator_`var'==.)
	gen rate_`var'_100000 = rate_`var'*100000

	sort disease mo_year_diagn measure_prev_any measure_inc_any rate_`var'_100000 
	by disease mo_year_diagn measure_prev_any measure_inc_any (rate_`var'_100000): replace rate_`var'_100000 = rate_`var'_100000[_n-1] if missing(rate_`var'_100000)
	sort disease mo_year_diagn measure_prev_any measure_inc_any numerator_`var'
	by disease mo_year_diagn measure_prev_any measure_inc_any (numerator_`var'): replace numerator_`var' = numerator_`var'[_n-1] if missing(numerator_`var')
	sort disease mo_year_diagn measure_prev_any measure_inc_any denominator_`var' 
	by disease mo_year_diagn measure_prev_any measure_inc_any (denominator_`var'): replace denominator_`var' = denominator_`var'[_n-1] if missing(denominator_`var')
}

*For Region
bys disease mo_year_diagn measure: egen numerator_east = sum(numerator) if region=="East"
bys disease mo_year_diagn measure: egen denominator_east = sum(denominator) if region=="East"

bys disease mo_year_diagn measure: egen numerator_eastmid = sum(numerator) if region=="East Midlands"
bys disease mo_year_diagn measure: egen denominator_eastmid = sum(denominator) if region=="East Midlands"

bys disease mo_year_diagn measure: egen numerator_london = sum(numerator) if region=="London"
bys disease mo_year_diagn measure: egen denominator_london = sum(denominator) if region=="London"

bys disease mo_year_diagn measure: egen numerator_northeast = sum(numerator) if region=="North East"
bys disease mo_year_diagn measure: egen denominator_northeast = sum(denominator) if region=="North East"

bys disease mo_year_diagn measure: egen numerator_northwest = sum(numerator) if region=="North West"
bys disease mo_year_diagn measure: egen denominator_northwest = sum(denominator) if region=="North West"

bys disease mo_year_diagn measure: egen numerator_southeast = sum(numerator) if region=="South East"
bys disease mo_year_diagn measure: egen denominator_southeast = sum(denominator) if region=="South East"

bys disease mo_year_diagn measure: egen numerator_southwest = sum(numerator) if region=="South West"
bys disease mo_year_diagn measure: egen denominator_southwest = sum(denominator) if region=="South West"

bys disease mo_year_diagn measure: egen numerator_westmid = sum(numerator) if region=="West Midlands"
bys disease mo_year_diagn measure: egen denominator_westmid = sum(denominator) if region=="West Midlands"

bys disease mo_year_diagn measure: egen numerator_yorkshire = sum(numerator) if region=="Yorkshire and The Humber"
bys disease mo_year_diagn measure: egen denominator_yorkshire = sum(denominator) if region=="Yorkshire and The Humber"

bys disease mo_year_diagn measure: egen numerator_regunk = sum(numerator) if region=="Unknown"
bys disease mo_year_diagn measure: egen denominator_regunk = sum(denominator) if region=="Unknown"

***Redact and round counts
foreach var in east eastmid london northeast northwest southeast southwest westmid yorkshire regunk {
	replace numerator_`var' =. if numerator_`var'<=7 | denominator_`var'<=7
	replace denominator_`var' =. if numerator_`var'<=7 | numerator_`var'==. | denominator_`var'<=7
	replace numerator_`var' = round(numerator_`var', 5)
	replace denominator_`var' = round(denominator_`var', 5)
	
	***Generate incidence rates per 100,000 population
	gen rate_`var' = (numerator_`var'/denominator_`var') if (numerator_`var'!=. & denominator_`var'!=.)
	replace rate_`var' =. if (numerator_`var'==. | denominator_`var'==.)
	gen rate_`var'_100000 = rate_`var'*100000

	sort disease mo_year_diagn measure_prev_any measure_inc_any rate_`var'_100000 
	by disease mo_year_diagn measure_prev_any measure_inc_any (rate_`var'_100000): replace rate_`var'_100000 = rate_`var'_100000[_n-1] if missing(rate_`var'_100000)
	sort disease mo_year_diagn measure_prev_any measure_inc_any numerator_`var'
	by disease mo_year_diagn measure_prev_any measure_inc_any (numerator_`var'): replace numerator_`var' = numerator_`var'[_n-1] if missing(numerator_`var')
	sort disease mo_year_diagn measure_prev_any measure_inc_any denominator_`var' 
	by disease mo_year_diagn measure_prev_any measure_inc_any (denominator_`var'): replace denominator_`var' = denominator_`var'[_n-1] if missing(denominator_`var')
}

save "$projectdir/output/data/processed_nonstandardised.dta", replace

*Calculate the age-standardized incidence rate using age-specific incidence data (European Standard Population 2013); amend if including 0-18 age group (total weight 80,700 currently)
use "$projectdir/output/data/processed_nonstandardised.dta", clear

*Append European Standard Population 2013
gen prop=14200 if age=="18_29"
replace prop=13500 if age=="30_39"
replace prop=14000 if age=="40_49"
replace prop=13500 if age=="50_59"
replace prop=11500 if age=="60_69"
replace prop=9000 if age=="70_79"
replace prop=5000 if age=="80"

*Apply standard population weights and generate standardised incidence and prevalence, overall and by sex
gen rate_100000 = ratio*100000
gen new_value = prop*rate_100000
drop rate_100000

bys disease mo_year_diagn measure: egen sum_new_value_male=sum(new_value) if sex=="male"
gen asr_male = sum_new_value_male/80700
replace asr_male =. if rate_male_100000 ==.
sort disease mo_year_diagn measure asr_male 
by disease mo_year_diagn measure (asr_male): replace asr_male = asr_male[_n-1] if missing(asr_male)

bys disease mo_year_diagn measure: egen sum_new_value_female=sum(new_value) if sex=="female" 
gen asr_female = sum_new_value_female/80700
replace asr_female =. if rate_female_100000 ==. 
sort disease mo_year_diagn measure asr_female 
by disease mo_year_diagn measure (asr_female): replace asr_female = asr_female[_n-1] if missing(asr_female)

bys disease mo_year_diagn measure: egen sum_new_value_all=sum(new_value)
gen asr_all = sum_new_value_all/161400
replace asr_all =. if rate_all_100000 ==. 

*Generate standardised incidence and prevalence, by age group
bys disease mo_year_diagn measure: egen sum_new_value_18_29=sum(new_value) if age=="18_29"
gen asr_18_29 = sum_new_value_18_29/28400
replace asr_18_29 =. if rate_18_29_100000 ==.
sort disease mo_year_diagn measure asr_18_29 
by disease mo_year_diagn measure (asr_18_29): replace asr_18_29 = asr_18_29[_n-1] if missing(asr_18_29)

bys disease mo_year_diagn measure: egen sum_new_value_30_39=sum(new_value) if age=="30_39"
gen asr_30_39 = sum_new_value_30_39/27000
replace asr_30_39 =. if rate_30_39_100000 ==.
sort disease mo_year_diagn measure asr_30_39 
by disease mo_year_diagn measure (asr_30_39): replace asr_30_39 = asr_30_39[_n-1] if missing(asr_30_39)

bys disease mo_year_diagn measure: egen sum_new_value_40_49=sum(new_value) if age=="40_49"
gen asr_40_49 = sum_new_value_40_49/28000
replace asr_40_49 =. if rate_40_49_100000 ==.
sort disease mo_year_diagn measure asr_40_49 
by disease mo_year_diagn measure (asr_40_49): replace asr_40_49 = asr_40_49[_n-1] if missing(asr_40_49)

bys disease mo_year_diagn measure: egen sum_new_value_50_59=sum(new_value) if age=="50_59"
gen asr_50_59 = sum_new_value_50_59/27000
replace asr_50_59 =. if rate_50_59_100000 ==.
sort disease mo_year_diagn measure asr_50_59 
by disease mo_year_diagn measure (asr_50_59): replace asr_50_59 = asr_50_59[_n-1] if missing(asr_50_59)

bys disease mo_year_diagn measure: egen sum_new_value_60_69=sum(new_value) if age=="60_69"
gen asr_60_69 = sum_new_value_60_69/23000
replace asr_60_69 =. if rate_60_69_100000 ==.
sort disease mo_year_diagn measure asr_60_69 
by disease mo_year_diagn measure (asr_60_69): replace asr_60_69 = asr_60_69[_n-1] if missing(asr_60_69)

bys disease mo_year_diagn measure: egen sum_new_value_70_79=sum(new_value) if age=="70_79"
gen asr_70_79 = sum_new_value_70_79/18000
replace asr_70_79 =. if rate_70_79_100000 ==.
sort disease mo_year_diagn measure asr_70_79 
by disease mo_year_diagn measure (asr_70_79): replace asr_70_79 = asr_70_79[_n-1] if missing(asr_70_79)

bys disease mo_year_diagn measure: egen sum_new_value_80=sum(new_value) if age=="80"
gen asr_80 = sum_new_value_80/10000
replace asr_80 =. if rate_80_100000 ==.
sort disease mo_year_diagn measure asr_80 
by disease mo_year_diagn measure (asr_80): replace asr_80 = asr_80[_n-1] if missing(asr_80)

sort disease mo_year_diagn measure age sex
bys measure interval_start: gen n=_n
keep if n==1
drop n

save "$projectdir/output/data/processed_standardised.dta", replace

*Output string version of incidence and prevalence (to stop conversion in excel for big numbers)
use "$projectdir/output/data/processed_standardised.dta", clear

keep if measure_inc==1 | measure_prev==1

foreach var in all male female {
	drop rate_`var'
	rename rate_`var'_100000 rate_`var' //unadjusted incidence rate per 100,000
	rename asr_`var' s_rate_`var' //age and sex-standardised incidence rate per 100,000
	order s_rate_`var', after(rate_`var')
	format s_rate_`var' %14.4f
	format rate_`var' %14.4f
	format numerator_`var' %14.0f
	format denominator_`var' %14.0f
}

foreach var in 18_29 30_39 40_49 50_59 60_69 70_79 80 white mixed black asian other ethunk imd1 imd2 imd3 imd4 imd5 imdunk east eastmid london northeast northwest southeast southwest westmid yorkshire regunk {
	drop rate_`var'
	rename rate_`var'_100000 rate_`var'
	format rate_`var' %14.4f
	order rate_`var', after(denominator_`var')
	format numerator_`var' %14.0f
	format denominator_`var' %14.0f
}

keep disease disease_full measure mo_year_diagn numerator_* denominator_* rate_* s_rate_*
replace measure = "Incidence" if substr(measure,-9,.) == "incidence"
replace measure = "Prevalence" if substr(measure,-10,.) == "prevalence"

foreach dis in $diseases {
	preserve
	keep if disease == "`dis'"
	save "$projectdir/output/tables/redacted_counts_`dis'.dta", replace
	export delimited using "$projectdir/output/tables/redacted_counts_`dis'.csv", datafmt replace
	restore
}

log close	
