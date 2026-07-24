version 16

/*==============================================================================
DO FILE NAME:			Survival model analyses
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
log using "$logdir/survival_models.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list (passed from yaml)
global arglist disease studyfup_date

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
	global studyfup_date "2026-06-30"
}
di "$disease"

set type double

set scheme plotplainblind

*Define programme to run Cox models and output values of interest =============

capture program drop cox_model

program define cox_model, rclass

	**Model arguments
	args model_terms focal_predictor outcome outlabel model_label
	di as txt "Model terms = `model_terms'"
	
	**Run model
	capture noisily stcox `model_terms', vce(cluster practice_id)
	
	**Skip if estimation failed
	if _rc {
		di as txt "Skipping model (estimation failure): `model_terms'"
		return scalar model_ok = 0
		exit
	}
	
	**Check to ensure model ran ok
	return scalar model_ok = 1
	
	**Store number of patients and practices, person-years of follow-up and n_events (all rounded)and degrees of freedom
	local n_patients = round(e(N), 5)

	local n_practices = .
	capture confirm scalar e(N_clust)
	if !_rc {
		local n_practices = round(e(N_clust), 5)
	}
	if missing(`n_practices') {
		tempvar tag_practice
		egen `tag_practice' = tag(practice_id) if e(sample)
		quietly count if `tag_practice'
		local n_practices = round(r(N), 5)
	}
	
	quietly count if e(sample) & _d == 1 
	local n_events = round(r(N), 5)
	
	quietly summarize _t if e(sample), meanonly 
	local person_years = round(r(sum), 5) //because stset scale is 365.25, this equates to years

	local df = e(df_m)
	
	**Strip factor prefix from focal predictor (if present)
    local focalvar "`focal_predictor'"
    local focalvar = subinstr("`focalvar'", "i.", "", .)
    local focalvar = subinstr("`focalvar'", "c.", "", .)
	
	**Store outputs from model
	matrix B = e(b)
	local cnames : colnames B
	
	**Cycle through column names
	foreach term of local cnames {

    **Skip intercepts (if present)
    if "`term'" == "_cons" continue

    **Store defaults
    local varname "`term'"
    local category "Continuous"
    local levelnum ""
    local omitted = 0
    local base = 0

    **Handle omitted terms
    if regexm("`term'", "^([0-9]+)o\.(.+)$") {
        local levelnum "`=regexs(1)'"
        local varname  "`=regexs(2)'"
        local omitted = 1
    }
    else if regexm("`term'", "^o\.(.+)$") {
        local varname "`=regexs(1)'"
        local category "Omitted"
        local omitted = 1
    }
	
    **Handle base factor terms
    else if regexm("`term'", "^([0-9]+)b\.(.+)$") {
        local levelnum "`=regexs(1)'"
        local varname  "`=regexs(2)'"
        local base = 1
    }

    **Handle regular terms
    else if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
        local levelnum "`=regexs(1)'"
        local varname  "`=regexs(3)'"
    }

    **Store factor level label, if applicable
    if "`levelnum'" != "" {
        local labname : value label `varname'
        if "`labname'" != "" {
            capture local category : label `labname' `levelnum'
            if _rc local category "`levelnum'"
        }
        else {
            local category "`levelnum'"
        }
    }

    **Annotate omitted terms
    if `omitted' == 1 {
        if "`category'" == "Continuous" local category "Omitted"
        else local category "`category' (omitted)"
    }

    **Restrict output to focal predictor variable
    if "`focalvar'" != "" {
        if "`varname'" != "`focalvar'" continue
    }

    **Store variable label
    local varlabel : variable label `varname'
    if "`varlabel'" == "" local varlabel "`varname'"

    **Post omitted terms
    if `omitted' == 1 {
        post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (.) (.) (.) (.)
        continue
    }

    **Post base factor terms
    if `base' == 1 {
        post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (1) (.) (.) (.)
        continue
    }

    **Output estimates for terms of interest
    capture scalar b = _b[`term']
    if _rc continue

    capture scalar se = _se[`term']
    if _rc continue
    if missing(se) continue
    if se == 0 continue

    ***Calculate HR, CI, p-values
    scalar hr = exp(b)
    scalar lo = exp(b - invnormal(0.975)*se)
    scalar hi = exp(b + invnormal(0.975)*se)
    scalar pv = 2*normal(-abs(b/se))

    local hazardratio = round(hr, 0.001)
    local lower95 = round(lo, 0.001)
    local upper95 = round(hi, 0.001)
    local pvalue = round(pv, 0.0001)

    **Post model results
    post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (`hazardratio') (`lower95') (`upper95') (`pvalue')
	}
end

*Load processed cohort ================================
use "$projectdir/output/data/cohort_processed.dta", clear

local n_km_graphs = 0
local n_loglog_graphs = 0

capture erase "$projectdir/output/figures/km_no_outputs.$img"
capture erase "$projectdir/output/figures/loglog_no_outputs.$img"

*Define key variables for landmark survival analysis ===============================================

**Define cohort entry date
local cohort_entry_date ult_first_date //date of first ULT drug

**Define landmark date
local landmark_date ult_landmark //date of first ULT drug + 12 months

**Define censor criteria
gen study_end = date("$studyfup_date", "YMD")
format study_end %td
local study_end_date study_end //end of study follow-up period
local death_date date_of_death //date of death
local dereg_date reg_end_date //end of practice registration
egen censor_date = rowmin(`study_end_date' `death_date' `dereg_date') //first of the above dates
format censor_date %td
label var censor_date "Censoring date"

**Primary exposure variable
local exposure_primary urate_12m_ult //urate checked and target attained vs. not attained within 12 months of ULT initiation (coded as 1/0/missing)

**Sensitivity exposure variables
local exposure_primary_3cat urate_12m_ult_cat //separate category coded if urate not checked (1/0/9)
local exposure_sens_nomiss urate_12m_ult_recode //recoded as not attained if urate not checked (coded as 1/0)

**Define exposure list to loop through
local exposures `exposure_primary' `exposure_primary_3cat' `exposure_sens_nomiss'

**Primary outcome
gen sec_ckd_egfr_land_date = second_egfr_ckd_date if (second_egfr_ckd_date > `landmark_date') & second_egfr_ckd_date !=. & `landmark_date' !=.
format sec_ckd_egfr_land_date %td
label var sec_ckd_egfr_land_date "Incident CKD by two eGFRs <60 after ULT landmark"

**Sensitivity outcomes
gen first_ckd_egfr_land_date = first_egfr_ckd_date if (first_egfr_ckd_date > `landmark_date') & first_egfr_ckd_date !=. & `landmark_date' !=.
format first_ckd_egfr_land_date %td
label var first_ckd_egfr_land_date "Incident CKD by one eGFR <60 after ULT landmark"

gen first_ckd_code_land_date = ckd_date if (ckd_date > `landmark_date') & ckd_date !=. & `landmark_date' !=.
format first_ckd_code_land_date %td
label var first_ckd_code_land_date "Incident CKD by one CKD code after ULT landmark"

**Define outcome list to loop through
local outcomes sec_ckd_egfr_land_date first_ckd_egfr_land_date first_ckd_code_land_date

**Outcome status at baseline/landmark variables
local outcome_free_baseline ckd_free_ult //CKD, defined using single eGFR <60 or CKD code at or before ULT initiation date
local outcome_free_landmark ckd_free_landmark //CKD, defined using single eGFR <60 or CKD code at or before ULT initiation date + 12 months

**Define patient-level predictors
local patient_predictors_core ///
    age_land_decile i.sex i.imd i.ethnicity i.bmicat i.smoke i.diabetes_land i.heart_failure_land i.chd_land i.cva_land i.hypertension_land i.alcohol_land i.diuretic_land i.sglt2_land i.ace_arb_land
    *rheum_appt_n_12m hosp_n_12m creatinine_n_12m
	
local patient_predictors_extra ///
	urate_before_ult_value egfr_before_ult_value

*Run landmark Cox models =======================================================

**Generate temporary file to store outputs
tempname cox_measures
postfile `cox_measures' str80(outcome) str80(outcome_label) str80(exposure) str80(exposure_category) str20(model) double n_patients n_practices n_events person_years df hazardratio lower95 upper95 pvalue ///
    using "$projectdir/output/data/landmark_cox_summary.dta", replace
	
global cox_measures `cox_measures'

capture stset, clear

/******Test criteria: remove********
local cohort_entry_date flare_overall_date_1
local landmark_date flare_overall_date_1
local exposure_primary ckd_comb
local exposure_primary_3cat imd
local exposure_sens_nomiss ethnicity
local exposures `exposure_primary' `exposure_primary_3cat' `exposure_sens_nomiss'
gen test = 1
local outcome_free_baseline test 
local outcome_free_landmark test 
local outcomes nsaid_last_date gout_adm_date_1     
*/

preserve

**Cohort of interest (exposure status set below)
keep if !missing(`landmark_date') & !missing(censor_date) & (censor_date > `landmark_date') //landmark date present and before censor date
keep if `outcome_free_baseline' ==1 //outcome not present before cohort entry
keep if `outcome_free_landmark' ==1 //outcome not present before landmark

**Loop through outcomes
foreach outcome of local outcomes {
	
	di as txt "Outcome = `outcome'"
			
	***Store outcome variable name
	local outlabel : variable label `outcome'
	if "`outlabel'" == "" local outlabel "`outcome'"
	
	***Clear previous values
	capture drop stop_date fail
	capture stset, clear
	
	***Assign fail and stop dates
	gen stop_date = censor_date
	replace stop_date = `outcome' if !missing(`outcome') & (`outcome' <= censor_date)
	format stop_date %td
	gen fail = !missing(`outcome') & (`outcome' <= censor_date)

	***Set survival model
	stset stop_date, origin(time `landmark_date') scale(365.25) failure(fail == 1)
		*id(patient_id)
	
	***Failsafe if no observations
	quietly count if _d == 1
	local nfail = r(N)

	if `nfail' == 0 {
		di as txt "No failures for `outcome'; skipping all analyses."
		continue
	}
	
	***Loop through exposures (primary vs. sensitivity analyses)
	foreach exposure of local exposures {
		
		di as txt "Exposure = `exposure'"
		
		***Failsafe if exposure has no observations
		quietly count if !missing(`exposure')
		if r(N) == 0 {
			di as txt "No non-missing observations for `exposure'; skipping."
			continue
		}

		***Failsafe if exposure-specific sample has no failures
		quietly count if !missing(`exposure') & _d == 1
		if r(N) == 0 {
			di as txt "No failures among patients with non-missing `exposure'; skipping."
			continue
		}

		****Run univariable model
		local model_terms i.`exposure' if !missing(`exposure')
		cox_model `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"Unadjusted"'
		
		****Run age and sex-adjusted model
		local model_terms i.`exposure' age_land_decile i.sex if !missing(`exposure')
		cox_model `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"Age/sex-adjusted"' 

		****Run multivariable model
		local model_terms i.`exposure' `patient_predictors_core' if !missing(`exposure')
		cox_model `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"Multivariable core"' 
		
		****Run multivariable model with baseline urate and eGFR (values closest to before ULT initiation, but within 12m)
		local model_terms i.`exposure' `patient_predictors_core' `patient_predictors_extra' if !missing(`exposure')
		cox_model `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"Multivariable extra"'
		
		****Output KM and loglog plots
		
		*****Store labels for graph
		levelsof `exposure' if !missing(`exposure'), local(levels)

		local colours "emerald orange red blue dkgreen cranberry navy maroon teal sienna purple"
		local legtitle : variable label `exposure'
		if "`legtitle'" == "" local legtitle "`exposure'"

		local i = 1
		local legorder
		local km_plotopts
		local loglog_plotopts

		foreach l of local levels {
			local lab : label (`exposure') `l'
			if "`lab'" == "" local lab "`l'"

			local legorder `legorder' `i' "`lab'"

			local col : word `i' of `colours'
			if "`col'" == "" local col "black"

			local km_plotopts `km_plotopts' plot`i'opts(lcolor(`col') lpattern(solid))
			local loglog_plotopts `loglog_plotopts' plot`i'opts(lcolor(`col') lpattern(solid) msymbol(i))

			local ++i
		}
		
		*****Naming of graphs
		local graphstub = substr("`exposure'_`outcome'", 1, 25)
		local kmname "km_`graphstub'"
		local loglogname "ll_`graphstub'"
		
		*****Survival plot
		sts graph if !missing(`exposure'), by(`exposure') survival `km_plotopts' ytitle("Survival probability", size(medsmall)) ylabel(, nogrid labsize(small)) xtitle("Years from landmark", size(medsmall) margin(medsmall)) xlabel(, nogrid labsize(small)) title("", size(medium) margin(b=2)) legend(order(`legorder') title("`legtitle'", size(small) margin(b=1))) name(`kmname', replace) saving("$projectdir/output/figures/km_`exposure'_`outcome'.gph", replace)
		capture graph export "$projectdir/output/figures/km_`exposure'_`outcome'.$img", replace
		
		if _rc == 0 {
			local ++n_km_graphs
		}
		
		*****Log-log plot
		stphplot if !missing(`exposure'), by(`exposure') `loglog_plotopts' ytitle("log{-log(Survival probability)}", size(medsmall)) ylabel(, nogrid labsize(small)) xtitle("log(Time)", size(medsmall) margin(medsmall)) xlabel(, nogrid labsize(small)) title("", size(medium) margin(b=2)) legend(order(`legorder') title("`legtitle'", size(small) margin(b=1))) name(`loglogname', replace) saving("$projectdir/output/figures/loglog_`exposure'_`outcome'.gph", replace)
		capture graph export "$projectdir/output/figures/loglog_`exposure'_`outcome'.$img", replace
		
		if _rc == 0 {
			local ++n_loglog_graphs
		}
	}
}
	
restore

*Close tempfile
postclose $cox_measures

*Output postfile to csv - with failsafes
capture use "$projectdir/output/data/landmark_cox_summary.dta", clear

if _rc {
	clear
	set obs 0
	gen str1 outcome = ""
}

format hazardratio lower95 upper95 %9.3f
format pvalue %9.4f

export delimited using "$projectdir/output/tables/landmark_cox_summary.csv", replace

*Create dummy KM graph only if no KM graphs were exported
if `n_km_graphs' == 0 {
    preserve
    clear
    set obs 1
    gen x = 1
    gen y = 1

    twoway scatter y x, msymbol(none) xlabel(none) ylabel(none) xtitle("") ytitle("") title("No Kaplan-Meier estimates available") legend(off)

    graph export "$projectdir/output/figures/km_no_outputs.$img", replace
    restore
}

*Create dummy log-log graph only if no log-log graphs were exported
if `n_loglog_graphs' == 0 {
    preserve
    clear
    set obs 1
    gen x = 1
    gen y = 1

    twoway scatter y x, msymbol(none) xlabel(none) ylabel(none) xtitle("") ytitle("") title("No log-log estimates available") legend(off)

    graph export "$projectdir/output/figures/loglog_no_outputs.$img", replace
    restore
}

log close
