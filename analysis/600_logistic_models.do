version 16

/*==============================================================================
DO FILE NAME:			Logistic regression analyses
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
log using "$logdir/logistic_models.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

*Set disease list (passed from yaml)
global arglist disease
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
}
di "$disease"

*Define age variable
global age_var age_decile

set type double

set scheme plotplainblind

*Define programme to run multi-level logistic models and output values of interest===================
capture program drop melogit_model

program define melogit_model, rclass

	**Model arguments
	args outcome model_terms focal_predictor model_label outlabel inclusion measures do_icc

	**Run model
	capture noisily melogit `outcome' `model_terms' || practice_id:, or iterate(100)
	
	**Skip if estimation failed
	if _rc {
		di as txt "Skipping model (estimation failure): `model_terms'"
		return scalar model_ok = 0
		exit
	}

	**Skip if model doesn't converge
	capture confirm scalar e(converged)
	if !_rc {
		if e(converged) == 0 {
			di as txt "Skipping model (no convergence): `model_terms'"
			return scalar model_ok = 0
			exit
		}
	}
	
	**Check to ensure model ran ok
	return scalar model_ok = 1
	
	**Store model-level N and degrees of freedom
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
	
	local df = e(df_m)
	
	**Strip factor prefix from focal predictor
    local focalvar "`focal_predictor'"
    local focalvar = subinstr("`focalvar'", "i.", "", .)
    local focalvar = subinstr("`focalvar'", "c.", "", .)
		
	**Store column names for terms
	local cnames : colnames e(b)
	
	**Cycle through terms
	foreach term of local cnames {
		
		***Skip intercept and random-effect variance parameters
		if "`term'" == "_cons" continue
		if strpos("`term'", "var(") continue

		***Local defaults
		local varname "`term'"
		local category "Continuous"
		local levelnum ""
		local omitted = 0
		local base = 0

		***Handle omitted terms
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

		***Handle base factor terms
		else if regexm("`term'", "^([0-9]+)b\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(2)'"
			local base = 1
		}

		***Handle regular terms
		else if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(3)'"
		}

		***Store factor level label, if applicable
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

		***Annotate omitted terms
		if `omitted' == 1 {
			if "`category'" == "Continuous" local category "Omitted"
			else local category "`category' (omitted)"
		}

		***Skip adjustment terms in age/sex-adjusted models unless they are the focal predictor
		if "`model_label'" == "Age/sex-adjusted" {
			if "`focal_predictor'" != "$age_var" & "`term'" == "$age_var" continue
			if "`focal_predictor'" != "i.sex" & "`varname'" == "sex" continue
		}
		
		***For calendar-time sensitivity models, output only focal predictor
		if "`model_label'" == "Calendar-time adjusted" {
			if "`varname'" != "`focalvar'" continue
		}

		***Store variable label
		local varlabel : variable label `varname'
		if "`varlabel'" == "" local varlabel "`varname'"

		***Post omitted terms
		if `omitted' == 1 {
			post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`n_patients') (`n_practices') (`df') (.) (.) (.) (.)
			continue
		}

		***Post base factor levels
		if `base' == 1 {
			post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`n_patients') (`n_practices') (`df') (1) (.) (.) (.)
			continue
		}

		***Output estimates for terms of interest
		capture scalar b = _b[`outcome':`term']
		if _rc continue

		capture scalar se = _se[`outcome':`term']
		if _rc continue
		if missing(se) continue
		if se == 0 continue
		
		***Calculate OR, CI, p-values
		scalar or = exp(b)
		scalar lo = exp(b - invnormal(0.975)*se)
		scalar hi = exp(b + invnormal(0.975)*se)
		scalar pv = 2*normal(-abs(b/se))

		local oddsratio = round(or, 0.001)
		local lower95 = round(lo, 0.001)
		local upper95 = round(hi, 0.001)
		local pvalue = round(pv, 0.0001)

		***Post model results
		post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`n_patients') (`n_practices') (`df') (`oddsratio') (`lower95') (`upper95') (`pvalue')
	}
	
	***Output ICC (Proportion of the total variation in the outcome attributable to differences between practices)
    
	***Passed argument to run ICC
	if "`do_icc'" == "1" {
	
		****Estimate ICC
		capture noisily estat icc

		****Check estimation/model hasn't failed
		if !_rc {
			capture confirm scalar r(icc2)
			if !_rc {
				local icc = round(r(icc2),0.001)
				
				capture matrix CI = r(ci2)
				if !_rc {
					local icc_lo = round(el(CI,1,1), 0.001)
					local icc_hi = round(el(CI,1,2), 0.001)
				}
				else {
					local icc_lo = .
					local icc_hi = .
				}
				
				*****Post ICC results
				post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("Intra-class correlation") ("Practice-level ICC") (`n_patients') (`n_practices') (`df') (`icc') (`icc_lo') (`icc_hi') (.)
			}
		}
	}
end

*Define programme to run standard logistic models and output values of interest ===================

capture program drop logistic_model

program define logistic_model

	args outcome model_terms focal_predictor model_label outlabel inclusion measures

	**Run logistic model
	capture noisily logistic `outcome' `model_terms', vce(cluster practice_id)
	
	**Skip if estimation failed
	if _rc {
		di as txt "Skipping model (estimation failure): `model_terms'"
		exit
	}
	
	**Store model-level N and degrees of freedom
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
	
	local df = e(df_m)
	
	/*	
	**Store outputs from model (different method)
	matrix T = r(table)
	local cnames : colnames T
	*/
	
	**Strip factor prefix from focal predictor (if present)
    local focalvar "`focal_predictor'"
    local focalvar = subinstr("`focalvar'", "i.", "", .)
    local focalvar = subinstr("`focalvar'", "c.", "", .)
	
	**Store outputs from model
	matrix B = e(b)
	local cnames : colnames B

	**Cycle through column names
	foreach term of local cnames {
		
		***Skip intercepts (if present)
		if "`term'" == "_cons" continue
		
		***Store defaults
		local varname "`term'"
		local category "Continuous"
		local levelnum ""
		local omitted = 0
		local base = 0
		
		***Handle omitted terms
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
		
		***Handle base factor terms
		else if regexm("`term'", "^([0-9]+)b\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(2)'"
			local base = 1
		}
		
		***Handle regular terms
		else if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname  "`=regexs(3)'"
		}

		***Store factor level label, if applicable
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

		***Annotate omitted terms
		if `omitted' == 1 {
			if "`category'" == "Continuous" local category "Omitted"
			else local category "`category' (omitted)"
		}
		
		***Select which variables to output, depending on model
		if "`model_label'" == "Multivariable" {
		}
		else {
			if "`varname'" != "`focalvar'" continue

			****Skip adjustment terms in age/sex-adjusted models unless they are the focal predictor
			if "`model_label'" == "Age/sex-adjusted" {
				if "`focal_predictor'" != "$age_var" & "`term'" == "$age_var" continue
				if "`focal_predictor'" != "i.sex" & "`varname'" == "sex" continue
			}
		}
		
		***For calendar-time sensitivity models, output only focal predictor
		if "`model_label'" == "Calendar-time adjusted" {
			if "`varname'" != "`focalvar'" continue
		}
		
		***Store variable label
		local varlabel : variable label `varname'
		if "`varlabel'" == "" local varlabel "`varname'"

		***Post omitted terms
		if `omitted' == 1 {
			post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`n_patients') (`n_practices') (`df') (.) (.) (.) (.)
			continue
		}

		***Post base factor levels
		if `base' == 1 {
			post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`n_patients') (`n_practices') (`df') (1) (.) (.) (.)
			continue
		}
		
		/*
		***Keep outputs of relevance (different method)
		local j = colnumb(T, "`term'")
		if missing(`j') continue

		local oddsratio = round(T[1,`j'], 0.001)
		local lower95 = round(T[5,`j'], 0.001)
		local upper95 = round(T[6,`j'], 0.001)
		local pvalue = round(T[4,`j'], 0.0001)
		*/
		
		***Output estimates for terms of interest
		capture scalar b = _b[`term']
		if _rc continue

		capture scalar se = _se[`term']
		if _rc continue
		if missing(se) continue
		if se == 0 continue

		***Calculate OR, CI, p-values
		scalar or = exp(b)
		scalar lo = exp(b - invnormal(0.975)*se)
		scalar hi = exp(b + invnormal(0.975)*se)
		scalar pv = 2*normal(-abs(b/se))

		local oddsratio = round(or, 0.001)
		local lower95 = round(lo, 0.001)
		local upper95 = round(hi, 0.001)
		local pvalue = round(pv, 0.0001)

		***Post model results
		post `measures' ("`inclusion'") ("`outcome'") ("`outlabel'") ("`model_label'") ("`varlabel'") ("`category'") (`n_patients') (`n_practices') (`df') (`oddsratio') (`lower95') (`upper95') (`pvalue')
	}
end

*Define common parameters for analyses ===============================================

**Inclusion criteria to loop through
local inclusions has_12m_fup has_12m_fup_ult

**Core patient-level predictors
local patient_predictors_base ///
    $age_var i.sex i.imd i.ethnicity i.bmicat i.smoke i.ckd_comb_bl i.chd_bl i.diabetes_bl i.heart_failure_bl i.cva_bl i.hypertension_bl i.alcohol_bl i.diuretic_bl i.sglt2_bl
	*i.ckd_comb_bl i.chd_bl //this is a test
	*consider baseline_urate too
		
**Core practice-level predictors
local practice_predictors practice_list_n practice_${disease}_ratio

*Run multi-level logistic models ================

**Generate temporary file to store outputs
tempname melogit_measures
postfile `melogit_measures' str80(inclusion) str80(outcome) str80(outcome_label) str80(model) str80(variable) str80(category) double n_patients n_practices df oddsratio lower95 upper95 pvalue ///
    using "$projectdir/output/data/melogit_summary.dta", replace

**Loop through inclusion criteria
foreach inclusion of local inclusions {

	***Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	***Keep those who meet inclusion criteria
	di "Inclusion criterion: `inclusion'"
    keep if `inclusion' == 1

    ***Define outcomes and additional predictors (specific to inclusion criteria); also need to think about interactions between i.post_nice and variables - Nb. programmes about won't output factors properly
	if "`inclusion'" == "has_12m_fup" {
		local outcomes ult_12m
		*local outcomes ckd_comb_bl //this is a test
		local post_nice_var i.post_nice_diag
		local patient_predictors_full `patient_predictors_base' `post_nice_var'
		local time_adjust diagnosis_year

    }
    else if "`inclusion'" == "has_12m_fup_ult" {
		local outcomes urate_12m_ult
        *local outcomes chd_bl //this is a test
		local post_nice_var i.post_nice_ult
		local patient_predictors_full `patient_predictors_base' `post_nice_var'
		local time_adjust ult_year
    }
    else {
        di as error "No outcomes defined for `inclusion'"
        continue
    }
		
	***Loop through outcomes
	foreach outcome of local outcomes {
		
		di "Outcome: `outcome'"
		
		****Store outcome variable name
		local outlabel : variable label `outcome'
		if "`outlabel'" == "" local outlabel "`outcome'"
				
		****Estimate practice-level variation/ICC without predictors ============
		melogit_model `outcome' `""' `""' `"ICC only"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 1

		****Univariable models with patient-level and practice-level predictors ==============
		
		*****Store predictors
		local predictors `patient_predictors_full' `practice_predictors'
		
		*****Loop through predictors
		foreach predictor of local predictors {
			melogit_model `outcome' `"`predictor'"' `"`predictor'"' `"Univariable"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 0
		}
		
		****Age and sex-adjusted models (Nb. don't usually need to present age/sex-adjusted practice-level variables) ============
		
		*****Store predictors
		local predictors `patient_predictors_full' 
		
		*****Loop through predictors
		foreach predictor of local predictors {
			
			local model_terms `predictor' age_decile i.sex
			if "`predictor'" == "$age_var" local model_terms $age_var i.sex
			if "`predictor'" == "i.sex" local model_terms i.sex $age_var
			
			melogit_model `outcome' `"`model_terms'"' `"`predictor'"' `"Age/sex-adjusted"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 0
		}
		
		****Multivariable model with patient-level predictors only ===========		
		local predictors `patient_predictors_full' 
		
		melogit_model `outcome' `"`predictors'"' `""' `"Multivariable patient"' `"`outlabel'"' `"`inclusion'"'  `melogit_measures' 1
		 
		****Multivariable model with patient-level predictors and practice-level predictors ============
		
		*****Store predictors
		local predictors `patient_predictors_full' `practice_predictors'
		
		melogit_model `outcome' `"`predictors'"' `""' `"Multivariable pt/practice"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 1
		 
		****Sensitivity multivariable model with patient-level predictors and practice-level predictors, as well as adjusting for calendar year ============
		
		*****Store predictors
		local predictors `patient_predictors_full' `practice_predictors' i.`time_adjust'
		
		melogit_model `outcome' `"`predictors'"' `"`post_nice_var'"' `"Calendar-time adjusted"' `"`outlabel'"' `"`inclusion'"' `melogit_measures' 0
	}
}

postclose `melogit_measures'

*Output postfiles to csv
use "$projectdir/output/data/melogit_summary.dta", clear
format oddsratio lower95 upper95 %9.3f
format pvalue %9.4f

export delimited using "$projectdir/output/tables/melogit_summary.csv", replace datafmt

*Run logistic regression models =====================================

**Generate temporary file to store outputs
tempname logistic_measures
postfile `logistic_measures' str80(inclusion) str80(outcome) str80(outcome_label) str80(model) str80(variable) str80(category) double n_patients n_practices df oddsratio lower95 upper95 pvalue ///
    using "$projectdir/output/data/logistic_summary.dta", replace
	
**Loop through inclusion criteria
foreach inclusion of local inclusions {

	***Import cleaned/processed cohort
	use "$projectdir/output/data/cohort_processed.dta", clear
	
	***Keep those who meet inclusion criteria
	di "Inclusion criterion: `inclusion'"
    keep if `inclusion' == 1

    ***Define outcomes and additional predictors (specific to inclusion criteria); also need to think about interactions between i.post_nice and variables - Nb. programmes about won't output factors properly
	if "`inclusion'" == "has_12m_fup" {
		local outcomes ult_12m
		*local outcomes ckd_comb_bl //this is a test
		local post_nice_var i.post_nice_diag
		local patient_predictors_full `patient_predictors_base' `post_nice_var'
		local time_adjust diagnosis_year

    }
    else if "`inclusion'" == "has_12m_fup_ult" {
		local outcomes urate_12m_ult
        *local outcomes chd_bl //this is a test
		local post_nice_var i.post_nice_ult
		local patient_predictors_full `patient_predictors_base' `post_nice_var'
		local time_adjust ult_year
    }
    else {
        di as error "No outcomes defined for `inclusion'"
        continue
    }
	
	***Loop through outcomes
	foreach outcome of local outcomes {
		
		****Store outcome variable name
		local outlabel : variable label `outcome'
		if "`outlabel'" == "" local outlabel "`outcome'"
		
		****Store predictors
		local predictors `patient_predictors_full'
		
		****Run univariable models ===============

		****Loop through predictors
		foreach predictor of local predictors {
			logistic_model `outcome' `"`predictor'"' `"`predictor'"' `"Univariable"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'
		}
		
		****Run age and sex-adjusted models ================
		
		****Loop through predictors
		foreach predictor of local predictors {
			
			local model_terms `predictor' age_decile i.sex
			if "`predictor'" == "$age_var" local model_terms $age_var i.sex
			if "`predictor'" == "i.sex" local model_terms i.sex $age_var
			
			logistic_model `outcome' `"`model_terms'"' `"`predictor'"' `"Age/sex-adjusted"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'
		}
		
		****Run multivariable model with all selected predictors ============
		
		logistic_model `outcome' `"`predictors'"' `""' `"Multivariable"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'
		
		****Sensitivity multivariable model with patient-level predictors and practice-level predictors, as well as adjusting for calendar year ============
		
		local predictors `patient_predictors_full' i.`time_adjust'
		
		logistic_model `outcome' `"`predictors'"' `"`post_nice_var'"' `"Calendar-time adjusted"' `"`outlabel'"' `"`inclusion'"' `logistic_measures'	
	}
}

*Close tempfile
postclose `logistic_measures'

*Output postfiles to csv
use "$projectdir/output/data/logistic_summary.dta", clear
format oddsratio lower95 upper95 %9.3f
format pvalue %9.4f

export delimited using "$projectdir/output/tables/logistic_summary.csv", replace datafmt

log close
