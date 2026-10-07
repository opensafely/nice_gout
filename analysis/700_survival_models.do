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
log off

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
	args model_terms focal_predictor outcome outlabel model_label run_ph
	di as txt "Model terms = `model_terms'"
	
	**Run model
	capture noisily stcox `model_terms', vce(cluster practice_id)
	
	**Skip if estimation failed
	if _rc {
		di as txt "Skipping model (estimation failure): `model_terms'"
		return scalar model_ok = 0
		exit
	}
	
	**Global proportional hazards test for primary fully adjusted model
	if "`run_ph'" == "1" {

		local ph_exposure = subinstr("`focal_predictor'", "i.", "", .)

		capture quietly estat phtest
		local ph_rc = _rc

		local ph_chi2 = .
		local ph_df = .
		local ph_pvalue = .

		if `ph_rc' == 0 {
			local ph_chi2 = r(chi2)
			local ph_df = r(df)
			local ph_pvalue = chi2tail(`ph_df', `ph_chi2')
		}

		post $cox_ph ("`outcome'") ("`ph_exposure'") ("`model_label'") (`ph_chi2') (`ph_df') (`ph_pvalue') (`ph_rc')
	}

	**Check to ensure model ran ok
	return scalar model_ok = 1
	
	**Store number of patients and practices, person-years of follow-up and n_events (all rounded) and degrees of freedom
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
    if missing(scalar(se)) continue
    if scalar(se) == 0 continue

    **Calculate HR, CI and p-values
    scalar hr = exp(scalar(b))
    scalar lo = exp(scalar(b) - invnormal(0.975)*scalar(se))
    scalar hi = exp(scalar(b) + invnormal(0.975)*scalar(se))
    scalar pv = 2*normal(-abs(scalar(b)/scalar(se)))

    local hazardratio = round(scalar(hr), 0.0001)
    local lower95 = round(scalar(lo), 0.0001)
    local upper95 = round(scalar(hi), 0.0001)
    local pvalue = round(scalar(pv), 0.0001)

    **Post model results
    post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (`hazardratio') (`lower95') (`upper95') (`pvalue')	
	}
	
	**Output model-specific risk table at specified follow-up times
	local risk_times 0 2 4 6 8 10

	**Store exposure value label
	local focal_labname : value label `focalvar'

	**Identify exposure categories included in this fitted model
	quietly levelsof `focalvar' if e(sample), local(risk_levels)

	foreach risk_level of local risk_levels {

		**Store exposure-category label
		local risk_category "`risk_level'"

		if "`focal_labname'" != "" {
			capture local risk_category : label `focal_labname' `risk_level'
			if _rc local risk_category "`risk_level'"
		}

		foreach risk_time of local risk_times {

			quietly count if e(sample) & `focalvar' == `risk_level' & _t >= `risk_time'

			local denominator_unrounded = r(N)

			**Cumulative events occurring by this time
			quietly count if e(sample) & `focalvar' == `risk_level' & _d == 1 & _t <= `risk_time'

			local events_unrounded = r(N)

			**Round both counts to nearest 5
			local denominator_rounded = round(`denominator_unrounded', 5)
			local events_rounded = round(`events_unrounded', 5)
			
			**At time zero, events are zero
			if `risk_time' == 0 {
				local events_output = 0

				**Only redact if denominator itself is less than 8
				if `denominator_unrounded' < 8 {
					local denominator_output = .
					local events_denominator "Redacted"
				}
				else {
					local denominator_output = `denominator_rounded'
					local events_denominator "0/`denominator_rounded'"
				}
			}

			**At later time points, redact if either underlying count is less than 8
			else if `events_unrounded' < 8 | `denominator_unrounded' < 8 {
				local events_output = .
				local denominator_output = .
				local events_denominator "Redacted"
			}
			else {
				local events_output = `events_rounded'
				local denominator_output = `denominator_rounded'
				local events_denominator "`events_rounded'/`denominator_rounded'"
			}

			**Post one row per category and time point
			post $cox_risk ("`outcome'") ("`outlabel'") ("`focalvar'") ("`risk_category'") ("`model_label'") (`risk_time') (`events_output') (`denominator_output') ("`events_denominator'")
		}
	}
	
	**Output absolute risk estimates at 5 years
	local absrisk_time 5

	**Store reference exposure category
	local ref_level : word 1 of `risk_levels'
	local ref_category "`ref_level'"

	if "`focal_labname'" != "" {
		capture local ref_category : label `focal_labname' `ref_level'
		if _rc local ref_category "`ref_level'"
	}

	**Crude Kaplan-Meier risk at 5 years
	if "`model_label'" == "Unadjusted" {

		tempvar km_survival km_lower km_upper

		quietly sts generate `km_survival' = s `km_lower' = lb(s) `km_upper' = ub(s) if e(sample), by(`focalvar')
		
		local reference_risk_crude = .
		
		foreach risk_level of local risk_levels {

			**Store exposure-category label
			local risk_category "`risk_level'"

			if "`focal_labname'" != "" {
				capture local risk_category : label `focal_labname' `risk_level'
				if _rc local risk_category "`risk_level'"
			}

			**Number remaining at risk at 5 years
			quietly count if e(sample) & `focalvar' == `risk_level' & _t >= `absrisk_time'
			local n_atrisk_5y = r(N)

			if `n_atrisk_5y' >= 8 {

				**Identify final observed analysis time at or before 5 years
				quietly summarize _t if e(sample) & `focalvar' == `risk_level' & _t <= `absrisk_time', meanonly
				local time_before_5 = r(max)
				
				**Kaplan-Meier risk estimate at 5 years
				quietly summarize `km_survival' if e(sample) & `focalvar' == `risk_level' & _t == `time_before_5', meanonly
	local risk_5y = 100 * (1 - r(mean))

				**Store reference crude risk
				if `risk_level' == `ref_level' {
					local reference_risk_crude = `risk_5y'
				}

				**95% confidence interval
				quietly summarize `km_upper' if e(sample) & `focalvar' == `risk_level' & _t == `time_before_5', meanonly
	local risk_5y_lower = 100 * (1 - r(mean))

				quietly summarize `km_lower' if e(sample) & `focalvar' == `risk_level' & _t == `time_before_5', meanonly
	local risk_5y_upper = 100 * (1 - r(mean))

				post $cox_absrisk ("`outcome'") ("`outlabel'") ("`focalvar'") ("`risk_category'") ("`ref_category'") ("`model_label'") ("Crude risk at 5 years") (`risk_5y') (`risk_5y_lower') (`risk_5y_upper')
				
				**Output crude absolute risk difference versus reference
				if `risk_level' != `ref_level' & !missing(`reference_risk_crude') {
					local risk_difference_crude = `risk_5y' - `reference_risk_crude'
					post $cox_absrisk ("`outcome'") ("`outlabel'") ("`focalvar'") ("`risk_category'") ("`ref_category'") ("`model_label'") ("Crude risk difference at 5 years") (`risk_difference_crude') (.) (.)
				}
			}
		else {
			post $cox_absrisk ("`outcome'") ("`outlabel'") ("`focalvar'") ("`risk_category'") ("`ref_category'") ("`model_label'") ("Crude risk at 5 years") (.) (.) (.)
		}
	}
}

	**Standardised adjusted risk at 5 years
	if "`model_label'" != "Unadjusted" {

		tempvar baseline_survival focal_original xb risk_cf

		**Store original exposure values
		gen double `focal_original' = `focalvar'

		**Estimate baseline survivor function
		quietly predict double `baseline_survival' if e(sample), basesurv

		**Check sufficient follow-up at 5 years
		quietly count if e(sample) & _t >= `absrisk_time'
		local n_atrisk_5y = r(N)

		if `n_atrisk_5y' >= 8 {

			**Baseline survival at 5 years
			quietly summarize `baseline_survival' if e(sample) & _t <= `absrisk_time', meanonly

			if r(N) > 0 {
				local s0_5 = r(min)
			}
			else {
				local s0_5 = 1
			}

			local reference_risk = .

			foreach risk_level of local risk_levels {

				**Store exposure-category label
				local risk_category "`risk_level'"

				if "`focal_labname'" != "" {
					capture local risk_category : label `focal_labname' `risk_level'
					if _rc local risk_category "`risk_level'"
				}

				**Set everyone to this exposure category
				replace `focalvar' = `risk_level' if e(sample)

				**Calculate predicted linear predictor
				capture drop `xb'
				quietly predict double `xb' if e(sample), xb

				**Calculate individual predicted 5-year risks
				capture drop `risk_cf'
				gen double `risk_cf' = 1 - (`s0_5'^exp(`xb')) if e(sample)

				**Average across estimation population
				quietly summarize `risk_cf' if e(sample), meanonly
				local adjusted_risk = r(mean)
				local adjusted_risk_pct = 100 * `adjusted_risk'

				**Store reference risk
				if `risk_level' == `ref_level' {
					local reference_risk = `adjusted_risk'
				}

				**Output adjusted absolute risk
				post $cox_absrisk ("`outcome'") ("`outlabel'") ("`focalvar'") ("`risk_category'") ("`ref_category'") ("`model_label'") ("Adjusted risk at 5 years") (`adjusted_risk_pct') (.) (.)

				**Output absolute risk difference versus reference
				if `risk_level' != `ref_level' & !missing(`reference_risk') {

					local risk_difference_pct = 100 * (`adjusted_risk' - `reference_risk')

					post $cox_absrisk ("`outcome'") ("`outlabel'") ("`focalvar'") ("`risk_category'") ("`ref_category'") ("`model_label'") ("Adjusted risk difference at 5 years") (`risk_difference_pct') (.) (.)
				}

				**Restore original exposure
				replace `focalvar' = `focal_original' if e(sample)
			}
		}

		drop `focal_original'
	}
end

*Multiple imputation models =============
capture program drop cox_model_mi

program define cox_model_mi, rclass

	**Model arguments
	args model_terms focal_predictor outcome outlabel model_label
	
	di as txt "MI model terms = `model_terms'"
	
	**Run multiply imputed Cox model
	capture noisily mi estimate, post noisily: stcox `model_terms', vce(cluster practice_id)
	
	local mi_rc = _rc

	if `mi_rc' {
		di as error "MI Cox estimation failed; return code `mi_rc'"
		return scalar model_ok = 0
		exit
	}

	di as result "MI Cox estimation completed"
	
	**Check model ran
	return scalar model_ok = 1
	
	**MI model sample descriptors, rounded to nearest 5
	local n_patients  = round(e(N_mi), 5)
	local n_practices = round(e(N_clust), 5)
	local n_events   = round(e(N_fail), 5)
	local person_years = round(e(risk), 5)

	**Model degrees of freedom, consistent with ordinary Cox rows
	local df = e(df_r_mi)
	
	**Strip factor prefix from focal predictor
	local focalvar "`focal_predictor'"
	local focalvar = subinstr("`focalvar'", "i.", "", .)
	local focalvar = subinstr("`focalvar'", "c.", "", .)
	
	**Store pooled MI estimates
	matrix B = e(b)
	matrix V = e(V)
	
	capture matrix drop D
	capture matrix D = e(df_mi)
	
	local cnames : colnames B
	
	**Cycle through column names
	foreach term of local cnames {

		**Skip intercepts
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
			local varname "`=regexs(2)'"
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
			local varname "`=regexs(2)'"
			local base = 1
		}

		**Handle regular factor terms
		else if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname "`=regexs(3)'"
		}

		**Store factor level label
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

		**Restrict output to focal predictor
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

		**Post base terms
		if `base' == 1 {
			post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (1) (.) (.) (.)
			continue
		}
		
		**Identify coefficient position
		local col = colnumb(B, "`term'")
		if missing(`col') continue

		**Pooled coefficient and SE
		scalar b = B[1, `col']
		scalar se = sqrt(V[`col', `col'])

		if missing(scalar(se)) continue
		if scalar(se) == 0 continue

		**Parameter-specific MI degrees of freedom
		scalar df_mi = .
		capture scalar df_mi = D[1, `col']

		**Calculate confidence limits and p-value
		if missing(scalar(df_mi)) {
			scalar crit = invnormal(0.975)
			scalar pv = 2*normal(-abs(scalar(b)/scalar(se)))
		}
		else {
			scalar crit = invttail(scalar(df_mi), 0.025)
			scalar pv = 2*ttail(scalar(df_mi), abs(scalar(b)/scalar(se)))
		}

		scalar hr = exp(scalar(b))
		scalar lo = exp(scalar(b) - scalar(crit)*scalar(se))
		scalar hi = exp(scalar(b) + scalar(crit)*scalar(se))

		local hazardratio = round(scalar(hr), 0.0001)
		local lower95 = round(scalar(lo), 0.0001)
		local upper95 = round(scalar(hi), 0.0001)
		local pvalue = round(scalar(pv), 0.0001)

		**Post pooled MI results to same Cox output
		post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (`hazardratio') (`lower95') (`upper95') (`pvalue')
	}
	
end

*Competing risk models (stcrprep) =========================*/
capture program drop competing_risk_model_stc

program define competing_risk_model_stc, rclass

	**Model arguments
    args model_terms focal_predictor outcome outlabel model_label compete_var

    **Save original patient-level dataset and survival settings
    tempfile fg_original
    quietly save "`fg_original'"

    log on
    di as text "Fine-Gray START: `c(current_date)' `c(current_time)'"

    **Prepare data and fit model; capture errors so data can be restored
    capture noisily {

        which stcrprep

        **Underlying variables for current main-effect model (model_terms must contain covariates only, without an if qualifier)
        local fg_vars = subinstr("`model_terms'", "i.", "", .)
        local fg_vars = subinstr("`fg_vars'", "c.", "", .)
        local fg_vars : list uniq fg_vars

        **Restrict to the complete-case analysis population
        keep if _st == 1

        foreach v of local fg_vars {
            quietly drop if missing(`v')
        }

        quietly drop if missing(practice_id)

        **Current implementation assumes one record per patient, with follow-up starting at landmark (analysis time zero)
        isid patient_id
        assert _t0 == 0
        assert !missing(`compete_var')
        assert _d == 0 if `compete_var' == 1

        **Store actual sample descriptors before expansion
        quietly count
        local n_patients = round(r(N), 5)

        tempvar fg_practice_tag
        egen `fg_practice_tag' = tag(practice_id)
        quietly count if `fg_practice_tag' == 1
        local n_practices = round(r(N), 5)

        quietly count if _d == 1
        local n_events = round(r(N), 5)

        quietly summarize _t, meanonly
        local person_years = round(r(sum), 5)

        **Event type: 0=censored, 1=CKD, 2=competing death
        tempvar fg_status fg_time
        gen byte `fg_status' = 0
        replace `fg_status' = 1 if _d == 1
        replace `fg_status' = 2 if `compete_var' == 1

        gen double `fg_time' = _t

        **Initial stset must identify all event types and patient ID
        quietly stset `fg_time', failure(`fg_status' == 1 2) id(patient_id)

        **Expand data for CKD only and calculate censoring weights
        stcrprep, events(`fg_status') trans(1) keep(`fg_vars' practice_id)

        **Identify the event of interest in expanded data
        tempvar fg_event
        gen byte `fg_event' = (failcode == `fg_status')

        **Weighted intervals: retain years as the time unit
        quietly stset tstop [pw=weight_c], enter(time tstart) failure(`fg_event' == 1)

        di as text "Weighted Cox START: `c(current_date)' `c(current_time)'"

        **This weighted Cox fit estimates Fine-Gray subhazard ratios
        stcox `model_terms', vce(cluster practice_id) breslow log iterate(20)

        **Do not post results from an unconverged model
        if e(converged) != 1 {
            error 430
        }

        local df = e(df_m)
    }

    local fg_rc = _rc

    di as text "Fine-Gray END: `c(current_date)' `c(current_time)'"
    di as text "Fine-Gray return code: `fg_rc'"
    log off

    if `fg_rc' {
        quietly use "`fg_original'", clear
        return scalar model_ok = 0
        exit
    }

	**Strip factor prefix from focal predictor
	local focalvar "`focal_predictor'"
	local focalvar = subinstr("`focalvar'", "i.", "", .)
	local focalvar = subinstr("`focalvar'", "c.", "", .)

	**Store outputs from model
	matrix B = e(b)
	local cnames : colnames B

	**Cycle through column names
	foreach term of local cnames {

		**Skip intercepts
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
			local varname "`=regexs(2)'"
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
			local varname "`=regexs(2)'"
			local base = 1
		}

		**Handle regular factor terms
		else if regexm("`term'", "^([0-9]+)([a-z]*)\.(.+)$") {
			local levelnum "`=regexs(1)'"
			local varname "`=regexs(3)'"
		}

		**Store factor level label
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

		**Restrict output to focal predictor
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

		**Post reference category
		if `base' == 1 {
			post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (1) (.) (.) (.)
			continue
		}

		**Extract coefficient and SE
		capture scalar b = _b[`term']
		if _rc continue

		capture scalar se = _se[`term']
		if _rc continue

		if missing(scalar(se)) continue
		if scalar(se) == 0 continue

		**Calculate subhazard ratio, CI and p-value
		scalar shr = exp(scalar(b))
		scalar lo = exp(scalar(b) - invnormal(0.975)*scalar(se))
		scalar hi = exp(scalar(b) + invnormal(0.975)*scalar(se))
		scalar pv = 2*normal(-abs(scalar(b)/scalar(se)))

		local subhazardratio = round(scalar(shr), 0.0001)
		local lower95 = round(scalar(lo), 0.0001)
		local upper95 = round(scalar(hi), 0.0001)
		local pvalue = round(scalar(pv), 0.0001)

		**Post model results
		post $cox_measures ("`outcome'") ("`outlabel'") ("`varlabel'") ("`category'") ("`model_label'") (`n_patients') (`n_practices') (`n_events') (`person_years') (`df') (`subhazardratio') (`lower95') (`upper95') (`pvalue')
	}
	
	**Restore patient-level data and original survival settings
    quietly use "`fg_original'", clear

    return scalar model_ok = 1
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
egen censor_date_death = rowmin(`study_end_date' `dereg_date')
format censor_date_death %td

**Primary exposure variable
local exposure_complete_360 urate_12m_ult //urate checked and target attained vs. not attained within 12 months of ULT initiation (coded as 1/0/missing)

**Sensitivity exposure variables
local exposure_nomiss_360 urate_12m_ult_recode //urate recoded as not attained if urate not checked (coded as 1/0)
local exposure_sens_misscat urate_12m_ult_cat //separate category coded if urate not checked (1/0/9)
local exposure_sens_300 urate_300_12m_ult
local exposure_sens_300_360 urate_targets_12m_ult

**Define exposure list to loop through
local primary_exposure `exposure_complete_360'
local secondary_exposures `exposure_nomiss_360'  `exposure_sens_misscat' `exposure_sens_300' `exposure_sens_300_360'
*local exposures `primary_exposure' `secondary_exposures'
local exposures `primary_exposure'

**Primary outcome
gen sec_ckd_egfr_land_date = second_egfr_ckd_date if (second_egfr_ckd_date > `landmark_date') & second_egfr_ckd_date !=. & `landmark_date' !=.
format sec_ckd_egfr_land_date %td
label var sec_ckd_egfr_land_date "Incident CKD by two eGFRs <60 after ULT landmark"

**Secondary/sensitivity outcomes
gen first_ckd_egfr_land_date = first_egfr_ckd_date if (first_egfr_ckd_date > `landmark_date') & first_egfr_ckd_date !=. & `landmark_date' !=.
format first_ckd_egfr_land_date %td
label var first_ckd_egfr_land_date "Incident CKD by one eGFR <60 after ULT landmark"

gen first_ckd_code_land_date = ckd_date if (ckd_date > `landmark_date') & ckd_date !=. & `landmark_date' !=.
format first_ckd_code_land_date %td
label var first_ckd_code_land_date "Incident CKD by one CKD code after ULT landmark"

**Death outcome
gen death_land_date = date_of_death if date_of_death > `landmark_date' & !missing(date_of_death) & !missing(`landmark_date')
format death_land_date %td
label var death_land_date "All-cause mortality after ULT landmark"

**Hernia (negative control) outcome
gen hernia_land_date = hernia_date if hernia_date > `landmark_date' & !missing(hernia_date) & !missing(`landmark_date')
format hernia_land_date %td
label var hernia_land_date "Incident inguinal hernia after ULT landmark"

**Define outcome list to loop through
local primary_outcome sec_ckd_egfr_land_date
local secondary_outcomes first_ckd_egfr_land_date first_ckd_code_land_date death_land_date hernia_land_date
*local outcomes `primary_outcome' `secondary_outcomes'
local outcomes `primary_outcome'

**Outcome status at baseline/landmark variables
local outcome_free_baseline ckd_free_ult //CKD, defined using single eGFR <60 or CKD code at or before ULT initiation date
local outcome_free_landmark ckd_free_landmark //CKD, defined using single eGFR <60 or CKD code at or before ULT initiation date + 12 months

**Define patient-level predictors
local patient_predictors_core ///
    age_land_decile i.sex i.imd i.ethnicity c.bmi_value i.smoke i.diabetes_land i.heart_failure_land i.chd_land i.cva_land i.hypertension_land i.alcohol_land i.diuretic_land i.sglt2_land i.ace_arb_land
    *rheum_appt_n_12m hosp_n_12m creatinine_n_12m
	
local patient_predictors_extra ///
	urate_before_ult_value egfr_before_ult_value
	
/***For BMI, combine obesity categories due to sparse data for Obese III
replace bmicat = 4 if inlist(bmicat, 5, 6)
label define bmicat_lab 1 "Underweight" 2 "Normal weight" 3 "Overweight" 4 "Obese" 9 "Missing", replace
label values bmicat bmicat_lab
*/

***Recode categorised missing covariates as missing	
foreach var in imd ethnicity bmicat smoke {
    replace `var' = . if `var' == 9
}
	
*Run landmark Cox models =======================================================

**Generate temporary file to store outputs
tempname cox_measures
postfile `cox_measures' str150(outcome) str150(outcome_label) str150(exposure) str150(exposure_category) str80(model) double n_patients n_practices n_events person_years df hazardratio lower95 upper95 pvalue ///
    using "$projectdir/output/data/landmark_cox_summary.dta", replace
	
global cox_measures `cox_measures'

**Generate temporary file to store model-specific risk tables
tempname cox_risk
postfile `cox_risk' str150(outcome) str150(outcome_label) str150(exposure) str150(exposure_category) str80(model) double time_years events denominator str30(events_denominator) ///
    using "$projectdir/output/data/landmark_cox_risk_table.dta", replace
	
global cox_risk `cox_risk'

tempname cox_absrisk
postfile `cox_absrisk' str150(outcome) str150(outcome_label) str150(exposure) str150(exposure_category) str150(reference_category) str80(model) str50(measure) double estimate_pct lower95_pct upper95_pct using "$projectdir/output/data/landmark_cox_absrisk.dta", replace	

global cox_absrisk `cox_absrisk'

**Global proportional hazards test results
tempname cox_ph

postfile `cox_ph' str150 outcome str150 exposure str80 model double chi2 df pvalue int return_code ///
    using "$projectdir/output/data/landmark_cox_ph.dta", replace

global cox_ph `cox_ph'

capture stset, clear

preserve

**Cohort of interest (exposure status set below)
keep if !missing(`landmark_date') & !missing(censor_date) & (censor_date > `landmark_date') //landmark date present and before censor date
keep if `outcome_free_baseline' ==1 //outcome not present before cohort entry
keep if `outcome_free_landmark' ==1 //outcome not present before landmark

**Save common cohort so each outcome starts with the same patients
tempfile landmark_cohort
quietly save `landmark_cohort'

**Loop through outcomes
foreach outcome of local outcomes {
	
	quietly use `landmark_cohort', clear
	
	**Exclude previous hernia only for the incident hernia analysis
    if "`outcome'" == "hernia_land_date" {
        keep if hernia_free_landmark == 1
    }
	
	di as txt "Outcome = `outcome'"
			
	***Store outcome variable name
	local outlabel : variable label `outcome'
	if "`outlabel'" == "" local outlabel "`outcome'"
	
	***Clear previous values
	capture drop stop_date fail
	capture stset, clear
	
	***Assign fail and stop dates
	if "`outcome'" == "death_land_date" {
		gen stop_date = censor_date_death
		replace stop_date = `outcome' if !missing(`outcome') & (`outcome' <= censor_date_death)
	}
	else {
		gen stop_date = censor_date
		replace stop_date = `outcome' if !missing(`outcome') & (`outcome' <= censor_date)
	}
	
	format stop_date %td
	
	gen fail = 0
	replace fail = 1 if !missing(`outcome') & `outcome' == stop_date
	
	***Define death as competing event (for sensitivity analyses)
	capture drop death_compete
	gen death_compete = 0
	replace death_compete = 1 if !missing(`death_date') & `death_date' <= `study_end_date' & (`death_date' <= `dereg_date' | missing(`dereg_date')) & (missing(`outcome') | `death_date' < `outcome')

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
		
		****Run multivariable model with baseline urate and eGFR (values closest to before ULT initiation, but within 12m); also output PH test for this if it is the primary model
		local model_terms i.`exposure' `patient_predictors_core' `patient_predictors_extra' if !missing(`exposure')
		local run_ph = ("`outcome'" == "`primary_outcome'" & "`exposure'" == "`primary_exposure'")
		cox_model `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"Multivariable extra"' "`run_ph'"
		
	    ****Run Fine-Gray competing-risk models (stcprep) for primary model
		if "`outcome'" == "`primary_outcome'" & "`exposure'" == "`primary_exposure'" {
			local model_terms i.`exposure' `patient_predictors_core' `patient_predictors_extra'
				competing_risk_model_stc `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"Fine-Gray extra"' `"death_compete"'
		}
		
		****Output KM and loglog plots
		
		*****Store labels for graph
		levelsof `exposure' if !missing(`exposure') & _st == 1, local(levels)

		local colours "emerald orange red blue dkgreen cranberry navy maroon teal sienna purple"
		
		local legtitle : variable label `exposure'
		if "`legtitle'" == "" local legtitle "`exposure'"

		if inlist("`exposure'", "urate_12m_ult", "urate_12m_ult_cat", "urate_12m_ult_recode", "urate_300_12m_ult", "urate_targets_12m_ult") {
			local legtitle "Urate target attained"
		}

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
		
		*****Determine latest non-redacted KM time point
		local risk_times 0 1 2 3 4 5 6 7 8 9 10
		local km_tmax 0

		foreach risk_time of local risk_times {
			local time_ok 1
			
			foreach risk_level of local levels {
				quietly count if _st == 1 & !missing(`exposure') & `exposure' == `risk_level' & _t >= `risk_time'
				if r(N) < 8 {
					local time_ok 0
				}
			}

			**Update maximum only when all exposure groups are non-redacted
			if `time_ok' == 1 {
				local km_tmax `risk_time'
			}
		}

		di as text "KM/log-log plots for `exposure' truncated at `km_tmax' years"
		
		*****Naming of graphs
		local graphstub = substr("`exposure'_`outcome'", 1, 25)
		local kmname "km_`graphstub'"
		local loglogname "ll_`graphstub'"
		
		*****Survival plot
		if `km_tmax' > 0 {
		
			capture noisily sts graph if !missing(`exposure') & _st==1, by(`exposure') survival tmax(`km_tmax') `km_plotopts' ytitle("Survival probability", size(medsmall)) ylabel(, nogrid labsize(small)) xtitle("Years from landmark", size(medsmall) margin(medsmall)) xlabel(0(1)`km_tmax', nogrid labsize(small)) title("", size(medium) margin(b=2)) legend(order(`legorder') title("`legtitle'", size(small) margin(b=1))) xsize(16) ysize(9) name(`kmname', replace) saving("$projectdir/output/figures/km_`exposure'_`outcome'.gph", replace)
			
			if _rc == 0 {
				capture graph export "$projectdir/output/figures/km_`exposure'_`outcome'.$img", replace
				
				if _rc == 0 {
					local ++n_km_graphs
				}
			}
			
			*****Log-log plot truncated at latest non-redacted time

			**Store truncated follow-up in years
			tempvar stop_truncated fail_truncated

			gen double `stop_truncated' = min(_t, `km_tmax')
			gen byte `fail_truncated' = _d == 1 & _t <= `km_tmax'

			**Temporarily reset survival data using truncated follow-up
			quietly stset `stop_truncated', failure(`fail_truncated' == 1)

			**Default: no log-log graph produced
			local loglog_graph_ok = 0

			**Check for events within truncated follow-up
			quietly summarize _t if !missing(`exposure') & _st == 1 & _d == 1 & _t > 0, meanonly

			if r(N) > 0 {

				local log_xmin = floor(ln(r(min)))

				quietly summarize _t if !missing(`exposure') & _st == 1 & _t > 0, meanonly

				local log_xmax = ceil(ln(r(max)))

				capture noisily stphplot if !missing(`exposure') & _st == 1, by(`exposure') `loglog_plotopts' ytitle("-log{-log(Survival probability)}", size(medsmall)) ylabel(, nogrid labsize(small)) xtitle("log(Time)", size(medsmall) margin(medsmall)) xscale(range(`log_xmin' `log_xmax')) xlabel(`log_xmin'(1)`log_xmax', nogrid labsize(small)) title("") legend(order(`legorder') title("`legtitle'", size(small) margin(b=1))) xsize(16) ysize(9) name(`loglogname', replace) saving("$projectdir/output/figures/loglog_`exposure'_`outcome'.gph", replace)

				local loglog_graph_ok = (_rc == 0)
			}

			**Restore original survival settings regardless of graph success
			stset stop_date, origin(time `landmark_date') scale(365.25) failure(fail == 1)

			drop `stop_truncated' `fail_truncated'

			**Export only if the graph was created successfully
			if `loglog_graph_ok' {

				capture graph export "$projectdir/output/figures/loglog_`exposure'_`outcome'.$img", name(`loglogname') replace

				if _rc == 0 {
					local ++n_loglog_graphs
				}
			}
		}
		else {
			di as text "No non-redacted follow-up beyond time zero; skipping KM and log-log plots."
		}

		****Run multiply imputed models for primary model
		if "`outcome'" == "`primary_outcome'" & "`exposure'" == "`primary_exposure'" {
			
			**Temporarily save current analysis dataset
			tempfile pre_mi
			quietly save `pre_mi'
			
			**Restrict to current exposure analysis population
			keep if !missing(`exposure')

			**Generate Nelson-Aalen cumulative hazard for imputation model
			capture drop na_hazard
			sts generate na_hazard = na
			
			**Ensure consistent patient order
			isid patient_id
			sort patient_id

			**Convert to MI data
			mi set mlong

			**Register variables to be imputed
			mi register imputed imd ethnicity bmi_value smoke urate_before_ult_value egfr_before_ult_value

			**Register regular variables
			mi register regular `exposure' `landmark_date' age_land_decile sex diabetes_land heart_failure_land chd_land cva_land hypertension_land alcohol_land diuretic_land sglt2_land ace_arb_land stop_date fail na_hazard practice_id
			
			log on
			di as text "MI START: `c(current_date)' `c(current_time)'"
			
			**Multiple imputation by chained equations // Limit categorical regression iterations; show ethnicity diagnostics
			capture noisily mi impute chained (ologit, iterate(20)) imd (mlogit, iterate(20) augment noisily) ethnicity (mlogit, iterate(20) augment) smoke (pmm, knn(5)) bmi_value urate_before_ult_value egfr_before_ult_value = i.`exposure' age_land_decile i.sex i.diabetes_land i.heart_failure_land i.chd_land i.cva_land i.hypertension_land i.alcohol_land i.diuretic_land i.sglt2_land i.ace_arb_land fail na_hazard, add(2) burnin(2) rseed(12345) showevery(1)
			
			**Skip MI models if imputation fails
			local mi_rc = _rc
			
			di as text "MI Imputation END: `c(current_date)' `c(current_time)'"
			tab ethnicity

			if `mi_rc' {
				di as error "MI imputation failed: `outcome' / `exposure'; return code `mi_rc'"
				quietly use `pre_mi', clear	
				log off
				continue
			}

			di as result "MI imputation completed: `outcome' / `exposure'"
			mi describe
			
			**Set survival data for MI analysis
			mi stset stop_date, origin(time `landmark_date') scale(365.25) failure(fail == 1)
			
			**MI multivariable model including baseline urate and eGFR
			local model_terms i.`exposure' `patient_predictors_core' `patient_predictors_extra'
			cox_model_mi `"`model_terms'"' `"i.`exposure'"' `"`outcome'"' `"`outlabel'"' `"MI multivariable extra"'
						
			**Restore dataset before MI
			quietly use `pre_mi', clear
			
			di as text "MI END: `c(current_date)' `c(current_time)'"
			log off
		}
	}
}
	
restore

*Close tempfile
postclose $cox_measures
postclose $cox_risk
postclose $cox_absrisk
postclose $cox_ph

*Output postfiles to csv - with failsafes
capture use "$projectdir/output/data/landmark_cox_summary.dta", clear

**Create empty risk table if no results were produced
if _rc {
	clear
	set obs 0

	foreach v in outcome outcome_label exposure exposure_category {
		gen str150 `v' = ""
	}
	gen str80 model = ""

	foreach v in n_patients n_practices n_events person_years df hazardratio lower95 upper95 pvalue {
		gen double `v' = .
	}
}

format hazardratio lower95 upper95 pvalue %12.4f

export delimited using "$projectdir/output/tables/landmark_cox_summary.csv", replace datafmt

*Output 5-year absolute risks
capture use "$projectdir/output/data/landmark_cox_absrisk.dta", clear

if _rc {

	**Create empty absolute risk table if no results were produced
	clear
	set obs 0
	gen str150 outcome = ""
	gen str150 outcome_label = ""
	gen str150 exposure = ""
	gen str150 exposure_category = ""
	gen str150 reference_category = ""
	gen str80 model = ""
	gen str50 measure = ""
	gen double estimate_pct = .
	gen double lower95_pct = .
	gen double upper95_pct = .
}
else {
	format estimate_pct lower95_pct upper95_pct %9.2f
	sort outcome exposure model measure exposure_category
}

export delimited using "$projectdir/output/tables/landmark_cox_absrisk.csv", replace datafmt

*Output Cox model results
capture use "$projectdir/output/data/landmark_cox_risk_table.dta", clear

if _rc {

	**Create empty risk table if no results were produced
	clear
	set obs 0

	gen str150 outcome = ""
	gen str150 outcome_label = ""
	gen str150 exposure = ""
	gen str150 exposure_category = ""
	gen str80 model = ""

	foreach t in 0 2 4 6 8 10 {
		gen double events_`t'y = .
		gen double atrisk_`t'y = .
	}
}
else {

	**Reshape risk-table results from long to wide
	keep outcome outcome_label exposure exposure_category model time_years events denominator

	***Only reshape if at least one valid time point exists
	quietly count if !missing(time_years)

	if r(N) > 0 {

		drop if missing(time_years)

		reshape wide events denominator, ///
			i(outcome outcome_label exposure exposure_category model) ///
			j(time_years)

		foreach t in 0 2 4 6 8 10 {
			capture rename events`t' events_`t'y
			capture rename denominator`t' atrisk_`t'y
		}

		sort outcome exposure model exposure_category
	}
	else {
		di as text "No valid risk-table times found; skipping risk-table output."
	}
}

export delimited using "$projectdir/output/tables/landmark_cox_risk_table.csv", replace

**Export PH results
use "$projectdir/output/data/landmark_cox_ph.dta", clear

format chi2 pvalue %12.4f
format df return_code %9.0f

export delimited using "$projectdir/output/tables/landmark_cox_ph.csv", datafmt replace

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
