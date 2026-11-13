version 16

/*==============================================================================
DO FILE NAME:           Rounded and redacted consultation summary tables
PROJECT:                OpenSAFELY NICE
AUTHOR:                 M Russell
DATASETS USED:          cohort_consults_processed.dta
DATASETS CREATED:       Practice-summary CSV
USER-INSTALLED ADO:     None
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
log using "$logdir/consults_tables.log", replace
log off
set more off
set type double

*Read variables passed from yaml
global arglist disease intervalstart_date
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
    global intervalstart_date "$studyend_date"
}

*Blood test and medication defined in the consultation dataset
global blood "urate"
global medication "ult"
local Medication = upper("$medication")
local Blood = upper("$blood")

*Inclusive end of the 12-calendar-month consultation interval
local interval_start = date("$intervalstart_date", "YMD")
local interval_end_month = mofd(`interval_start') + 12
local interval_end = dofm(`interval_end_month') + min(day(`interval_start'), day(dofm(`interval_end_month' + 1) - 1)) - 2
global intervalend_date = string(`interval_end', "%tdCCYY-NN-DD")

*Function to summarise practices, round N and redact outputs ==============*

capture program drop rounded_practice_summary
program define rounded_practice_summary
    syntax varlist(min=1 max=1 numeric), denominator(varname) outfile(string) [COUNTS ICC LABEL(string) COMPLEMENTS(varlist)]

    local outcome `varlist'
    local outcome_desc : variable label `outcome'
    if `"`label'"' != "" local outcome_desc `"`label'"'

    preserve

        **A complement of 1-7 patients suppresses this whole output row, including N.
        local small_complement = 0
        foreach flag of local complements {
            quietly count if `flag' == 1 & !missing(consult_practice_id)
            if inrange(r(N), 1, 7) local small_complement = 1
        }

        **Apply the patient denominator for this row
        keep if `denominator' == 1 & !missing(consult_practice_id, `outcome')
        quietly count
        local patient_n = r(N)
        quietly count if `outcome' == 1
        local events = r(N)
        local nonevents = `patient_n' - `events'

        **Default release values; no unrounded data are saved to the results file
        local n_text "0"
        local mean_text "Redacted"
        local sd_text "Redacted"
        local median_text "Redacted"
        local p25_text "Redacted"
        local p75_text "Redacted"
        local icc_text "Not required"
        local se_text "Not required"
        if "`icc'" != "" {
            local icc_text "Redacted"
            local se_text "Redacted"
        }

        if `patient_n' > 0 {

            **Count contributing practices before rounding
            egen practice_group = group(consult_practice_id)
            egen practice_tag = tag(practice_group)
            quietly count if practice_tag == 1
            local practice_n = r(N)

            local n_text "<8"
            if `practice_n' >= 8 local n_text = strtrim(string(round(`practice_n', 5), "%12.0f"))

            **Primary suppression: require sufficient patients and practices
            local release = `patient_n' >= 8 & `practice_n' >= 8

            **For proportions, protect both numerator and complementary outcome count
            **COUNTS summarises practice patient numbers instead of a binary proportion.
            if "`counts'" == "" & (`events' < 8 | `nonevents' < 8) local release = 0
            if `small_complement' local release = 0
            if !`release' & `practice_n' >= 8 local n_text "Redacted"

            **Optional observed-scale ICC and asymptotic SE, using PATIENT outcomes
            **Save returned values before tabstat overwrites r().
            if "`icc'" != "" & `release' {
                local icc_text "Not estimable"
                local se_text "Not estimable"
                capture quietly loneway `outcome' practice_group
                if !_rc {
                    local rho = r(rho)
                    local rho_se = r(se)
                    if !missing(`rho') local icc_text = strtrim(string(round(`rho', .001), "%9.3f"))
                    if !missing(`rho_se') local se_text = strtrim(string(round(`rho_se', .001), "%9.3f"))
                }
            }

            **Collapse to practice counts and proportions using exact internal data
            collapse (count) denominator_n=`outcome' (sum) numerator_n=`outcome', by(consult_practice_id)
            gen practice_value = numerator_n / denominator_n
            if "`counts'" != "" replace practice_value = denominator_n

            **Summarise across practices, giving each practice equal weight
            if `release' {
                quietly tabstat practice_value, stats(n mean sd median p25 p75) save
                matrix table = r(StatTotal)'
                local mean_text = strtrim(string(round(table[1,2], .001), "%9.3f"))
                local sd_text = strtrim(string(round(table[1,3], .001), "%9.3f"))
                local median_text = strtrim(string(round(table[1,4], .001), "%9.3f"))
                local p25_text = strtrim(string(round(table[1,5], .001), "%9.3f"))
                local p75_text = strtrim(string(round(table[1,6], .001), "%9.3f"))
            }
        }

        if `small_complement' local n_text "Redacted"

        **Create one rounded/redacted summary row
        clear
        set obs 1
        gen str10 interval_start = "$intervalstart_date"
        gen str10 interval_end = "$intervalend_date"
        gen str160 variable = `"`outcome_desc'"'
        gen str12 n_practices = "`n_text'"
        gen str20 mean = "`mean_text'"
        gen str20 sd = "`sd_text'"
        gen str20 median = "`median_text'"
        gen str20 p25 = "`p25_text'"
        gen str20 p75 = "`p75_text'"
        gen str20 icc = "`icc_text'"
        gen str20 se_icc = "`se_text'"

        **Append earlier summary rows; do not save any patient/practice records
        capture confirm file `"`outfile'"'
        if !_rc append using `"`outfile'"'
        save `"`outfile'"', replace

    restore
end

*Load processed patient dataset =========================================*

use "$projectdir/output/data/cohort_consults_processed.dta", clear

*Define patient outcomes and denominators ================================*

gen all_consults = 1
gen baseline_test = had_${blood}_bl == 1
gen baseline_below = ${blood}_bl_cat == 1
gen full_followup = has_12m_fup == 1
gen eligible_followup = eligible == 1 & full_followup == 1
gen eligible_test_followup = eligible_followup == 1 & ${blood}_within_12m == 1

**The prescription window was already applied in the dataset definition
gen medication_within_12m = !missing(${medication}_after_date)

**New prescribing requires no prescription in the baseline six-month window
**Its denominator below remains ALL eligible patients with follow-up
gen medication_new_12m = ${medication}_bl == 0 & medication_within_12m == 1

**Any prescription: baseline six-month window or the twelve months after consultation
gen medication_any_12m = ${medication}_bl == 1 | medication_within_12m == 1

*Define differences between overlapping table populations ================*

**These patient flags allow the same function to protect comparisons between rows.
gen no_baseline_test = baseline_test == 0
gen eligible_no_test = eligible_followup == 1 & eligible_test_followup == 0
gen not_eligible_followup = eligible_followup == 0
gen controlled_followup = full_followup == 1 & eligible == 0

*Generate the practice-summary table =====================================*

tempfile results

**Patient numbers, not repeated consultation episodes
rounded_practice_summary all_consults, denominator(all_consults) outfile("`results'") counts label("Number of patients consulting per practice")

**Baseline and follow-up proportions
rounded_practice_summary full_followup, denominator(all_consults) outfile("`results'") label("Proportion with 12m+ follow-up")
rounded_practice_summary ${medication}_bl, denominator(all_consults) outfile("`results'") label("Baseline prescribed `Medication'")
rounded_practice_summary baseline_test, denominator(all_consults) outfile("`results'") label("Baseline recorded ${blood} (24-month window)")
rounded_practice_summary baseline_below, denominator(all_consults) outfile("`results'") label("Baseline ${blood} <360 (all cases)")
rounded_practice_summary baseline_below, denominator(baseline_test) outfile("`results'") label("Baseline ${blood} <360 (with test)") complements(no_baseline_test)
rounded_practice_summary controlled_bl, denominator(all_consults) outfile("`results'") label("Baseline ${blood} <360 and on `Medication'")
rounded_practice_summary eligible, denominator(all_consults) outfile("`results'") label("Eligible proportion")

**Actual practice cohort size for the eligible-with-follow-up outcome analyses
rounded_practice_summary all_consults, denominator(eligible_followup) outfile("`results'") counts label("Number eligible with 12m+ follow-up per practice") complements(not_eligible_followup controlled_followup)

**Outcome rows: ICC and SE only where requested
rounded_practice_summary medication_new_12m, denominator(eligible_followup) outfile("`results'") icc label("`Medication' initiated within 12m (eligible with follow-up)")
rounded_practice_summary ${blood}_within_12m, denominator(eligible_followup) outfile("`results'") label("`Blood' recorded within 12m (eligible with follow-up)")
rounded_practice_summary ${blood}_target_12m, denominator(eligible_followup) outfile("`results'") icc label("`Blood' <360 on `Medication' within 12m (eligible with follow-up, all)")
rounded_practice_summary ${blood}_target_12m, denominator(eligible_test_followup) outfile("`results'") icc label("`Blood' <360 on `Medication' within 12m (eligible with follow-up and test)") complements(eligible_no_test)
rounded_practice_summary medication_any_12m, denominator(full_followup) outfile("`results'") label("`Medication' before or within 12m after consultation (all cases with follow-up)")

*Export the rounded and redacted aggregate table ==========================*

use "`results'", clear

**Appending places the latest row first; restore the order of calls above
gen row_order = -_n
sort row_order
drop row_order

export delimited using "$projectdir/output/tables/consults_table.csv", replace
log close
