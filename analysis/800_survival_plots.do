version 16

/*==============================================================================
DO FILE NAME:			Survival plots
PROJECT:				OpenSAFELY NICE 
AUTHOR:					M Russell								
DATASETS USED:			Survival tables (rounded)
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
log using "$logdir/survival_plots.log", replace

*Set Ado file path
adopath + "$projectdir/analysis/extra_ados"

set type double

set scheme plotplainblind

local n_km_graphs = 0
local n_loglog_graphs = 0

capture erase "$projectdir/output/figures/km_recon_no_outputs.$img"
capture erase "$projectdir/output/figures/loglog_recon_no_outputs.$img"

*Generate reconstructed survival plots from rounded counts ==========================
import delimited "$projectdir/output/tables/km_survival_tables_rounded.csv", clear varnames(1) stringcols(1 2 3 4 6 7)

quietly count

if r(N) > 0 {

    *****Store combined table
    tempfile survival_tables
    save "`survival_tables'", replace

    *****Identify exposures
    quietly levelsof exposure_name, local(exposures)

    foreach exposure of local exposures {

        use "`survival_tables'", clear
        quietly levelsof outcome if exposure_name == "`exposure'", local(outcomes)

        foreach outcome of local outcomes {

            *****Load current exposure/outcome combination
            use "`survival_tables'", clear
            keep if exposure_name == "`exposure'" & outcome == "`outcome'"

            *****Check one row per category and interval
            isid exposure_group end_week
            assert !missing(exposure_group)

            *****Set legend title
            local legtitle = exposure[1]
            if `"`legtitle'"' == "" local legtitle "`exposure'"

            if inlist("`exposure'", "urate_12m_ult", "urate_12m_ult_cat", "urate_12m_ult_recode", "urate_300_12m_ult", "urate_targets_12m_ult") {
                local legtitle "Urate target attained"
            }

            *****Determine whole-year upper axis limit
            quietly summarize end_years, meanonly
            local km_tmax = ceil(r(max))

            *****Calculate interval survival from rounded counts
            sort exposure_group end_week

            assert !missing(n_at_risk, n_events, n_censored)
            assert n_at_risk > 0
            assert n_events >= 0 & n_censored >= 0
            assert n_at_risk - n_censored/2 > 0
            assert n_events <= n_at_risk - n_censored/2

            gen double interval_survival = 1 - n_events / (n_at_risk - n_censored/2)
            assert inrange(interval_survival, 0, 1)

            *****Calculate cumulative survival
            by exposure_group: gen double survival = exp(sum(ln(interval_survival)))

            *****Keep survival at zero after any zero-survival interval
            by exposure_group: gen long zero_intervals = sum(interval_survival == 0)
            replace survival = 0 if zero_intervals > 0

            *****Log-minus-log transformation
            gen double loglog_survival = -ln(-ln(survival)) if survival > 0 & survival < 1

            *****Build curves and legend from available categories
            quietly levelsof exposure_group, local(groups)

            local colours "emerald orange red blue dkgreen cranberry navy maroon teal sienna purple"
            local legorder
            local km_curves
            local loglog_curves
            local i = 0

            foreach g of local groups {

                local ++i

                *****Retrieve this group's category label
                gen byte current_group = exposure_group == `g'
                gsort -current_group end_week
                local category = exposure_category[1]
                drop current_group

                if `"`category'"' == "" local category "`g'"

                *****Assign colours in ascending category-code order
                local colour : word `i' of `colours'
                if "`colour'" == "" local colour "black"

                local legorder `"`legorder' `i' "`category'""'
                local km_curves `"`km_curves' (line survival end_years if exposure_group == `g', sort connect(stairstep) lcolor(`colour') lpattern(solid))"'
                local loglog_curves `"`loglog_curves' (line loglog_survival end_years if exposure_group == `g', sort lcolor(`colour') lpattern(solid))"'
            }

            *****Reconstructed survival plot
            capture noisily twoway `km_curves', yscale(range(0 1)) ylabel(0(0.2)1, format(%3.1f) nogrid labsize(small)) ytitle("Survival probability", size(medsmall)) xtitle("Years from landmark", size(medsmall) margin(medsmall)) xscale(range(0 `km_tmax')) xlabel(0(1)`km_tmax', nogrid labsize(small)) title("") legend(order(`legorder') title(`"`legtitle'"', size(small) margin(b=1))) xsize(16) ysize(9) name(km_reconstructed, replace)

            if _rc == 0 {
                capture noisily graph export "$projectdir/output/figures/km_recon_`exposure'_`outcome'.$img", name(km_reconstructed) replace
                if _rc == 0 local ++n_km_graphs
            }

            *****Reconstructed log-minus-log plot
            quietly count if !missing(loglog_survival)

            if r(N) > 0 {

                capture noisily twoway `loglog_curves', ytitle("-log{-log(Survival probability)}", size(medsmall)) ylabel(, nogrid labsize(small)) xtitle("Years from landmark", size(medsmall) margin(medsmall)) xscale(range(0 `km_tmax')) xlabel(0(1)`km_tmax', nogrid labsize(small)) title("") legend(order(`legorder') title(`"`legtitle'"', size(small) margin(b=1))) xsize(16) ysize(9) name(loglog_reconstructed, replace)

                if _rc == 0 {
                    capture noisily graph export "$projectdir/output/figures/loglog_recon_`exposure'_`outcome'.$img", name(loglog_reconstructed) replace
                    if _rc == 0 local ++n_loglog_graphs
                }
            }
            else {
                di as text "No finite log-minus-log values for `exposure' / `outcome'; log-log plot skipped."
            }
        }
    }
}
else {
    di as text "No weekly survival data available; reconstructed graphs skipped."
}

*Create dummy KM graph only if no KM graphs were exported
if `n_km_graphs' == 0 {
    preserve
    clear
    set obs 1
    gen x = 1
    gen y = 1

    twoway scatter y x, msymbol(none) xlabel(none) ylabel(none) xtitle("") ytitle("") title("No Kaplan-Meier estimates available") legend(off)

    graph export "$projectdir/output/figures/km_recon_no_outputs.$img", replace
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

    graph export "$projectdir/output/figures/loglog_recon_no_outputs.$img", replace
    restore
}

capture log close