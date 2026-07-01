********************************************************************************
* Replication Code for Main Results
* "Do Power Outages Impact Mental Health? Empirical Evidence from Maryland"
*
* Revised version (June 2026).
* This file consolidates the Stata code used to produce the main tables and
* figures in the manuscript and the supplemental information, organized by
* manuscript order. 
*
* Required packages: ppmlhdfe, reghdfe, estout, coefplot, winsor2, outreg2, ftools
*
* Data availability: the underlying patient-level HCUP Maryland State Inpatient
* Database (SID) is licensed and cannot be redistributed. The raw-data
* processing pipeline (SID cleaning, outage/weather/ACS merges, collapsing to
* the ZIP-quarter panel) is therefore not part of this replication file; the
* Appendix at the end documents how the secondary-diagnosis outcomes are
* constructed from the SID so that the outcome definitions are fully
* transparent. Details of all variable constructions are in the STAR Methods
* and Methods S1-S6 of the Supplemental Information.
*
* Notes:
*   - Figures 5 and 6 (infrastructure density maps) were produced in ArcGIS Pro.
********************************************************************************

clear all
set more off

* >>> EDIT THIS PATH before running: set it to the folder containing the two
* analysis datasets (collapsed_data_for_regression.dta and
* collapsed_data_with_mh_secondary.dta). These datasets are derived from the
* licensed HCUP Maryland SID and are NOT included in this repository; 
cd "YOUR_PROJECT_FOLDER"

* Optional long-running / restricted-data sections (0 = skip, 1 = run)
local run_balanced_zip   0   // PART 3: 500-replication balanced-ZIP check (slow)
local build_secondary    0   // APPENDIX: rebuild secondary outcomes from raw SID

* Install required packages if missing
foreach pkg in ppmlhdfe reghdfe estout coefplot winsor2 outreg2 ftools {
    cap which `pkg'
    if _rc ssc install `pkg'
}

* All tables and figures are written to this subfolder
cap mkdir "replication_outputs"


********************************************************************************
*                                                                              *
*                    PART 1: MAIN ANALYSES                                     *
*                    Data: collapsed_data_for_regression.dta                   *
*                                                                              *
********************************************************************************

use "collapsed_data_for_regression.dta", clear

* Winsorize psychiatrist supply at the 0-99th percentile by year
cap drop psych_w
gen psych_w = psych_per_100k
winsor2 psych_w, replace cuts(0 99) by(AYEAR)

* Control sets used repeatedly below
* Model 3 (main specification) controls:
local ctrl_main "avg_cmr_clean avg_cmr_mortality_clean poverty_rate pct_male median_income_1k black_pct asian_pct pct_medicaid pct_uninsured psych_w"
* Model 5 (disease-specific comorbidity) controls:
local ctrl_com  "pct_psychoses pct_alcohol pct_diabetes_cx pct_hypertension_cx pct_copd avg_cmr_mortality_clean poverty_rate pct_male median_income_1k pct_medicaid pct_uninsured psych_w"
* Common estimation options:
local opts      "absorb(zip_num i.AYEAR#i.DQTR) exposure(total_population) vce(cluster zip_num) eform"


*=============================================================================
* TABLE S5: Summary Statistics - Full Variable List (ZIP-Quarter Level)
*=============================================================================

label var mh_primary_cases "Primary mental health hospitalizations"
label var new_1hr_episodes_p75 "Outage episodes >= 1hr (P75)"
label var ln_outage_1hr_p75 "ln(Outage episodes + 1)"
label var extreme_heat_days_p90 "Extreme heat days (>= P90)"
label var avg_cmr_clean "Comorbidity index (readmission)"
label var avg_cmr_mortality_clean "Comorbidity index (mortality)"
label var median_income_1k "Median income ($1,000s)"
label var poverty_rate "Poverty rate (%)"
label var pct_male "Male patients (%)"
label var pct_medicaid "Medicaid (%)"
label var pct_uninsured "Uninsured (%)"
label var psych_w "Psychiatrists per 100,000 (winsorized)"
label var total_population "ZIP code population"

local s5_vars "mh_primary_cases new_1hr_episodes_p75 ln_outage_1hr_p75 extreme_heat_days_p90 avg_cmr_clean avg_cmr_mortality_clean pct_psychoses pct_alcohol pct_drug_abuse pct_diabetes_cx pct_hypertension_cx pct_copd pct_renal_disease median_income_1k poverty_rate total_population mean_age pct_male black_pct asian_pct pct_medicaid pct_uninsured psych_w"

estpost summarize `s5_vars', detail
esttab using "replication_outputs/TableS5_summary_statistics.rtf", replace ///
    cells("count(fmt(0)) mean(fmt(2)) sd(fmt(2)) p25(fmt(2)) p75(fmt(2))") ///
    noobs nonumber nomtitle ///
    collabels("N" "Mean" "SD" "P25" "P75") ///
    title("Table S5: Summary Statistics - Full Variable List") ///
    addnotes("Unit of observation: ZIP code x year x quarter, Maryland 2018-2023." ///
             "Descriptive statistics use all available data; the estimation sample" ///
             "for the main models (Table S9) is N = 8,152 after merging and" ///
             "dropping fixed-effect singletons.")


*=============================================================================
* FIGURE 1 + TABLE S9: Baseline PPML-HDFE Estimates (5 Model Specifications)
*=============================================================================

// Model 1: Basic socioeconomic controls + comorbidity indices
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 ///
    poverty_rate median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    , `opts'
estimates store model1

// Model 2: + Extreme weather (heat, precipitation, cold at P90)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 extreme_precip_days_p90 extreme_cold_days_p90 ///
    poverty_rate median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    , `opts'
estimates store model2

// Model 3: Full demographic model (main specification)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    `ctrl_main' ///
    , `opts'
estimates store model3

// Model 4: Alternative heat definition (P75 threshold)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 unusual_heat_days_p75 ///
    `ctrl_main' ///
    , `opts'
estimates store model4

// Model 5: Disease-specific comorbidity indicators
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p95 ///
    `ctrl_com' ///
    , `opts'
estimates store model5

* Export Table S9
esttab model1 model2 model3 model4 model5 using "replication_outputs/TableS9_main_results.rtf", replace ///
    eform b(3) se(3) ///
    title("Table S9: PPML Results - Mental Health Hospitalizations and Power Outages") ///
    mtitle("Basic" "+ Weather" "Full Demo" "Alt Heat" "Comorbid") ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) compress ///
    addnotes("Fixed effects: ZIP code and Year x Quarter. Exposure: total population." ///
             "Standard errors clustered at ZIP code level.")

*--- Figure 1 Panel (a): Power Outage Exposure across 5 models ---

* Extract IRR labels with significance stars for each model
foreach m in 1 2 3 4 5 {
    estimates restore model`m'
    local b = exp(_b[ln_outage_1hr_p75])
    local se = _se[ln_outage_1hr_p75]
    local z = _b[ln_outage_1hr_p75]/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_model`m' : di %5.3f `b'
    local lab_model`m' "`lab_model`m''`stars'"
}

coefplot (model1, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Model 1: Basic") mcolor(blue) ciopts(lcolor(blue)) mlabcolor(blue) mlabel("`lab_model1'")) ///
         (model2, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Model 2: + Weather") mcolor(cranberry) ciopts(lcolor(cranberry)) mlabcolor(cranberry) mlabel("`lab_model2'")) ///
         (model3, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Model 3: Main") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_model3'")) ///
         (model4, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Model 4: Alt. Weather") mcolor(gold) ciopts(lcolor(gold)) mlabcolor(gold) mlabel("`lab_model4'")) ///
         (model5, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Model 5: All Controls") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_model5'")), ///
    xline(1, lcolor(gs10) lpattern(dash)) eform ///
    xlabel(1.000(0.010)1.040) xtitle("Incidence Rate Ratio (IRR)") ///
    title("a. Power Outage Exposure", size(medium)) ///
    mlabposition(12) mlabsize(medium) ///
    legend(off) ///
    graphregion(color(white)) bgcolor(white)
graph export "replication_outputs/panel_a_outage.png", replace width(1200)

*--- Figure 1 Panel (b): Weather Controls ---

estimates restore model3
local b = exp(_b[extreme_heat_days_p90])
local se = _se[extreme_heat_days_p90]
local z = _b[extreme_heat_days_p90]/`se'
local p = 2*(1-normal(abs(`z')))
if `p' < 0.01      local stars "***"
else if `p' < 0.05 local stars "**"
else if `p' < 0.10 local stars "*"
else                local stars ""
local lab_heat_m3 : di %5.3f `b'
local lab_heat_m3 "`lab_heat_m3'`stars'"

estimates restore model2
foreach var in extreme_precip_days_p90 extreme_cold_days_p90 {
    local b = exp(_b[`var'])
    local se = _se[`var']
    local z = _b[`var']/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_`var' : di %5.3f `b'
    local lab_`var' "`lab_`var''`stars'"
}

estimates restore model4
local b = exp(_b[unusual_heat_days_p75])
local se = _se[unusual_heat_days_p75]
local z = _b[unusual_heat_days_p75]/`se'
local p = 2*(1-normal(abs(`z')))
if `p' < 0.01      local stars "***"
else if `p' < 0.05 local stars "**"
else if `p' < 0.10 local stars "*"
else                local stars ""
local lab_uheat_m4 : di %5.3f `b'
local lab_uheat_m4 "`lab_uheat_m4'`stars'"

estimates restore model5
local b = exp(_b[extreme_heat_days_p95])
local se = _se[extreme_heat_days_p95]
local z = _b[extreme_heat_days_p95]/`se'
local p = 2*(1-normal(abs(`z')))
if `p' < 0.01      local stars "***"
else if `p' < 0.05 local stars "**"
else if `p' < 0.10 local stars "*"
else                local stars ""
local lab_heat_m5 : di %5.3f `b'
local lab_heat_m5 "`lab_heat_m5'`stars'"

coefplot (model3, keep(extreme_heat_days_p90) rename(extreme_heat_days_p90 = "Extreme Heat (P90) [M3]") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_heat_m3'")) ///
         (model2, keep(extreme_precip_days_p90) rename(extreme_precip_days_p90 = "Extreme Precip (P90) [M2]") mcolor(cranberry) ciopts(lcolor(cranberry)) mlabcolor(cranberry) mlabel("`lab_extreme_precip_days_p90'")) ///
         (model2, keep(extreme_cold_days_p90) rename(extreme_cold_days_p90 = "Extreme Cold (P90) [M2]") mcolor(cranberry) ciopts(lcolor(cranberry)) mlabcolor(cranberry) mlabel("`lab_extreme_cold_days_p90'")) ///
         (model4, keep(unusual_heat_days_p75) rename(unusual_heat_days_p75 = "Unusual Heat (P75) [M4]") mcolor(gold) ciopts(lcolor(gold)) mlabcolor(gold) mlabel("`lab_uheat_m4'")) ///
         (model5, keep(extreme_heat_days_p95) rename(extreme_heat_days_p95 = "Extreme Heat (P95) [M5]") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_heat_m5'")), ///
    xline(1, lcolor(gs10) lpattern(dash)) eform ///
    xlabel(0.996(0.002)1.008) xtitle("Incidence Rate Ratio (IRR)") ///
    title("b. Weather Controls", size(medium)) ///
    mlabposition(12) mlabsize(medium) ///
    legend(off) ///
    graphregion(color(white)) bgcolor(white)
graph export "replication_outputs/panel_b_weather.png", replace width(1200)

*--- Figure 1 Panel (c): Socioeconomic, Demographic & Healthcare Access (Model 3) ---

estimates restore model3
foreach var in poverty_rate median_income_1k pct_male pct_medicaid pct_uninsured psych_w {
    local b = exp(_b[`var'])
    local se = _se[`var']
    local z = _b[`var']/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_`var' : di %5.3f `b'
    local lab_`var' "`lab_`var''`stars'"
}

coefplot (model3, keep(poverty_rate) rename(poverty_rate = "Poverty Rate") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_poverty_rate'")) ///
         (model3, keep(median_income_1k) rename(median_income_1k = "Median Income (1k)") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_median_income_1k'")) ///
         (model3, keep(pct_male) rename(pct_male = "Pct Male") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_pct_male'")) ///
         (model3, keep(pct_medicaid) rename(pct_medicaid = "Pct Medicaid") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_pct_medicaid'")) ///
         (model3, keep(pct_uninsured) rename(pct_uninsured = "Pct Uninsured") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_pct_uninsured'")) ///
         (model3, keep(psych_w) rename(psych_w = "Psychiatrist Supply") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_psych_w'")), ///
    xline(1, lcolor(gs10) lpattern(dash)) eform ///
    xlabel(1.00(0.01)1.04) xtitle("Incidence Rate Ratio (IRR)") ///
    title("c. Socioeconomic, Demographic & Healthcare Access (Model 3)", size(medium)) ///
    mlabposition(12) mlabsize(medium) ///
    legend(off) ///
    graphregion(color(white)) bgcolor(white)
graph export "replication_outputs/panel_c_socioeconomic.png", replace width(1200)

*--- Figure 1 Panel (d): Health Comorbidity Controls ---

estimates restore model3
foreach var in avg_cmr_clean avg_cmr_mortality_clean {
    local b = exp(_b[`var'])
    local se = _se[`var']
    local z = _b[`var']/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_`var' : di %5.3f `b'
    local lab_`var' "`lab_`var''`stars'"
}

estimates restore model5
foreach var in pct_psychoses pct_alcohol pct_diabetes_cx pct_hypertension_cx pct_copd {
    local b = exp(_b[`var'])
    local se = _se[`var']
    local z = _b[`var']/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_`var' : di %5.3f `b'
    local lab_`var' "`lab_`var''`stars'"
}

coefplot (model3, keep(avg_cmr_clean) rename(avg_cmr_clean = "CMR Readmission [M3]") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_avg_cmr_clean'")) ///
         (model3, keep(avg_cmr_mortality_clean) rename(avg_cmr_mortality_clean = "CMR Mortality [M3]") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_avg_cmr_mortality_clean'")) ///
         (model5, keep(pct_psychoses) rename(pct_psychoses = "Psychoses [M5]") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_pct_psychoses'")) ///
         (model5, keep(pct_alcohol) rename(pct_alcohol = "Alcohol Use [M5]") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_pct_alcohol'")) ///
         (model5, keep(pct_diabetes_cx) rename(pct_diabetes_cx = "Diabetes [M5]") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_pct_diabetes_cx'")) ///
         (model5, keep(pct_hypertension_cx) rename(pct_hypertension_cx = "Hypertension [M5]") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_pct_hypertension_cx'")) ///
         (model5, keep(pct_copd) rename(pct_copd = "COPD [M5]") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_pct_copd'")), ///
    xline(1, lcolor(gs10) lpattern(dash)) eform ///
    xlabel(0.925(0.025)1.075) xtitle("Incidence Rate Ratio (IRR)") ///
    title("d. Health Comorbidity Controls", size(medium)) ///
    mlabposition(12) mlabsize(medium) ///
    legend(off) ///
    graphregion(color(white)) bgcolor(white)
graph export "replication_outputs/panel_d_comorbidity.png", replace width(1200)


*=============================================================================
* FIGURE 2 + TABLE S10: Robustness - Alternative Outage Threshold Definitions
*=============================================================================

// P50 threshold (1hr, 50th percentile)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p50 extreme_heat_days_p90 `ctrl_main', `opts'
estimates store outage_p50

// P90 threshold (1hr, 90th percentile)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p90 extreme_heat_days_p90 `ctrl_main', `opts'
estimates store outage_p90

// 2hr P50 threshold
ppmlhdfe mh_primary_cases ///
    ln_outage_2hr_p50 extreme_heat_days_p90 `ctrl_main', `opts'
estimates store outage_2hr_p50

* Export Table S10
esttab outage_p50 outage_p90 outage_2hr_p50 using "replication_outputs/TableS10_robustness_outage.rtf", replace ///
    eform b(3) se(3) ///
    title("Table S10: Robustness - Different Outage Thresholds") ///
    mtitle("1hr P50" "1hr P90" "2hr P50") ///
    stats(N N_clust, labels("N" "Clusters") fmt(0 0)) ///
    star(* 0.10 ** 0.05 *** 0.01) compress

* Figure 2: Robustness coefplot
foreach est in outage_p50 outage_p90 outage_2hr_p50 {
    estimates restore `est'
    if "`est'" == "outage_p50"     local vname "ln_outage_1hr_p50"
    if "`est'" == "outage_p90"     local vname "ln_outage_1hr_p90"
    if "`est'" == "outage_2hr_p50" local vname "ln_outage_2hr_p50"

    local b = exp(_b[`vname'])
    local se = _se[`vname']
    local z = _b[`vname']/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_`est' : di %5.3f `b'
    local lab_`est' "`lab_`est''`stars'"
}

coefplot (outage_p50, keep(ln_outage_1hr_p50) rename(ln_outage_1hr_p50 = ">= 1 hour, 50th pctl") mcolor(blue) ciopts(lcolor(blue)) mlabcolor(blue) mlabel("`lab_outage_p50'")) ///
         (outage_p90, keep(ln_outage_1hr_p90) rename(ln_outage_1hr_p90 = ">= 1 hour, 90th pctl") mcolor(blue) ciopts(lcolor(blue)) mlabcolor(blue) mlabel("`lab_outage_p90'")) ///
         (outage_2hr_p50, keep(ln_outage_2hr_p50) rename(ln_outage_2hr_p50 = ">= 2 hours, 50th pctl") mcolor(blue) ciopts(lcolor(blue)) mlabcolor(blue) mlabel("`lab_outage_2hr_p50'")), ///
    xline(1, lcolor(gs10) lpattern(dash)) eform ///
    xlabel(1.00(0.01)1.04) xtitle("Incidence Rate Ratio (IRR)") ///
    title("Robustness: Power Outage Effect Across Alternative Threshold Definitions", size(medsmall)) ///
    mlabposition(12) mlabsize(medium) ///
    legend(off) ///
    ysize(6) xsize(10) ///
    graphregion(color(white)) bgcolor(white)
graph export "replication_outputs/figure2_robustness_outage.png", replace width(1200)


*=============================================================================
* FIGURE 3 + TABLE S11: Sensitivity - Duration Thresholds (1hr-8hr, P50)
*=============================================================================

* Table S11: descriptive distribution of the log-transformed outage measures
estpost summarize ln_outage_1hr_p50 ln_outage_2hr_p50 ln_outage_3hr_p50 ///
    ln_outage_4hr_p50 ln_outage_5hr_p50 ln_outage_6hr_p50 ///
    ln_outage_7hr_p50 ln_outage_8hr_p50, detail
esttab using "replication_outputs/TableS11_duration_distribution.rtf", replace ///
    cells("mean(fmt(2)) sd(fmt(2)) min(fmt(2)) p50(fmt(2)) max(fmt(2))") ///
    noobs nonumber nomtitle ///
    collabels("Mean" "SD" "Min" "p50" "Max") ///
    title("Table S11: Outage Episodes by Duration Threshold (Log-transformed)")

* Run PPML for each duration threshold (plotted in Figure 3)
foreach dur in 1hr 2hr 3hr 4hr 5hr 6hr 7hr 8hr {
    ppmlhdfe mh_primary_cases ln_outage_`dur'_p50 extreme_heat_days_p90 ///
        `ctrl_main', `opts'
    estimates store dur_`dur'
}

* Supporting table for Figure 3
esttab dur_1hr dur_2hr dur_3hr dur_4hr dur_5hr dur_6hr dur_7hr dur_8hr ///
    using "replication_outputs/Figure3_duration_regressions.rtf", replace ///
    eform b(3) se(3) ///
    title("Sensitivity Analysis: Outage Duration Thresholds (P50)") ///
    mtitle(">=1hr" ">=2hr" ">=3hr" ">=4hr" ">=5hr" ">=6hr" ">=7hr" ">=8hr") ///
    keep(ln_outage_*) ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    addnotes("Robust standard errors in parentheses, clustered at ZIP code level." ///
             "All models include ZIP code and Year x Quarter fixed effects." ///
             "Exposure: total population.")

* Figure 3: Duration sensitivity coefplot
coefplot dur_1hr dur_2hr dur_3hr dur_4hr dur_5hr dur_6hr dur_7hr dur_8hr, ///
    keep(ln_outage_*) ///
    eform ///
    vertical ///
    yline(1, lcolor(red) lpattern(dash)) ///
    ylabel(0.98(0.01)1.05, format(%4.2f)) ///
    xtitle("Minimum Outage Duration Threshold") ///
    ytitle("Incidence Rate Ratio (IRR)") ///
    title("Effect of Power Outages on Mental Health Hospitalizations" ///
          "by Outage Duration Threshold") ///
    coeflabels(ln_outage_1hr_p50 = ">=1hr***" ///
               ln_outage_2hr_p50 = ">=2hr**" ///
               ln_outage_3hr_p50 = ">=3hr**" ///
               ln_outage_4hr_p50 = ">=4hr" ///
               ln_outage_5hr_p50 = ">=5hr" ///
               ln_outage_6hr_p50 = ">=6hr" ///
               ln_outage_7hr_p50 = ">=7hr" ///
               ln_outage_8hr_p50 = ">=8hr") ///
    msymbol(D) mcolor(navy) ///
    ciopts(lcolor(navy)) ///
    graphregion(color(white)) ///
    note("Notes: Each point shows the IRR from a separate PPML regression." ///
         "Outage measure: ln(1 + episodes), where episodes = outage events >= X hours" ///
         "and above the 50th percentile of severity in each ZIP-quarter." ///
         "95% CI shown. *** p<0.01, ** p<0.05, * p<0.1. All models include ZIP and Year x Quarter FE.", size(vsmall))
graph export "replication_outputs/figure3_sensitivity_duration.png", replace width(1200)


*=============================================================================
* TABLE 1 (Manuscript): Effect Modification by Extreme Heat
*                       (Outage x Heat Interaction)
*=============================================================================

* Create mean-centered variables
capture drop ln_outage_c heat_c mean_outage mean_heat
egen double mean_outage = mean(ln_outage_1hr_p75)
egen double mean_heat = mean(extreme_heat_days_p90)
gen double ln_outage_c = ln_outage_1hr_p75 - mean_outage
gen double heat_c = extreme_heat_days_p90 - mean_heat
label var ln_outage_c "Log outage episodes (centered)"
label var heat_c "Extreme heat days (centered)"

* Run interaction model
ppmlhdfe mh_primary_cases c.ln_outage_c##c.heat_c `ctrl_main', `opts'
eststo heat_interact
test c.ln_outage_c#c.heat_c

* Export Table 1
esttab heat_interact using "replication_outputs/Table1_heat_interaction.rtf", replace ///
    eform b(4) se(3) ///
    keep(ln_outage_c heat_c c.ln_outage_c#c.heat_c) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N N_clust, labels("N" "Clusters") fmt(0 0))


*=============================================================================
* TABLE 2 (Manuscript): IV Poisson Estimates
*                       3 Instrument Sets: Transformer, Substation+Transmission, All 3
*=============================================================================

* Controls used in the IV models (exogenous regressors)
local ivctrl "extreme_heat_days_p90 pct_psychoses pct_alcohol poverty_rate pct_male median_income_1k avg_cmr_clean avg_cmr_mortality_clean black_pct asian_pct pct_medicaid pct_uninsured psych_w"

*--- Column (1): Transformer density as sole IV ---

* First stage
reghdfe ln_outage_1hr_p75 dens_transformer `ivctrl' ///
    , absorb(i.AYEAR) cluster(zip_num)
test dens_transformer
estadd scalar FirstStage_F = r(F)
estimates store fs_trans

* Second stage
ivpoisson gmm mh_primary_cases `ivctrl' i.AYEAR ///
    (ln_outage_1hr_p75 = dens_transformer) ///
    , exposure(total_population) vce(cluster zip_num)
estimates store ss_trans

*--- Column (2): Substation + Transmission line density as IVs ---

* First stage
reghdfe ln_outage_1hr_p75 dens_substation dens_transmission_line `ivctrl' ///
    , absorb(i.AYEAR) cluster(zip_num)
test dens_substation dens_transmission_line
estadd scalar FirstStage_F = r(F)
estimates store fs_sub_trans

* Second stage
ivpoisson gmm mh_primary_cases `ivctrl' i.AYEAR ///
    (ln_outage_1hr_p75 = dens_substation dens_transmission_line) ///
    , exposure(total_population) vce(cluster zip_num)
estat overid
estadd scalar hansen_p = chi2tail(e(J_df), e(J))
estimates store ss_sub_trans

*--- Column (3): All 3 IVs ---

* First stage
reghdfe ln_outage_1hr_p75 dens_transformer dens_transmission_line dens_substation `ivctrl' ///
    , absorb(i.AYEAR) cluster(zip_num)
test dens_transformer dens_transmission_line dens_substation
estadd scalar FirstStage_F = r(F)
estimates store fs_all3

* Second stage
ivpoisson gmm mh_primary_cases `ivctrl' i.AYEAR ///
    (ln_outage_1hr_p75 = dens_transformer dens_transmission_line dens_substation) ///
    , exposure(total_population) vce(cluster zip_num)
estat overid
estadd scalar hansen_p = chi2tail(e(J_df), e(J))
estimates store ss_all3

* Export Table 2 (Second Stage)
esttab ss_trans ss_sub_trans ss_all3 using "replication_outputs/Table2_IV_comparison.rtf", replace ///
    b(3) se(3) ///
    title("Table 2: IV Estimates with Alternative Instrument Sets") ///
    mtitle("Transformer Only" "Substation + Trans" "All 3 IVs") ///
    keep(ln_outage_1hr_p75) ///
    coeflabels(ln_outage_1hr_p75 "Log outage exposure") ///
    stats(N hansen_p, labels("Observations" "Hansen J p-value") fmt(0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    addnotes("Column 1: Transformer density as single IV (exactly identified)." ///
             "Column 2: Substation and transmission line density." ///
             "Column 3: All three infrastructure measures." ///
             "Hansen J test not available for exactly identified model.")

* Export First Stage results
esttab fs_trans fs_sub_trans fs_all3 using "replication_outputs/Table2_IV_firststage.rtf", replace ///
    b(3) se(3) ///
    title("Table 2 Panel B: First Stage Results") ///
    mtitle("Transformer Only" "Substation + Trans" "All 3 IVs") ///
    keep(dens_transformer dens_transmission_line dens_substation) ///
    stats(N FirstStage_F, labels("Observations" "F-statistic") fmt(0 2)) ///
    star(* 0.10 ** 0.05 *** 0.01)


*=============================================================================
* FIGURE 4 + TABLE S13: Heterogeneity Analysis
*                       (Urban/Rural, Poverty, Medicaid, Season, Gender)
* Table S13 reports, for each dimension: subgroup IRRs and p-values, the
* z-test for the subgroup difference, and the interaction-model p-value.
*=============================================================================

*--- 1. Urban vs Rural (Large Metro RUCC=1 vs Rest) ---

cap drop large_metro
gen large_metro = (PL_RUCC == 1) if !missing(PL_RUCC)

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if large_metro == 1, `opts'
estimates store large_metro_grp

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if large_metro == 0, `opts'
estimates store non_large_metro_grp

* Interaction test
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 c.ln_outage_1hr_p75#i.large_metro ///
    i.large_metro ///
    extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if !missing(large_metro) ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num)
testparm c.ln_outage_1hr_p75#i.large_metro
local pint_metro = r(p)

* Z-test for subgroup difference
estimates restore large_metro_grp
local b1 = _b[ln_outage_1hr_p75]
local se1 = _se[ln_outage_1hr_p75]
estimates restore non_large_metro_grp
local b2 = _b[ln_outage_1hr_p75]
local se2 = _se[ln_outage_1hr_p75]
local z_diff = (`b1' - `b2') / sqrt(`se1'^2 + `se2'^2)
local p_metro = 2*(1-normal(abs(`z_diff')))
di "TABLE S13 [Urban/Rural]: z-test p=" %6.4f `p_metro' "  interaction p=" %6.4f `pint_metro'

*--- 2. Poverty Level (Census Bureau 20% threshold) ---

cap drop high_poverty_census
gen high_poverty_census = (poverty_rate >= 20) if !missing(poverty_rate)

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_poverty_census == 0, `opts'
estimates store low_pov_census

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_poverty_census == 1, `opts'
estimates store high_pov_census

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 c.ln_outage_1hr_p75#i.high_poverty_census ///
    i.high_poverty_census ///
    extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if !missing(high_poverty_census) ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num)
testparm c.ln_outage_1hr_p75#i.high_poverty_census
local pint_pov = r(p)

estimates restore low_pov_census
local b1 = _b[ln_outage_1hr_p75]
local se1 = _se[ln_outage_1hr_p75]
estimates restore high_pov_census
local b2 = _b[ln_outage_1hr_p75]
local se2 = _se[ln_outage_1hr_p75]
local z_diff = (`b1' - `b2') / sqrt(`se1'^2 + `se2'^2)
local p_pov = 2*(1-normal(abs(`z_diff')))
di "TABLE S13 [Poverty]: z-test p=" %6.4f `p_pov' "  interaction p=" %6.4f `pint_pov'

*--- 3. Medicaid Coverage (25% threshold) ---

cap drop high_medicaid_25
gen high_medicaid_25 = (pct_medicaid >= 25) if !missing(pct_medicaid)

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    median_income_1k psych_w ///
    if high_medicaid_25 == 0, `opts'
estimates store low_medicaid

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    median_income_1k psych_w ///
    if high_medicaid_25 == 1, `opts'
estimates store high_medicaid_grp

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 c.ln_outage_1hr_p75#i.high_medicaid_25 ///
    i.high_medicaid_25 ///
    extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    median_income_1k psych_w ///
    if !missing(high_medicaid_25) ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num)
testparm c.ln_outage_1hr_p75#i.high_medicaid_25
local pint_mcd = r(p)

estimates restore low_medicaid
local b1 = _b[ln_outage_1hr_p75]
local se1 = _se[ln_outage_1hr_p75]
estimates restore high_medicaid_grp
local b2 = _b[ln_outage_1hr_p75]
local se2 = _se[ln_outage_1hr_p75]
local z_diff = (`b1' - `b2') / sqrt(`se1'^2 + `se2'^2)
local p_mcd = 2*(1-normal(abs(`z_diff')))
di "TABLE S13 [Medicaid]: z-test p=" %6.4f `p_mcd' "  interaction p=" %6.4f `pint_mcd'

*--- 4. Season (Warm Q2-Q3 vs Cold Q1&Q4) ---

cap drop warm_season
gen warm_season = (DQTR == 2 | DQTR == 3) if !missing(DQTR)

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if warm_season == 0, `opts'
estimates store cold_season

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if warm_season == 1, `opts'
estimates store warm_season_grp

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 c.ln_outage_1hr_p75#i.warm_season ///
    i.warm_season ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if !missing(warm_season) ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num)
testparm c.ln_outage_1hr_p75#i.warm_season
local pint_season = r(p)

estimates restore cold_season
local b1 = _b[ln_outage_1hr_p75]
local se1 = _se[ln_outage_1hr_p75]
estimates restore warm_season_grp
local b2 = _b[ln_outage_1hr_p75]
local se2 = _se[ln_outage_1hr_p75]
local z_diff = (`b1' - `b2') / sqrt(`se1'^2 + `se2'^2)
local p_season = 2*(1-normal(abs(`z_diff')))
di "TABLE S13 [Season]: z-test p=" %6.4f `p_season' "  interaction p=" %6.4f `pint_season'

*--- 5. Gender Composition (median split) ---

cap drop high_male
sum pct_male, detail
gen high_male = (pct_male > r(p50)) if !missing(pct_male)

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_male == 0, `opts'
estimates store low_male

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_male == 1, `opts'
estimates store high_male

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 c.ln_outage_1hr_p75#i.high_male ///
    i.high_male ///
    extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if !missing(high_male) ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num)
testparm c.ln_outage_1hr_p75#i.high_male
local pint_male = r(p)

estimates restore low_male
local b1 = _b[ln_outage_1hr_p75]
local se1 = _se[ln_outage_1hr_p75]
estimates restore high_male
local b2 = _b[ln_outage_1hr_p75]
local se2 = _se[ln_outage_1hr_p75]
local z_diff = (`b1' - `b2') / sqrt(`se1'^2 + `se2'^2)
local p_male = 2*(1-normal(abs(`z_diff')))
di "TABLE S13 [Gender]: z-test p=" %6.4f `p_male' "  interaction p=" %6.4f `pint_male'

*--- Table S13 summary (printed to log/results window) ---
di _newline "==================== TABLE S13 SUMMARY ===================="
di "Urban/Rural: z-test p=" %6.4f `p_metro'  "  interaction p=" %6.4f `pint_metro'
di "Poverty:     z-test p=" %6.4f `p_pov'    "  interaction p=" %6.4f `pint_pov'
di "Medicaid:    z-test p=" %6.4f `p_mcd'    "  interaction p=" %6.4f `pint_mcd'
di "Season:      z-test p=" %6.4f `p_season' "  interaction p=" %6.4f `pint_season'
di "Gender:      z-test p=" %6.4f `p_male'   "  interaction p=" %6.4f `pint_male'
di "==========================================================="

*--- Export full subgroup regression tables (supporting Figure 4 / Table S13) ---

outreg2 [large_metro_grp non_large_metro_grp low_pov_census high_pov_census] ///
    using "replication_outputs/TableS13_hetero_partA.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Heterogeneity Analysis - Part A") ///
    ctitle("Large Metro", "Non-Large Metro", "Low Poverty", "High Poverty") ///
    addtext(ZIP FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at ZIP code level.")

outreg2 [low_medicaid high_medicaid_grp cold_season warm_season_grp] ///
    using "replication_outputs/TableS13_hetero_partB.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Heterogeneity Analysis - Part B") ///
    ctitle("Low Medicaid", "High Medicaid", "Cold Season", "Warm Season") ///
    addtext(ZIP FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at ZIP code level.")

outreg2 [low_male high_male] ///
    using "replication_outputs/TableS13_hetero_partC.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Heterogeneity Analysis - Part C") ///
    ctitle("Lower Male %", "Higher Male %") ///
    addtext(ZIP FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at ZIP code level.")

*--- Figure 4: Heterogeneity coefplot ---

foreach est in large_metro_grp non_large_metro_grp low_pov_census high_pov_census ///
               low_medicaid high_medicaid_grp cold_season warm_season_grp low_male high_male {
    estimates restore `est'
    local b = exp(_b[ln_outage_1hr_p75])
    local se = _se[ln_outage_1hr_p75]
    local z = _b[ln_outage_1hr_p75]/`se'
    local p = 2*(1-normal(abs(`z')))
    if `p' < 0.01      local stars "***"
    else if `p' < 0.05 local stars "**"
    else if `p' < 0.10 local stars "*"
    else                local stars ""
    local lab_`est' : di %5.3f `b'
    local lab_`est' "`lab_`est''`stars'"
}

coefplot (large_metro_grp, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Large Metro") mcolor(blue) ciopts(lcolor(blue)) mlabcolor(blue) mlabel("`lab_large_metro_grp'")) ///
         (non_large_metro_grp, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Non-Large Metro") mcolor(blue) ciopts(lcolor(blue)) mlabcolor(blue) mlabel("`lab_non_large_metro_grp'")) ///
         (low_pov_census, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Low Poverty") mcolor(cranberry) ciopts(lcolor(cranberry)) mlabcolor(cranberry) mlabel("`lab_low_pov_census'")) ///
         (high_pov_census, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "High Poverty") mcolor(cranberry) ciopts(lcolor(cranberry)) mlabcolor(cranberry) mlabel("`lab_high_pov_census'")) ///
         (low_medicaid, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Low Medicaid") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_low_medicaid'")) ///
         (high_medicaid_grp, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "High Medicaid") mcolor(forest_green) ciopts(lcolor(forest_green)) mlabcolor(forest_green) mlabel("`lab_high_medicaid_grp'")) ///
         (cold_season, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Cold Season") mcolor(gold) ciopts(lcolor(gold)) mlabcolor(gold) mlabel("`lab_cold_season'")) ///
         (warm_season_grp, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Warm Season") mcolor(gold) ciopts(lcolor(gold)) mlabcolor(gold) mlabel("`lab_warm_season_grp'")) ///
         (low_male, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "Low Male %") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_low_male'")) ///
         (high_male, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = "High Male %") mcolor(purple) ciopts(lcolor(purple)) mlabcolor(purple) mlabel("`lab_high_male'")), ///
    xline(1, lcolor(gs10) lpattern(dash)) eform ///
    xlabel(0.96 0.98 1.00 1.02 1.04 1.06, labsize(small)) ///
    xtitle("Incidence Rate Ratio (IRR)") ///
    mlabposition(12) mlabsize(medsmall) ///
    legend(off) ///
    headings("Large Metro" = "{bf:Urban/Rural}" ///
             "Low Poverty" = "{bf:Poverty Level}" ///
             "Low Medicaid" = "{bf:Medicaid Coverage}" ///
             "Cold Season" = "{bf:Seasonality}" ///
             "Low Male %" = "{bf:Gender Composition}") ///
    graphregion(color(white)) bgcolor(white) ///
    xsize(12) ysize(10)
graph export "replication_outputs/figure4_heterogeneity.png", replace width(1600)


*=============================================================================
* TABLE S12: City x Year FE Robustness (aggregation robustness check)
*=============================================================================

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_main' ///
    , absorb(city_id#i.AYEAR i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster city_id) eform
estimates store m_city_year_fe

outreg2 [m_city_year_fe] using "replication_outputs/TableS12_cityFE.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Table S12: PPML Results - City x Year FE Robustness") ///
    ctitle("City x Year FE + Year x Quarter FE") ///
    addtext(City x Year FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at city level.")


*=============================================================================
* TABLE S14: Raw Outage Episode Count as Alternative Exposure Measure
*=============================================================================

ppmlhdfe mh_primary_cases ///
    new_1hr_episodes_p75 extreme_heat_days_p90 `ctrl_main', `opts'
estimates store model_raw_count

esttab model_raw_count using "replication_outputs/TableS14_raw_count.rtf", replace ///
    eform b(3) se(3) ///
    title("Table S14: Robustness - Raw Outage Episode Count as Exposure") ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) compress


*=============================================================================
* TABLE S15 (NEW, iScience revision): Outage effect under each extreme-heat
* threshold (P75/P90/P95) within the same specification
* Columns (1)-(3): Model 3 controls; columns (4)-(6): Model 5 controls.
* Column (1) = manuscript Model 4; (2) = Model 3; (6) = Model 5.
*=============================================================================

eststo clear
eststo m3_p75: ppmlhdfe mh_primary_cases ln_outage_1hr_p75 unusual_heat_days_p75  `ctrl_main', `opts'
eststo m3_p90: ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_main', `opts'
eststo m3_p95: ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p95 `ctrl_main', `opts'
eststo m5_p75: ppmlhdfe mh_primary_cases ln_outage_1hr_p75 unusual_heat_days_p75  `ctrl_com', `opts'
eststo m5_p90: ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_com', `opts'
eststo m5_p95: ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p95 `ctrl_com', `opts'

esttab m3_p75 m3_p90 m3_p95 m5_p75 m5_p90 m5_p95 using "replication_outputs/TableS15_heat_threshold.rtf", replace ///
    eform b(3) se(3) ///
    keep(ln_outage_1hr_p75 unusual_heat_days_p75 extreme_heat_days_p90 extreme_heat_days_p95) ///
    title("Table S15: Robustness - power outage effect under each extreme-heat threshold within the same specification") ///
    mtitle("Main P75" "Main P90" "Main P95" "Comorb P75" "Comorb P90" "Comorb P95") ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) compress ///
    addnotes("Fixed effects: ZIP code and Year x Quarter. Exposure: total population." ///
             "Standard errors clustered at ZIP code level." ///
             "Columns (1)-(3): Model 3 controls. Columns (4)-(6): Model 5 controls." ///
             "Column (1) = manuscript Model 4; column (2) = Model 3; column (6) = Model 5.")


*=============================================================================
* FIGURE S1: IV Robustness - Alternative Outage Definitions (all 3 IVs)
*=============================================================================

foreach dur_var in ln_outage_1hr_p75 ln_outage_2hr_p75 ln_outage_3hr_p75 ///
                   ln_outage_4hr_p75 ln_outage_5hr_p50 ln_outage_6hr_p50 ln_outage_8hr_p50 {

    * First stage
    reghdfe `dur_var' ///
        dens_transformer dens_transmission_line dens_substation `ivctrl' ///
        , absorb(i.AYEAR) cluster(zip_num)
    test dens_transformer dens_transmission_line dens_substation

    * Second stage
    ivpoisson gmm mh_primary_cases `ivctrl' i.AYEAR ///
        (`dur_var' = dens_transformer dens_transmission_line dens_substation) ///
        , exposure(total_population) vce(cluster zip_num)
    estimates store iv_`dur_var'
    estat overid
}

coefplot (iv_ln_outage_1hr_p75, keep(ln_outage_1hr_p75) rename(ln_outage_1hr_p75 = ">=1hr, P75")) ///
         (iv_ln_outage_2hr_p75, keep(ln_outage_2hr_p75) rename(ln_outage_2hr_p75 = ">=2hr, P75")) ///
         (iv_ln_outage_3hr_p75, keep(ln_outage_3hr_p75) rename(ln_outage_3hr_p75 = ">=3hr, P75")) ///
         (iv_ln_outage_4hr_p75, keep(ln_outage_4hr_p75) rename(ln_outage_4hr_p75 = ">=4hr, P75")) ///
         (iv_ln_outage_5hr_p50, keep(ln_outage_5hr_p50) rename(ln_outage_5hr_p50 = ">=5hr, P50")) ///
         (iv_ln_outage_6hr_p50, keep(ln_outage_6hr_p50) rename(ln_outage_6hr_p50 = ">=6hr, P50")) ///
         (iv_ln_outage_8hr_p50, keep(ln_outage_8hr_p50) rename(ln_outage_8hr_p50 = ">=8hr, P50")), ///
    xline(0, lcolor(gs10) lpattern(dash)) ///
    xtitle("IV coefficient (log outage exposure)") ///
    title("Figure S1: IV Second-Stage Estimates - Alternative Outage Definitions", size(medsmall)) ///
    legend(off) mcolor(navy) ciopts(lcolor(navy)) ///
    graphregion(color(white)) bgcolor(white)
graph export "replication_outputs/figureS1_iv_robustness.png", replace width(1200)


*=============================================================================
* TABLES S6-S8: Correlation Matrices and VIF
*=============================================================================

* Table S6: correlation matrix for main variables
estpost correlate ///
    mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 ///
    poverty_rate median_income_1k avg_cmr_clean ///
    pct_medicaid pct_male psych_w, ///
    matrix listwise
esttab using "replication_outputs/TableS6_correlation.rtf", ///
    unstack not noobs compress ///
    star(* 0.05 ** 0.01 *** 0.001) ///
    title("Table S6: Pairwise Correlation Matrix") ///
    replace

* Table S7: IV correlation matrix
estpost correlate ///
    ln_outage_1hr_p75 dens_substation dens_transmission_line dens_transformer, ///
    matrix listwise
esttab using "replication_outputs/TableS7_correlation_IV.rtf", ///
    unstack not noobs compress ///
    star(* 0.05 ** 0.01 *** 0.001) ///
    title("Table S7: Correlation Matrix for Instrumental Variables") ///
    replace

* Table S8: VIF test (results printed to log)
regress mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_main'
estat vif


********************************************************************************
*                                                                              *
*    PART 2: MENTAL HEALTH AS A SECONDARY (COMORBID) DIAGNOSIS                 *
*            Tables S16-S17 and Note S4               *
*            Data: collapsed_data_with_mh_secondary.dta                        *
*                                                                              *
*  Outcomes:                   
*    mh_secondary_cases   = admissions with a non-MH principal diagnosis but   *
*                           >=1 ICD-10 F-code in DX2-DX101                     *
*    nonmh_nocomorb_cases = admissions with no MH code in any position         *
*    <grp>_mh/_nomh_cases = vulnerable groups (diabetes E08-E13, COPD/asthma   *
*                           J40-J47, heart I20-I25/I50, dementia G30-G31,      *
*                           CKD N18) split by coexisting depression/anxiety    *
*                           (F32/F33/F41/F43) in secondary positions           *
*                                                                              *
*  Headline specification excludes the two Elixhauser comorbidity indices      *
*  (avg_cmr_clean, avg_cmr_mortality_clean) because they are constructed       *
*  from the same secondary diagnosis fields as these outcomes. The with-       *
*  Elixhauser versions are reported as a sensitivity check in the S16 note.    *
********************************************************************************

use "collapsed_data_with_mh_secondary.dta", clear

cap drop psych_w
gen psych_w = psych_per_100k
winsor2 psych_w, replace cuts(0 99) by(AYEAR)

* Model 3 controls, and Model 3 minus the two Elixhauser indices (headline)
local ctrl_main "avg_cmr_clean avg_cmr_mortality_clean poverty_rate pct_male median_income_1k black_pct asian_pct pct_medicaid pct_uninsured psych_w"
local ctrl_red  "poverty_rate pct_male median_income_1k black_pct asian_pct pct_medicaid pct_uninsured psych_w"
local opts      "absorb(zip_num i.AYEAR#i.DQTR) exposure(total_population) vce(cluster zip_num) eform"

*=============================================================================
* TABLE S16: Power outage effects by diagnostic position of mental health
*=============================================================================

* Headline (no-Elixhauser) specification
eststo clear
eststo primary:   ppmlhdfe mh_primary_cases     ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_red', `opts'
eststo secondary: ppmlhdfe mh_secondary_cases   ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_red', `opts'
eststo contrast:  ppmlhdfe nonmh_nocomorb_cases ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_red', `opts'

esttab primary secondary contrast using "replication_outputs/TableS16_diagnostic_position.rtf", replace ///
    eform b(3) se(3) keep(ln_outage_1hr_p75 extreme_heat_days_p90) ///
    title("Table S16: Power outage effects by diagnostic position of mental health") ///
    mtitle("Primary MH" "MH secondary" "No MH code") ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    addnotes("Model 3 controls excluding the two Elixhauser comorbidity indices," ///
             "which are constructed from the same secondary diagnosis fields as" ///
             "the outcomes in columns (2)-(3).")

* Sensitivity: retain the Elixhauser indices (cited in the Table S16 note:
* column 1 IRR = 1.015, column 2 = 0.999, column 3 = 1.004)
eststo clear
eststo primary_f:   ppmlhdfe mh_primary_cases     ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_main', `opts'
eststo secondary_f: ppmlhdfe mh_secondary_cases   ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_main', `opts'
eststo contrast_f:  ppmlhdfe nonmh_nocomorb_cases ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_main', `opts'

esttab primary_f secondary_f contrast_f using "replication_outputs/TableS16_withElixhauser_sensitivity.rtf", replace ///
    eform b(3) se(3) keep(ln_outage_1hr_p75 extreme_heat_days_p90) ///
    title("Table S16 sensitivity: full Model 3 control set (with Elixhauser indices)") ///
    mtitle("Primary MH" "MH secondary" "No MH code") ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01)

*=============================================================================
* NOTE S4 / TABLE S16 NOTE: 24 alternative outage definitions
* Duration thresholds 1-8 hours x severity percentiles P50/P75/P90,
* outcome = mh_secondary_cases, headline (no-Elixhauser) controls.
* Reported result: IRRs range 0.996-1.001, none significant.
*=============================================================================

foreach p in 50 75 90 {
    foreach h in 1 2 3 4 5 6 7 8 {
        cap confirm variable new_`h'hr_episodes_p`p'
        if _rc == 0 {
            cap drop ln_out_`h'hr_p`p'
            gen ln_out_`h'hr_p`p' = ln(new_`h'hr_episodes_p`p' + 1)
        }
        else di as error ">>> MISSING variable: new_`h'hr_episodes_p`p'"
    }
}

local irr_min = 99
local irr_max = -99
di _newline "==== 24-definition grid, mh_secondary_cases (no-Elixhauser controls) ===="
foreach p in 50 75 90 {
    foreach h in 1 2 3 4 5 6 7 8 {
        cap confirm variable ln_out_`h'hr_p`p'
        if _rc == 0 {
            quietly ppmlhdfe mh_secondary_cases ln_out_`h'hr_p`p' ///
                extreme_heat_days_p90 `ctrl_red', `opts'
            local b = exp(_b[ln_out_`h'hr_p`p'])
            local pval = 2*(1-normal(abs(_b[ln_out_`h'hr_p`p']/_se[ln_out_`h'hr_p`p'])))
            di "  `h'hr x p`p':  IRR=" %6.4f `b' "   p=" %6.4f `pval'
            if `b' < `irr_min' local irr_min = `b'
            if `b' > `irr_max' local irr_max = `b'
        }
    }
}
di "IRR RANGE ACROSS ALL 24 DEFINITIONS: " %6.4f `irr_min' " to " %6.4f `irr_max'
di "=========================================================================="

*=============================================================================
* TABLE S17: Power outage effects in vulnerable populations,
*            by coexisting depression or anxiety
*=============================================================================

eststo clear
foreach g in diab copd cvd dem ckd {
    eststo `g'_mh:   ppmlhdfe `g'_mh_cases   ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_red', `opts'
    eststo `g'_nomh: ppmlhdfe `g'_nomh_cases ln_outage_1hr_p75 extreme_heat_days_p90 `ctrl_red', `opts'
}

esttab diab_mh diab_nomh copd_mh copd_nomh cvd_mh cvd_nomh dem_mh dem_nomh ckd_mh ckd_nomh ///
    using "replication_outputs/TableS17_vulnerable_groups.rtf", replace ///
    eform b(3) se(3) keep(ln_outage_1hr_p75) ///
    title("Table S17: Power outage effects in vulnerable populations, by coexisting depression or anxiety") ///
    mtitle("Diab+D/A" "Diab-D/A" "COPD+D/A" "COPD-D/A" "Heart+D/A" "Heart-D/A" "Dem+D/A" "Dem-D/A" "CKD+D/A" "CKD-D/A") ///
    stats(N N_clust r2_p, labels("Observations" "Clusters" "Pseudo R2") fmt(0 0 3)) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    addnotes("Main specification excluding the two Elixhauser comorbidity indices (Note S4)." ///
             "+D/A: coexisting depression or anxiety (F32, F33, F41, F43) in secondary" ///
             "diagnosis fields; -D/A: the remainder.")


********************************************************************************
*                                                                              *
*    PART 3: Balanced-ZIP poverty sensitivity check                 *
*    Fixed ZIP-level poverty classification          *
*    (study-period mean poverty rate >= 20%), lower-poverty group randomly     *
*    down-sampled to the higher-poverty ZIP count over 500 replications.       *
*    Slow (~500 x 2 PPML models); enable via run_balanced_zip = 1 above.       *
*                                                                              *
********************************************************************************

if `run_balanced_zip' == 1 {

    use "collapsed_data_for_regression.dta", clear

    cap drop psych_w
    gen psych_w = psych_per_100k
    winsor2 psych_w, replace cuts(0 99) by(AYEAR)

    * Fixed ZIP-level poverty status
    bysort zip_num: egen zip_poverty_mean = mean(poverty_rate)
    gen byte high_pov_zip = (zip_poverty_mean >= 20) if !missing(zip_poverty_mean)

    * Same controls as the manuscript poverty subgroup analysis
    local controls pct_alcohol pct_drug_abuse avg_cmr_clean ///
                   unemployment_rate pct_medicaid pct_uninsured psych_w
    local bopts absorb(zip_num i.AYEAR#i.DQTR) exposure(total_population) vce(cluster zip_num)

    * Mutually exclusive ZIP counts + sampling frame
    preserve
        keep zip_num high_pov_zip
        bysort zip_num: keep if _n == 1
        drop if missing(high_pov_zip)
        count if high_pov_zip == 1
        local C_high = r(N)
        count if high_pov_zip == 0
        local C_low = r(N)
        tempfile ziplist
        save `ziplist'
    restore
    di ">>> Fixed ZIP-level: High-poverty ZIPs = `C_high'; Low-poverty ZIPs = `C_low'"

    * Full-sample estimates under the fixed classification
    ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 ///
        `controls' if high_pov_zip == 0, `bopts' eform
    di "Low-pov  (fixed-ZIP, full): IRR=" %5.3f exp(_b[ln_outage_1hr_p75])
    ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 ///
        `controls' if high_pov_zip == 1, `bopts' eform
    di "High-pov (fixed-ZIP, full): IRR=" %5.3f exp(_b[ln_outage_1hr_p75])

    * Balanced down-sampling: low-poverty ZIPs -> same count as high-poverty ZIPs
    tempname M
    postfile `M' iter b_low se_low p_interact using "replication_outputs/balanced_fixedzip_results.dta", replace

    set seed 12345
    forvalues r = 1/500 {
        preserve
            use `ziplist', clear
            keep if high_pov_zip == 0
            gen double u = runiform()
            sort u
            keep in 1/`C_high'
            gen byte selected_low = 1
            keep zip_num selected_low
            tempfile sellow
            save `sellow'
        restore

        preserve
            merge m:1 zip_num using `sellow', keep(master match) nogen
            replace selected_low = 0 if missing(selected_low)
            gen byte balanced_sample = (high_pov_zip == 1) | (high_pov_zip == 0 & selected_low == 1)

            cap ppmlhdfe mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 ///
                `controls' if high_pov_zip == 0 & selected_low == 1, `bopts'
            if _rc == 0 {
                scalar b_l  = _b[ln_outage_1hr_p75]
                scalar se_l = _se[ln_outage_1hr_p75]
                scalar pint = .
                cap ppmlhdfe mh_primary_cases c.ln_outage_1hr_p75##i.high_pov_zip ///
                    extreme_heat_days_p90 `controls' if balanced_sample == 1, `bopts'
                if _rc == 0 {
                    cap testparm 1.high_pov_zip#c.ln_outage_1hr_p75
                    if _rc == 0 scalar pint = r(p)
                }
                post `M' (`r') (b_l) (se_l) (pint)
            }
        restore
    }
    postclose `M'

    * Summary (Table R1)
    use "replication_outputs/balanced_fixedzip_results.dta", clear
    gen irr_low = exp(b_low)
    gen byte same_dir = (irr_low > 1)
    gen byte sig_low = (2*(1-normal(abs(b_low/se_low))) < 0.05)
    gen byte no_interaction = (p_interact > 0.05)

    di _newline "============ BALANCED-ZIP SUMMARY (R = " _N " usable draws) ============"
    summ irr_low, detail
    centile irr_low, centile(2.5 50 97.5)
    summ same_dir sig_low no_interaction
    di "======================================================================="
}


********************************************************************************
*                                                                              *
*    APPENDIX: Construction of the secondary-diagnosis outcomes     *
*    from the restricted HCUP Maryland SID individual-level file.              *
*    Requires md_sid_complete_merged_data.dta (licensed, NOT distributable);   *
*    included for transparency of the outcome definitions (STAR Methods and    *
*    Note S4). Enable via build_secondary = 1 above.                           *
*                                                                              *
********************************************************************************

if `build_secondary' == 1 {

    * Load only needed columns from the individual-level SID file
    use ZIP AYEAR DQTR PSTATE I10_DX* using "md_sid_complete_merged_data.dta", clear

    * Maryland filter (same as the main data preparation)
    cap confirm variable PSTATE
    if _rc == 0 keep if PSTATE == "MD"

    * Primary MH flag: any ICD-10 F-code as principal diagnosis (DX1)
    gen byte prim_mh = regexm(upper(I10_DX1), "^F")

    * Secondary-position flags over DX2-DX101
    gen byte comorb_any_mh  = 0     // any F-code in a secondary position
    gen byte comorb_depanx  = 0     // depression/anxiety (F32 F33 F41 F43)
    forvalues _k = 2/101 {
        capture confirm string variable I10_DX`_k'
        if _rc continue
        replace comorb_any_mh = 1 if regexm(upper(I10_DX`_k'), "^F")
        replace comorb_depanx = 1 if regexm(upper(I10_DX`_k'), "^F32|^F33|^F41|^F43")
    }

    * Outcomes for Table S16
    gen byte mh_secondary    = (prim_mh == 0) & (comorb_any_mh == 1)
    gen byte nonmh_nocomorb  = (prim_mh == 0) & (comorb_any_mh == 0)

    * Vulnerable groups by principal diagnosis (Table S17)
    gen byte vp_diab = regexm(upper(I10_DX1), "^E0[89]|^E1[0-3]")
    gen byte vp_copd = regexm(upper(I10_DX1), "^J4[0-7]")
    gen byte vp_cvd  = regexm(upper(I10_DX1), "^I2[0-5]|^I50")
    gen byte vp_dem  = regexm(upper(I10_DX1), "^G3[01]")
    gen byte vp_ckd  = regexm(upper(I10_DX1), "^N18")

    foreach g in diab copd cvd dem ckd {
        gen byte `g'_mh   = vp_`g' & comorb_depanx
        gen byte `g'_nomh = vp_`g' & !comorb_depanx
    }

    * Collapse to ZIP-quarter
    collapse (sum) mh_secondary_cases   = mh_secondary   ///
                   nonmh_nocomorb_cases = nonmh_nocomorb ///
                   diab_mh_cases  = diab_mh   diab_nomh_cases = diab_nomh ///
                   copd_mh_cases  = copd_mh   copd_nomh_cases = copd_nomh ///
                   cvd_mh_cases   = cvd_mh    cvd_nomh_cases  = cvd_nomh  ///
                   dem_mh_cases   = dem_mh    dem_nomh_cases  = dem_nomh  ///
                   ckd_mh_cases   = ckd_mh    ckd_nomh_cases  = ckd_nomh, ///
             by(ZIP AYEAR DQTR)
    save "mh_secondary_byzipq.dta", replace

    * Merge into a copy of the regression panel; unmatched ZIP-quarters get
    * zero counts (no such admissions there)
    use "collapsed_data_for_regression.dta", clear
    merge 1:1 ZIP AYEAR DQTR using "mh_secondary_byzipq.dta", keep(master match)
    foreach v of varlist mh_secondary_cases nonmh_nocomorb_cases ///
        diab_mh_cases diab_nomh_cases copd_mh_cases copd_nomh_cases ///
        cvd_mh_cases cvd_nomh_cases dem_mh_cases dem_nomh_cases ///
        ckd_mh_cases ckd_nomh_cases {
        replace `v' = 0 if _merge == 1
    }
    drop _merge
    save "collapsed_data_with_mh_secondary.dta", replace
}



