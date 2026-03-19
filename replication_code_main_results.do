********************************************************************************
* Replication Code for Main Results
* "Do Power Outages Impact Mental Health? Empirical Evidence from Maryland"
*
* This file consolidates the Stata code used to produce the main tables
* and figures in the manuscript. Code is organized by manuscript order.
*
* Required packages: ppmlhdfe, estout, coefplot, reghdfe, ivpoisson, winsor2, outreg2
* Required data: collapsed_data_for_regression.dta
*                (with density variables merged from 05_merge_density.do)
*
* Note: Figures 5 and 6 (infrastructure density maps) were produced in ArcGIS Pro.
********************************************************************************

clear all
set more off

cd "/Users/irenefeng/Downloads/workingfolder/Jan11_2026"

* Install required packages
foreach pkg in ppmlhdfe estout coefplot reghdfe ivpoisson2 winsor2 outreg2 {
    cap which `pkg'
    if _rc ssc install `pkg'
}

* Load data
use "collapsed_data_for_regression.dta", clear

* Winsorize healthcare resource variables
cap drop psych_w
gen psych_w = psych_per_100k
winsor2 psych_w, replace cuts(0 99) by(AYEAR)


********************************************************************************
*                                                                              *
*                   SUMMARY STATISTICS (Supplementary Table)                   *
*                                                                              *
********************************************************************************

* Variable labels
label var mh_primary_cases "Mental health hospitalizations"
label var new_1hr_episodes_p75 "Outage episodes >= 1hr (P75)"
label var ln_outage_1hr_p75 "ln(Outage episodes + 1)"
label var extreme_heat_days_p95 "Extreme heat days (>= P95)"
label var extreme_heat_days_p90 "Extreme heat days (>= P90)"
label var avg_cmr_clean "Comorbidity index (readmission)"
label var avg_cmr_mortality_clean "Comorbidity index (mortality)"
label var median_income_1k "Median income ($1,000s)"
label var poverty_rate "Poverty rate (%)"
label var pct_male "Male patients (%)"
label var pct_medicaid "Medicaid (%)"
label var pct_uninsured "Uninsured (%)"
label var psych_per_100k "Psychiatrists per 100,000"
label var total_population "ZIP code population"

* Table 1 (Summary Statistics)
local main_vars "mh_primary_cases new_1hr_episodes_p75 extreme_heat_days_p95 total_population poverty_rate median_income_1k pct_medicaid pct_uninsured psych_per_100k avg_cmr_clean avg_cmr_mortality_clean pct_male"

estpost summarize `main_vars', detail

esttab using "Table1_Summary_Statistics.rtf", replace ///
    cells("count(fmt(0)) mean(fmt(2)) sd(fmt(2)) min(fmt(2)) p50(fmt(2)) max(fmt(2))") ///
    noobs nonumber nomtitle ///
    collabels("N" "Mean" "SD" "Min" "Median" "Max") ///
    title("Table 1: Summary Statistics") ///
    addnotes("Notes: Sample includes Maryland ZIP codes from 2018-2023." ///
             "Unit of observation: ZIP code x year x quarter.")


********************************************************************************
*                                                                              *
*     FIGURE 1: Baseline PPML-HDFE Estimates (5 Model Specifications)          *
*                                                                              *
********************************************************************************

*--- Run Models 1-5 ---

// Model 1: Basic socioeconomic controls + comorbidity indices
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 ///
    poverty_rate median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store model1

// Model 2: + Extreme weather (heat, precipitation, cold at P90)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 extreme_precip_days_p90 extreme_cold_days_p90 ///
    poverty_rate median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store model2

// Model 3: Full demographic model (main specification)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store model3

// Model 4: Alternative heat definition (P75 threshold)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 unusual_heat_days_p75 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store model4

// Model 5: Disease-specific comorbidity indicators
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p95 ///
    pct_psychoses pct_alcohol pct_diabetes_cx pct_hypertension_cx pct_copd ///
    avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store model5

* Export regression table (Supplementary Table S14)
esttab model1 model2 model3 model4 model5 using "Table1_main_results.rtf", replace ///
    eform b(3) se(3) ///
    title("Table 1: PPML Results - Mental Health Hospitalizations and Power Outages") ///
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
graph save "panel_a.gph", replace
graph export "panel_a_outage.png", replace width(1200)

*--- Figure 1 Panel (b): Weather Controls ---

* Extract labels for weather variables
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
graph export "panel_b_weather.png", replace width(1200)

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
graph export "panel_c_socioeconomic.png", replace width(1200)

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
graph export "panel_d_comorbidity.png", replace width(1200)


********************************************************************************
*                                                                              *
*   FIGURE 2: Robustness - Alternative Outage Threshold Definitions            *
*                                                                              *
********************************************************************************

// P50 threshold (1hr, 50th percentile)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p50 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store outage_p50

// P90 threshold (1hr, 90th percentile)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p90 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store outage_p90

// 2hr P50 threshold
ppmlhdfe mh_primary_cases ///
    ln_outage_2hr_p50 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store outage_2hr_p50

* Export table (Supplementary Table S15)
esttab outage_p50 outage_p90 outage_2hr_p50 using "Table2_robustness_outage.rtf", replace ///
    eform b(3) se(3) ///
    title("Table 2: Robustness - Different Outage Thresholds") ///
    mtitle("1hr P50" "1hr P90" "2hr P50") ///
    stats(N N_clust, labels("N" "Clusters") fmt(0 0)) ///
    star(* 0.10 ** 0.05 *** 0.01) compress

* Figure 2: Robustness coefplot
foreach est in outage_p50 outage_p90 outage_2hr_p50 {
    estimates restore `est'
    * Get the outage variable name from the stored estimates
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
graph export "figure_robustness_outage.png", replace width(1200)


********************************************************************************
*                                                                              *
*     FIGURE 3: Sensitivity Test - Duration Thresholds (1hr-8hr, P50)          *
*                                                                              *
********************************************************************************

* Run PPML for each duration threshold
foreach dur in 1hr 2hr 3hr 4hr 5hr 6hr 7hr 8hr {
    ppmlhdfe mh_primary_cases ln_outage_`dur'_p50 extreme_heat_days_p90 ///
        avg_cmr_clean avg_cmr_mortality_clean ///
        poverty_rate pct_male median_income_1k ///
        black_pct asian_pct ///
        pct_medicaid pct_uninsured psych_w ///
        , absorb(zip_num i.AYEAR#i.DQTR) ///
        exposure(total_population) vce(cluster zip_num) eform
    estimates store dur_`dur'
}

* Export table (Supplementary Table S16)
esttab dur_1hr dur_2hr dur_3hr dur_4hr dur_5hr dur_6hr dur_7hr dur_8hr using "Table_sensitivity_duration.rtf", replace ///
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
graph export "Figure_sensitivity_duration.png", replace width(1200)


********************************************************************************
*                                                                              *
*     TABLE 1 (Manuscript): Effect Modification by Extreme Heat                *
*                            (Outage x Heat Interaction)                       *
*                                                                              *
********************************************************************************

* Create mean-centered variables
capture drop ln_outage_c heat_c mean_outage mean_heat
egen double mean_outage = mean(ln_outage_1hr_p75)
egen double mean_heat = mean(extreme_heat_days_p90)
gen double ln_outage_c = ln_outage_1hr_p75 - mean_outage
gen double heat_c = extreme_heat_days_p90 - mean_heat
label var ln_outage_c "Log outage hours (centered)"
label var heat_c "Extreme heat days (centered)"

* Run interaction model
ppmlhdfe mh_primary_cases c.ln_outage_c##c.heat_c ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
eststo heat_interact
test c.ln_outage_c#c.heat_c

* Export Table 1
esttab heat_interact using "heat_interaction.rtf", replace ///
    eform b(4) se(3) ///
    keep(ln_outage_c heat_c c.ln_outage_c#c.heat_c) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(N N_clust, labels("N" "Clusters") fmt(0 0))


********************************************************************************
*                                                                              *
*     TABLE 2 (Manuscript): IV Poisson Estimates                               *
*     3 Instrument Sets: Transformer, Substation+Transmission, All 3           *
*                                                                              *
********************************************************************************

*--- Column (1): Transformer density as sole IV ---

* First stage
reghdfe ln_outage_1hr_p75 ///
    dens_transformer ///
    extreme_heat_days_p90 ///
    pct_psychoses pct_alcohol ///
    poverty_rate pct_male median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(i.AYEAR) cluster(zip_num)
estimates store fs_trans
test dens_transformer
local F_trans = r(F)
estadd scalar FirstStage_F = `F_trans'

* Second stage
ivpoisson gmm mh_primary_cases ///
    extreme_heat_days_p90 ///
    pct_psychoses pct_alcohol ///
    poverty_rate pct_male median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    i.AYEAR ///
    (ln_outage_1hr_p75 = dens_transformer) ///
    , exposure(total_population) vce(cluster zip_num)
estimates store ss_trans

*--- Column (2): Substation + Transmission line density as IVs ---

* First stage
reghdfe ln_outage_1hr_p75 ///
    dens_substation dens_transmission_line ///
    extreme_heat_days_p90 ///
    pct_psychoses pct_alcohol ///
    poverty_rate pct_male median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(i.AYEAR) cluster(zip_num)
estimates store fs_sub_trans
test dens_substation dens_transmission_line
local F_sub_trans = r(F)
estadd scalar FirstStage_F = `F_sub_trans'

* Second stage
ivpoisson gmm mh_primary_cases ///
    extreme_heat_days_p90 ///
    pct_psychoses pct_alcohol ///
    poverty_rate pct_male median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    i.AYEAR ///
    (ln_outage_1hr_p75 = dens_substation dens_transmission_line) ///
    , exposure(total_population) vce(cluster zip_num)
estimates store ss_sub_trans
estat overid
local hansen_sub_trans = chi2tail(e(J_df), e(J))
estadd scalar hansen_p = `hansen_sub_trans'

*--- Column (3): All 3 IVs ---

* First stage
reghdfe ln_outage_1hr_p75 ///
    dens_transformer dens_transmission_line dens_substation ///
    extreme_heat_days_p90 ///
    pct_psychoses pct_alcohol ///
    poverty_rate pct_male median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(i.AYEAR) cluster(zip_num)
estimates store fs_all3
test dens_transformer dens_transmission_line dens_substation
local F_all3 = r(F)
estadd scalar FirstStage_F = `F_all3'

* Second stage
ivpoisson gmm mh_primary_cases ///
    extreme_heat_days_p90 ///
    pct_psychoses pct_alcohol ///
    poverty_rate pct_male median_income_1k ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    i.AYEAR ///
    (ln_outage_1hr_p75 = dens_transformer dens_transmission_line dens_substation) ///
    , exposure(total_population) vce(cluster zip_num)
estimates store ss_all3
estat overid
local hansen_all3 = chi2tail(e(J_df), e(J))
estadd scalar hansen_p = `hansen_all3'

* Export Table 2 (Second Stage)
esttab ss_trans ss_sub_trans ss_all3 using "Table_IV_comparison.rtf", replace ///
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
esttab fs_trans fs_sub_trans fs_all3 using "Table_IV_firststage_comparison.rtf", replace ///
    b(3) se(3) ///
    title("Table 2 Panel B: First Stage Results") ///
    mtitle("Transformer Only" "Substation + Trans" "All 3 IVs") ///
    keep(dens_transformer dens_transmission_line dens_substation) ///
    stats(N FirstStage_F, labels("Observations" "F-statistic") fmt(0 2)) ///
    star(* 0.10 ** 0.05 *** 0.01)


********************************************************************************
*                                                                              *
*     FIGURE 4: Heterogeneity Analysis                                         *
*     (Urban/Rural, Poverty, Medicaid, Season, Gender)                         *
*                                                                              *
********************************************************************************

*--- 1. Urban vs Rural (Large Metro RUCC=1 vs Rest) ---

cap drop large_metro
gen large_metro = (PL_RUCC == 1) if !missing(PL_RUCC)

* Large metro
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if large_metro == 1 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store large_metro_grp

* Non-large-metro
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if large_metro == 0 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
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
estimates store interact_large_metro

* Z-test for subgroup difference
estimates restore large_metro_grp
local b1 = _b[ln_outage_1hr_p75]
local se1 = _se[ln_outage_1hr_p75]
estimates restore non_large_metro_grp
local b2 = _b[ln_outage_1hr_p75]
local se2 = _se[ln_outage_1hr_p75]
local z_diff = (`b1' - `b2') / sqrt(`se1'^2 + `se2'^2)
local p_diff = 2*(1-normal(abs(`z_diff')))
di "Large Metro vs Rest: z=" %6.3f `z_diff' "  p=" %6.4f `p_diff'

*--- 2. Poverty Level (Census Bureau 20% threshold) ---

cap drop high_poverty_census
gen high_poverty_census = (poverty_rate >= 20) if !missing(poverty_rate)

* Low poverty (<20%)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_poverty_census == 0 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store low_pov_census

* High poverty (>=20%)
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_poverty_census == 1 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store high_pov_census

* Interaction test
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
estimates store interact_pov_census

*--- 3. Medicaid Coverage (25% threshold) ---

cap drop high_medicaid_25
gen high_medicaid_25 = (pct_medicaid >= 25) if !missing(pct_medicaid)

* Low Medicaid
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    median_income_1k psych_w ///
    if high_medicaid_25 == 0 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store low_medicaid

* High Medicaid
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    median_income_1k psych_w ///
    if high_medicaid_25 == 1 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store high_medicaid_grp

* Interaction test
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
estimates store interact_medicaid

*--- 4. Season (Warm Q2-Q3 vs Cold Q1&Q4) ---

cap drop warm_season
gen warm_season = (DQTR == 2 | DQTR == 3) if !missing(DQTR)

* Cold season
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if warm_season == 0 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store cold_season

* Warm season
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if warm_season == 1 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store warm_season_grp

* Interaction test
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 c.ln_outage_1hr_p75#i.warm_season ///
    i.warm_season ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if !missing(warm_season) ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num)
estimates store interact_season

*--- 5. Gender Composition (median split) ---

cap drop high_male
sum pct_male, detail
gen high_male = (pct_male > r(p50)) if !missing(pct_male)

* Low male %
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_male == 0 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store low_male

* High male %
ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    pct_alcohol pct_drug_abuse avg_cmr_clean ///
    poverty_rate unemployment_rate ///
    pct_medicaid pct_uninsured psych_w ///
    if high_male == 1 ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store high_male

* Interaction test
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
estimates store interact_gender

*--- Export heterogeneity tables ---

outreg2 [large_metro_grp non_large_metro_grp low_pov_census high_pov_census] ///
    using "Table1_hetero_partA_results.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Table 4: Heterogeneity Analysis - Part A") ///
    ctitle("Large Metro", "Non-Large Metro", "Low Poverty", "High Poverty") ///
    addtext(ZIP FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at ZIP code level.")

outreg2 [low_medicaid high_medicaid_grp cold_season warm_season_grp] ///
    using "Table1_hetero_partB_results.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Table 4: Heterogeneity Analysis - Part B") ///
    ctitle("Low Medicaid", "High Medicaid", "Cold Season", "Warm Season") ///
    addtext(ZIP FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at ZIP code level.")

outreg2 [low_male high_male] ///
    using "Table1_hetero_partC_results.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("Table 4: Heterogeneity Analysis - Part C") ///
    ctitle("Lower Male %", "Higher Male %") ///
    addtext(ZIP FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at ZIP code level.")

*--- Figure 4: Heterogeneity coefplot ---

* Extract labels with significance stars
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
graph export "figure_heterogeneity.png", replace width(1600)


********************************************************************************
*                                                                              *
*     SUPPLEMENTARY: City-by-Year FE Robustness (Table S17)                    *
*                                                                              *
********************************************************************************

ppmlhdfe mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured psych_w ///
    , absorb(city_id#i.AYEAR i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster city_id) eform
estimates store m2_city_year_fe

outreg2 [m2_city_year_fe] using "Table_cityFE_results.doc", replace eform se ///
    bdec(3) sdec(3) ///
    title("PPML Results - City x Year FE Robustness") ///
    ctitle("City x Year FE + Year x Quarter FE") ///
    addtext(City x Year FE, YES, Year x Quarter FE, YES, Exposure, Population) ///
    addnote("Standard errors clustered at city level.")


********************************************************************************
*                                                                              *
*     SUPPLEMENTARY: Raw Count Specification (Table S19)                       *
*                                                                              *
********************************************************************************

ppmlhdfe mh_primary_cases ///
    new_1hr_episodes_p75 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured ///
    psych_w ///
    , absorb(zip_num i.AYEAR#i.DQTR) ///
    exposure(total_population) vce(cluster zip_num) eform
estimates store model_raw_count


********************************************************************************
*                                                                              *
*     SUPPLEMENTARY: IV Robustness - Alternative Outage Definitions            *
*     (Supplementary Figure N1)                                                *
*                                                                              *
********************************************************************************

* IV estimates across alternative outage definitions (all using 3 IVs)
foreach dur_var in ln_outage_1hr_p75 ln_outage_2hr_p75 ln_outage_3hr_p75 ///
                   ln_outage_4hr_p75 ln_outage_5hr_p50 ln_outage_6hr_p50 ln_outage_8hr_p50 {

    local i = 0
    local i = `i' + 1

    * First stage
    reghdfe `dur_var' ///
        dens_transformer dens_transmission_line dens_substation ///
        extreme_heat_days_p90 ///
        pct_psychoses pct_alcohol ///
        poverty_rate pct_male median_income_1k ///
        avg_cmr_clean avg_cmr_mortality_clean ///
        black_pct asian_pct ///
        pct_medicaid pct_uninsured ///
        psych_w ///
        , absorb(i.AYEAR) cluster(zip_num)
    test dens_transformer dens_transmission_line dens_substation

    * Second stage
    ivpoisson gmm mh_primary_cases ///
        extreme_heat_days_p90 ///
        pct_psychoses pct_alcohol ///
        poverty_rate pct_male median_income_1k ///
        avg_cmr_clean avg_cmr_mortality_clean ///
        black_pct asian_pct ///
        pct_medicaid pct_uninsured ///
        psych_w ///
        i.AYEAR ///
        (`dur_var' = dens_transformer dens_transmission_line dens_substation) ///
        , exposure(total_population) vce(cluster zip_num)
    estimates store iv_`dur_var'
    estat overid
}


********************************************************************************
*                                                                              *
*     SUPPLEMENTARY: Correlation Matrices and VIF                              *
*                                                                              *
********************************************************************************

* Correlation matrix for main variables
estpost correlate ///
    mh_primary_cases ln_outage_1hr_p75 extreme_heat_days_p90 ///
    poverty_rate median_income_1k avg_cmr_clean ///
    pct_medicaid pct_male psych_w, ///
    matrix listwise

esttab using "correlation_table.rtf", ///
    unstack not noobs compress ///
    star(* 0.05 ** 0.01 *** 0.001) ///
    title("Pairwise Correlation Matrix") ///
    replace

* IV correlation matrix
estpost correlate ///
    ln_outage_1hr_p75 dens_substation dens_transmission_line dens_transformer, ///
    matrix listwise

esttab using "correlation_IV.rtf", ///
    unstack not noobs compress ///
    star(* 0.05 ** 0.01 *** 0.001) ///
    title("Supplementary Table S11: Correlation Matrix for Instrumental Variables") ///
    replace

* VIF test
regress mh_primary_cases ///
    ln_outage_1hr_p75 extreme_heat_days_p90 ///
    avg_cmr_clean avg_cmr_mortality_clean ///
    poverty_rate pct_male median_income_1k ///
    black_pct asian_pct ///
    pct_medicaid pct_uninsured psych_w
estat vif


********************************************************************************
di "Replication code complete."
di "Note: Figures 5 and 6 (infrastructure density maps) were produced in ArcGIS Pro."
********************************************************************************
