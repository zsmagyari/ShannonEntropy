# Dynamic World Shannon Entropy as a Scale-Sensitive Indicator of Surface Urban Heat Island Intensity: Evidence from Seven Romanian Cities
*Authored by Zsolt Magyari-Sáska and Ionel Haidu*

The repository has two main folders:

- **scripts** – contains the R scripts corresponding to each methodological step.  
- **results** – contains the principal output files, organized according to the structure described below.  

---

### `RO_<CITY_KEY>_sensitivity_results_LONG_2021-2025.csv`
The files in this folder contain the full long-format results of the sensitivity analysis, documenting for each tested parameter combination and window size the model-fitting outcomes and performance diagnostics used to evaluate the SUHI–entropy relationship.
Each row corresponds to one *city × scenario × entropy window* combination.

*Column Description*

- identifiers and analysis context (`city`, `interval`, `scenario_id`, `window`)

- urban–rural delineation settings (`TH_URBAN`, `TH_RURAL`, `BUFFER_M`, `K_RURAL_MAX`, `RURAL_MAX_M`, `ELEV_TOL_M`)

- sample sizes (`n_urban_px`, `n_rural_px`, `n_points`)

- thermal summaries (`lst_urban_median`, `lst_rural_median`, `suhi_delta_median_C`)

- correlation metrics (`spearman_rho`, `pearson_r`)

- candidate-model AIC values (`aic_null`, `aic_lin`, `aic_quad`, `aic_ns3`)

- candidate-model AICc values (`aicc_null`, `aicc_lin`, `aicc_quad`, `aicc_ns3`)

- best-model selection metrics (`best_model`, `best_aicc`, `delta_aicc_vs_null`, `best_r2`)


### `RO_<CITY_KEY>_main_GAM_summary_2021-2025.csv`
It contains GAM results across all selected entropy windows using:

- the fixed best-supported urban-side setting,

- the chosen reference rural setting,

- and a common-point comparison framework across windows.

Each row represents one tested entropy window.
Main fields include:
- analysis settings (`TH_CORE`, `TH_WATER`, `TH_URBAN`, `BUFFER_M`, `TH_RURAL`, `K_RURAL_MAX`, `ELEV_TOL_M`)

- sample sizes (`n_urban\_px`, `n_rural_px`, `n_points`)

- thermal summaries (`lst_urban_median`, `lst_rural_median`, `suhi_delta_median_C`)

- correlation metrics (`spearman_rho`, `pearson_r`)

- variability indicators for entropy and SUHI (`entropy_iqr`, `entropy_sd`, `entropy_range`, `suhi_iqr`, `suhi_sd`, `suhi_range`)

- GAM1 metrics (`gam1_aic`, `gam1_aicc`, `gam1_dev_expl`, `gam1_r2_adj`, `gam1_edf_entropy`, `gam1_p_smooth_entropy`)

- GAM2 metrics (`gam2_aic`, `gam2_aicc`, `gam2_dev_expl`, `gam2_r2_adj`, `gam2_edf_entropy`, `gam2_edf_xy`, `gam2_p_smooth_entropy`, `gam2_p_smooth_xy`)

- spatial smooth settings (`k_xy_used`, `k_xy_mode`, `k_xy_auto`)

- city-level ranking fields (`delta_gam2_aicc_city`, `delta_gam1_aicc_city`, `rank_gam2_aicc_city`, `rank_abs_rho_city`, `rank_gam2_devexpl_city`, `is_gam2_best_city`)


### `RO_<CITY_KEY>_robustness_GAM_summary_2021-2025.csv`

Contains the robustness-analysis GAM results for alternative tied-best rural delineations while keeping the urban-side setting fixed.
Each row represents one *entropy window × alternative tied-best rural setting* combination.
Its structure is the same as that of the main GAM summary table.


