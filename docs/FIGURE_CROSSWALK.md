# Figure script crosswalk

The current manuscript has five main figures and thirteen supplementary figures. This repository releases their plotting code, not the empirical plotting datasets or rendered images. The older `SF` identifiers in script names reflect the original plotting bundle, not current journal numbering.

| Current figure | Final plotting script | Input status |
| --- | --- | --- |
| Figure 1 | `figures/main/scripts/Figure1.R` | Diagram inputs/configuration must be supplied locally. |
| Figure 2 | `figures/main/scripts/Figure2.R` | Restricted aggregate sources required. |
| Figure 3 | `figures/main/scripts/Figure3.R` | Restricted aggregate sources required; `prepare_Figure3_shape_summary.R` creates its grouped shape panel. |
| Figure 4 | `figures/D3/scripts/Figure4.R` | Unchanged plotting calculations: A, label fidelity; B, corrected hospital funnel; C, corrected Base + K2 and completed corrected Raw-EN/Raw-RF, retained Base and source-internal CV; D, corrected Base + K2/retained Base held-out mortality recalibration. Raw completion 04 stages a complete source directory with new panel C and unchanged A/B/D. |
| Figure 5 | `figures/main/scripts/Figure5.R` | Restricted aggregate sources required; the published presentation omits individual weight extrema. |
| Supplementary Figure S1 | `figures/main/scripts/SF_RiskSet.R` | `build_public_riskset_flow.R` combines pre-landmark RRT exclusions; the small component remains restricted. |
| Supplementary Figure S2 | `figures/supplement/R/figures/supplement/SF8_oracle_accuracy_dose_response.R` | Restricted aggregate sources required. |
| Supplementary Figure S3 | `figures/supplement/R/figures/supplement/SF6_domain4_gate_overlap_v4_countpanels.R` | Restricted aggregate sources required. |
| Supplementary Figure S4 | `figures/supplement/R/figures/supplement/SF7_D1_split_scenario_ids_v4_countpanels.R` | Restricted aggregate sources required. |
| Supplementary Figure S5 | `figures/supplement/R/figures/supplement/SF12_D3_scenario_ids.R` | Historical L2 control. |
| Supplementary Figure S6 | `figures/supplement/R/figures/supplement/SF1_mice_diagnostics_v4_countpanels.R` | Screened aggregate density curves required; individual density input is excluded. |
| Supplementary Figure S7 | `figures/supplement/R/figures/supplement/SF2_phenotype_harmonization.R` | Historical reconstructed-rule profile. |
| Supplementary Figure S8 | `figures/supplement/R/figures/supplement/SF3_D1_lineage_fisher_ecdf_v5.R` | Aggregate plotting input required. |
| Supplementary Figure S9 | `figures/supplement/R/figures/supplement/SF4_domain2_label_lineage_v4_countpanels.R` | Aggregate plotting input required. |
| Supplementary Figure S10 | `figures/supplement/R/figures/supplement/SF5_external_validation_sensitivity.R` | Historical L2 external sensitivity. |
| Supplementary Figure S11 | `figures/supplement/R/figures/supplement/SF9_subtype_separation_reproducibility.R` | Aggregate plotting input required. |
| Supplementary Figure S12 | `figures/supplement/R/figures/supplement/SF10_missingness_scenario_ids.R` | Aggregate plotting input required. |
| Supplementary Figure S13 | `figures/supplement/R/figures/supplement/SF11_algorithm_scenario_ids.R` | Aggregate plotting input required. |

No small-hospital row, patient-level histogram bin, or ungrouped pre-landmark RRT exclusion is distributed here. The current Figure 4 script draws 96 eligible hospitals individually and 103 in grouped markers; this is not the older rank-plot design.

Date: 2026-10-08. Only one Figure 4 renderer is retained. Correction 08 is the earlier Base + K2-only updater and requires existing fidelity, internal-performance and reference-range inputs. Raw completion 01 provides the new panel-C transport source, 03 finalizes historical comparisons without changing corrected estimates, 02 verifies completion, and 04 stages all eleven renderer inputs by byte-preserving copy.

Main sources and generated figures/QA are external. The shared adapter is `R/common/02_figure_paths.R`; main sources default to `EHR_AUDIT_OUTPUT_ROOT/main_figures/source_data/`. Current Figure 4 sources default to `EHR_AUDIT_OUTPUT_ROOT/domain3_raw_temperature_completion/figure4/source_data/`. `EHR_AUDIT_FIGURE4_SOURCE_ROOT` overrides that default; one explicit directory argument overrides the environment in both run_all and Figure4. Correction-stage sources use `EHR_AUDIT_CORRECTION_FIGURE_SOURCE_ROOT` and remain separate.

Every main renderer records actual consumed CSV paths/checksums outside the repository. The total entry checks them against its unchanged selected-source snapshot, including source membership and the distributed main style specification. The final figure was rendered in the parent workflow with R 4.5.3; this packaging did not render or inspect pixels. Raw model recovery/evaluation used R 4.6.0 separately. No empirical source CSV or pixels are distributed.
