# Plasmix proteomics analysis scripts

Run every script from the `01_integration` project root. The scripts use paths relative to that root and expect `data/`, `results/`, `figures/`, `tables/`, `cache/` and `utils/` to exist. They do not create directories automatically.

By default, internal inputs are read from the sibling `00_data/` directory. Set `PLASMIX_RELEASE_ROOT` or `PLASMIX_UPSTREAM_DATA_ROOT` only when the project or upstream data are stored elsewhere.

## Execution order

1. `export_software_environment.R` records the R, Bioconductor, package and optional Python/FreeSASA environment.
2. `prepare_profiles_and_detection.R` creates the complete assay-level public profile and feature metadata. Figure 1 summarizes protein coverage from the mapped assay records. For protein-level downstream analyses, multiple SOMAmer measurements assigned to the same single accession are averaged on the log2 scale; named analytes from the other platforms remain distinct. Unmapped and multiple-accession features are excluded through `utils/feature_mapping.R`.
3. `prepare_physchem_annotation.R` and `prepare_differential_expression.R` create the reusable physicochemical, differential-expression and MAPD results.
4. `fig1_profiles_overview.R` and `fig2_titration_benchmark.R` create the first two main figures and the titration/CV results.
5. `ed1_titration_tolerance.R` and `ed3_som_reshaping_compact.R` use the Figure 2 outputs.
6. `ed4_external_sex_consistency.R` creates the standardized external-cohort sex-effect tables; `ed5_plasmix_sex_gradient.R` uses them.
7. `ed6_mfnp_contrast_difference.R` creates the cross-platform differential-expression consensus used by `fig5_sample_reference_ratio.R`.
8. `fig3_cross_setting_concordance.R` uses the Figure 2 and external-cohort results.
9. `fig4_distortion_cause.R` creates the model objects used by `ed7_distortion_features.R`.
10. `fig6_integration_guidance.R` creates the reusable integration object used by `ed8_reference_background_integration.R` and `ed9_abundance_integration_outcomes.R`.
11. `ed2_normalization_detection.R` reports assay-level normalization factors and stage-specific detection and may run any time after step 2.

Source-data workbooks are written only to `tables/`. Reusable computational outputs and caches are written to `results/` or `cache/`.

Extended Data figure numbers in script names, output filenames and source-data sheet labels follow the manuscript order. Execution follows the dependencies above, not numerical order.
