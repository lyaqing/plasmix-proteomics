# Plasmix proteomics integration

This repository contains the R analysis workflows and shared utilities used for the Plasmix plasma-proteomics integration study. The analyses benchmark technical precision, preservation of a predefined plasma titration gradient, quantitative comparability across analytical settings, and reference-based data integration.

## Repository contents

- `scripts/`: data-preparation, main-figure, extended-data-figure, and software-environment workflows.
- `utils/`: shared functions for feature mapping, batch correction, benchmarking, differential analysis, imputation, and figure styling.
- `reference_manifest.tsv`: versions, citations, checksums, access conditions, and analysis roles of external resources.
- `data/`: processed public inputs downloaded from Figshare; data files are not tracked by Git.
- `results/`: reusable derived inputs and analysis outputs; generated files are not tracked by Git.
- `figures/`: generated figures.
- `tables/`: generated source-data workbooks.
- `cache/`: reusable computational caches.

## Public data

The processed dataset is available from Figshare:

<https://doi.org/10.6084/m9.figshare.32797509>

Download the files and place them at the following paths:

| Repository path | Description |
| --- | --- |
| `data/study_metadata.xlsx` | Batch, sample, and analytical-variance metadata. |
| `data/protein_profiles_long.tsv.gz` | Harmonized assay-level abundance profiles across platforms, batches, samples, and processing levels. |
| `data/feature_metadata.tsv.gz` | Assay-feature, target and UniProt annotations used to derive batch-specific analysis sets. |
| `data/physchem_matrix.tsv.gz` | Protein-level physicochemical, structural, localization, and circulating-abundance annotations. |
| `data/physchem_dictionary.tsv` | Definitions, computation rules, quality summaries, and feature-retention decisions for the physicochemical annotations. |
| `results/detection_status.tsv.gz` | Assay-level, feature-by-batch detection summaries; SOM protein-level calls are derived from this file when needed. |

Source Data workbooks are provided with the associated article and are not duplicated in this repository. Raw and third-party inputs are not redistributed here. Their versions, sources, access conditions, and uses are documented in `reference_manifest.tsv`.

## Running the analysis

Run scripts from the repository root. The project expects `data/`, `results/`, `figures/`, `tables/`, `cache/`, `utils/`, and `scripts/` to be present.

The default local layout assumes that non-public or externally downloaded upstream inputs are stored in a sibling directory named `00_data`. Alternative locations can be configured with:

- `PLASMIX_RELEASE_ROOT`: path to this repository.
- `PLASMIX_UPSTREAM_DATA_ROOT`: path to the upstream data directory.

The complete execution order and script dependencies are documented in [`scripts/README.md`](scripts/README.md).

The public long profile retains harmonized assay records, including protein-group measurements. Figure 1 summarizes protein coverage from mapped assay records. For protein-level downstream analyses, multiple SOMAmer measurements assigned to the same single accession are averaged on the log2 scale, while named analytes from the other platforms remain distinct; unmapped and multiple-accession features are excluded from those analyses. Protein-level SOM detection is derived from the assay-level detection file when needed. Assay-level normalization and stage-specific detection remain available in Extended Data Figure 2. Mapping rules are implemented in `utils/feature_mapping.R`.

The public downstream workflow starts from the six Figshare files listed above. Preparation and external-cohort scripts that require raw, vendor-licensed, publisher-hosted, or restricted resources are retained for methodological transparency but cannot be reproduced from the Figshare files alone.

## Software environment

Package availability is checked within the analysis scripts. The following optional command records the local R, Bioconductor, package, and optional Python/FreeSASA environment; it is not required to run the analyses:

```bash
Rscript scripts/export_software_environment.R
```

The resulting environment records are written to `results/`.

## Citation

Please cite the associated article and the Figshare dataset when using the code or processed data. Publication details will be added when available.
