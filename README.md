# HCC Vascular Invasion Continuum — Analysis Code & Data

Reproducible analysis code and derived data for the manuscript
*"Integrative molecular dissection of the vascular invasion continuum in
hepatocellular carcinoma"* (manuscript in preparation).

> **Note on scope.** This repository contains the **analysis source code** and the
> **derived data tables** needed to verify every reported number, table, and figure.
> It does **not** contain the manuscript draft (`.docx` / `.md`) — the manuscript is
> authored and version-controlled separately to avoid bundling multiple draft
> versions. The table/figure generation script below runs against the derived data only.

## What is in this repository

```
analysis_repo/
├── 01_data_acquisition/   # scripts that fetch raw data from GEO / TCGA / cBioPortal
├── 02_analysis/           # main analysis pipeline (_step1 → _step4)
├── 03_figures/            # figure-generation scripts
├── 04_reproducibility/    # table/figure generation scripts
├── utils/                 # optional manuscript-QA helpers (not required for science)
├── data/tcga/             # small derived cBioPortal JSON (3 files)
├── results/               # derived tables (30 CSV + fig6_manifest.json)
├── requirements.txt
├── LICENSE
└── README.md
```

## Data sources (raw inputs)

Raw data are **not** committed (size). They are publicly available and can be
re-fetched with the scripts in `01_data_acquisition/`, or downloaded directly:

| Source | Accession |
|--------|-----------|
| GEO (MVI paired) | GSE77509, GSE69164 |
| GEO (MVI external validation) | GSE14520 |
| GEO (single-cell) | GSE149614 |
| GEO (immunotherapy) | GSE202069, GSE36376 |
| TCGA | TCGA-LIHC (STAR TPM) |
| cBioPortal | LIHC cohort (patient / VI attributes) |

The single-cell processed object `results/gse149614_processed_qc.h5ad`
(1.1 GB) is deposited on Zenodo (see DOI in the manuscript Data availability).

## How to reproduce

```bash
pip install -r requirements.txt
```

**Stage 1 — acquire raw data.** Run the scripts in `01_data_acquisition/`
(each fetches one source), or download the accessions above.

**Stage 2 — analysis pipeline.** Run `02_analysis/_step1_3_batch.py` through
`02_analysis/_step4_2_external.py` in numeric order. Outputs land in `results/`.

**Stage 3 — figures.** Run the scripts in `03_figures/`.

**Stage 4 — table/figure regeneration (runs out-of-the-box, no large downloads).**
From the repository root:

```bash
python 04_reproducibility/generate_tables.py       # regenerates Table 1 & 2 from source data
```

- `generate_tables.py` recomputes Table 1 (clinical characteristics) and Table 2
  (Cox regression) directly from the committed `results/` tables and prints them.
  It also saves `results/table1_clinical.csv` and `results/table2_cox.csv`. (If the
  manuscript file is present locally it can additionally inline the tables into the
  draft; otherwise it only prints and saves the CSVs.)

> `04_reproducibility/_fig6_reorder.py` regenerates the Fig. 6 manifest
> (`results/fig6_manifest.json`) from the single-cell object. It requires
> `results/gse149614_processed_qc.h5ad` (1.1 GB, on Zenodo) and is provided as
> reference code rather than a run-out-of-the-box step.

## Scripts by stage

### 01 — Data acquisition
- `_cbio_attr.py`
- `_cbio_vi.py`
- `_download_gse14520.py`
- `_download_remaining.py`
- `_download_sc.py`
- `_download_sc_fixed.py`
- `_geo_download.py`
- `_geo_download_suppl.py`
- `_geo_download_suppl2.py`
- `_geo_fetch_meta.py`
- `_geo_gse14520.py`
- `_tcga_clinical.py`
- `_tcga_download.py`
- `_tcga_mvi_check.py`
- `_tide_analysis.py`
- `_tide_download_expr.py`

### 02 — Analysis pipeline
- `_step1_3_batch.py`
- `_step2_1_diff.py`
- `_step2_1_enrichment.py`
- `_step2_1_volcano.py`
- `_step2_2_finish.py`
- `_step2_2_score.py`
- `_step2_2_ssgsea.py`
- `_step2_3_cluster.py`
- `_step2_3_heatmap.py`
- `_step2_4_cox_fix.py`
- `_step2_4_survival_immune.py`
- `_step2_4_table1_cox.py`
- `_step3_1_preprocess.py`
- `_step3_1_qc.py`
- `_step3_2_cluster.py`
- `_step3_2_invasive.py`
- `_step3_3_program.py`
- `_step3_3_review.py`
- `_step3_3_trajectory.py`
- `_step4_2_external.py`

### 03 — Figures
- `_fig1_fix.py`
- `_fig2_fix.py`
- `_figs_clinical.py`
- `_figs_core.py`
- `_figs_figs2_fix.py`
- `_figs_supplementary.py`

### 04 — Table / figure generation
- `generate_tables.py` — regenerates Table 1 & 2 from source data (runs without the manuscript).
- `_fig6_reorder.py` — Fig. 6 manifest generator (requires the 1.1 GB h5ad on Zenodo).

### utils — optional manuscript QA (not required for science)
- `_apply_new_numbers.py`
- `_diag_imgs.py`
- `_docx_tables.py`
- `_fix_v5_docx.py`
- `_fix_v5_final.py`
- `_ground_truth.py`
- `_inspect_docx_images.py`
- `_map_docx_images.py`
- `_md_to_docx.py`
- `_pubmed_check.py`
- `_pubmed_check2.py`
- `_pubmed_check3.py`
- `_recompute.py`
- `_refcheck_body.py`
- `_refcheck_consistency.py`
- `_refcheck_extract.py`
- `_resume_count_only.py`
- `_verify_v5fix.py`

## Dependencies

lifelines, matplotlib, numpy, pandas, python-docx, requests, scanpy, scikit-learn, scipy

## License

MIT — see `LICENSE`.

## Citation

Chen Q, et al. Integrative molecular dissection of the vascular invasion
continuum in hepatocellular carcinoma. Manuscript in preparation (2026).
If you use this code or these data, please cite the manuscript above.
