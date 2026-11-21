# Higher cortical morphometric similarity in women with fibromyalgia

This repository contains the analysis code for a cross-sectional secondary study of 33 women with fibromyalgia and 33 healthy control women in OpenNeuro ds004144. Individual morphometric inverse divergence (MIND) networks were estimated from T1-weighted MRI. The workflow covers the primary network comparisons and the available component, regional, sensitivity, clinical, voxelwise, and spatial annotation analyses.

**Scope:** The code includes the PET reference map analyses for the NBS components. The whole-cortex MIND–PET correlation analysis reported in the manuscript is outside the scope of this repository.

### Repository layout

- [`scripts/`](scripts/): workflow entry point, MRI quality control, FreeSurfer processing, and MIND network construction.
- [`analysis/`](analysis/): statistical models, sensitivity analyses, and spatial annotation.
- [`config/`](config/) and [`environment/`](environment/): configuration, reference resources, and software dependencies.
- [`docs/outputs.md`](docs/outputs.md): mapping from generated outputs to manuscript tables and figures.

### System requirements

- Linux, Bash, Git, `sha256sum`, Apptainer, and network access for public resources and TemplateFlow. Use local scratch space and avoid whitespace in paths supplied to MATLAB or container bind commands. Start with 64 GB RAM and 100 GB of free disk space; adjust job counts downward if needed. Full FreeSurfer, MRIQC, and VBM runs can take many hours.
- MRIQC `24.1.0.dev0+gd5b13cb5.d20240826`; FreeSurfer `7.3.2` (used in the fMRIPrep `25.2.5` container); MATLAB `R2026a` Update 4; NBS `1.2`; SPM `25.01.02`; R `4.6.0`, `permuco` `1.1.3`. See `environment/` for library versions. The MIND subject-level stage uses the Python shipped with the FreeSurfer container.
- A valid local FreeSurfer license and MATLAB installation are required. No licensed software, container image, or license key is bundled.

### Installation and inputs

1. Obtain the [OpenNeuro ds004144 snapshot 1.0.2](https://openneuro.org/datasets/ds004144/versions/1.0.2), including the T1w images and both `phenotype/participants.tsv` and `phenotype/clinical.tsv`. Use OpenNeuro's download instructions and preserve the directory structure. Place it at `data/ds004144` or set `BIDS_DIR`. Do not commit dataset contents.
2. Use [MRIQC](https://github.com/nipreps/mriqc), [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/), and the [fMRIPrep distribution](https://github.com/nipreps/fmriprep) to obtain container images with the specified versions. Verify the software versions inside them. Place the images and your FreeSurfer license under `external/` or configure their locations. The MRIQC stage loads the bundled QI2 compatibility patch; a vanilla MRIQC run is not the exact original run configuration. Visual inspection of T1-weighted images and FreeSurfer reconstructions is performed separately from the automated workflow.
3. Install MATLAB. Extract [NBS](https://www.nitrc.org/projects/nbs/) and [SPM](https://www.fil.ion.ucl.ac.uk/spm/software/) into `external/NBS1.2` and `external/spm` or set `NBS_HOME` and `SPM_HOME`. These directories must contain `NBSrun.m` and `spm.m`, respectively.
4. Install the official MIND implementation outside the tracked source tree:

   ```bash
   mkdir -p external
   git clone https://github.com/isebenius/MIND.git external/MIND
   git -C external/MIND checkout --detach 0d334454eaac62e49a197801446ce7750256c882
   ```

5. Install R and the recorded R packages. For example, use CRAN source archives with `install.packages(..., repos = NULL, type = "source")` for the versions in `environment/R-session.txt`; install `permute` and `Rcpp` before `permuco`. Use a clean R library and specify `R_LIBS` if necessary.
6. Install the Python packages in separate environments. The original Neurosynth and spin analyses used Python `3.8.7rc1` on Windows. Use Python 3.8 for these environments. Numerical equivalence across platforms and Python builds has not been established.

   ```bash
   python3.12 -m venv .venv
   .venv/bin/python -m pip install -r environment/requirements.txt
   python3.8 -m venv .venv-nimare
   .venv-nimare/bin/python -m pip install -r environment/requirements-neurosynth.txt
   python3.8 -m venv .venv-spins
   .venv-spins/bin/python -m pip install -r environment/requirements-spins.txt
   cp config/local.example.sh config/local.sh
   ```

   In `config/local.sh`, set `PYTHON_BIN=.venv/bin/python`, `NEUROSYNTH_PYTHON=.venv-nimare/bin/python`, and `SPIN_PYTHON=.venv-spins/bin/python`. Set any non-default paths there. Paths in the example are relative to the repository root. Set `FM_WORK` to a dedicated directory for generated outputs.
7. Download the [Neurosynth v7 database](https://github.com/neurosynth/neurosynth-data) (e.g. using NiMARE's `fetch_neurosynth` for version 7). Place the four inputs under `external/neurosynth-v7/`: `data-neurosynth_version-7_coordinates.tsv.gz`, `data-neurosynth_version-7_metadata.tsv.gz`, `data-neurosynth_version-7_vocab-terms_source-abstract_type-tfidf_features.npz`, and `data-neurosynth_version-7_vocab-terms_vocabulary.txt`. CBIG atlas assets and the Hansen PET reference maps are fetched by the relevant steps. Schaefer 200 surface files and the MNI atlas are checked against the recorded input hashes. The VBM step needs the spatial step's component maps for its descriptive overlap crosswalk.

### Run in order

```bash
bash scripts/run.sh help
for stage in qc freesurfer atlas mind primary nbs reference nodes robustness clinical surface parcellation maps terms decode spatial vbm-preprocess vbm; do
  bash scripts/run.sh "$stage" || break
done
```

Runs write to `FM_WORK` (default `work/`). To start from already generated data, set `FS_DIR` and `MRIQC_OUT` to the matching existing derivatives; then start at `atlas`. Check that the two derivatives match the original BIDS subject IDs and software. To avoid reusing stale partial artifacts, use a fresh `FM_WORK` for independent runs. Approximate times below are planning estimates; throughput depends on hardware, software, and cached inputs.

| Stage | Scripts, in order | Manuscript output | Approximate time |
|---|---|---|---|
| `qc` | `scripts/01_qc/run_mriqc.sh`, `patches/sitecustomize.py` | QC and S1 input | 1–4 h |
| `freesurfer` | `scripts/02_freesurfer/run_cohort.sh`, `run_subject.sh` | FreeSurfer reconstructions, eTIV | 1–3 days at several workers |
| `atlas` | `scripts/03_mind/atlas/01_prepare_schaefer200_atlas.py`, `02_validate_cohort_inputs.py`, `03_write_analysis_specification.py` | Cohort mapping, Schaefer/Yeo labels, S2/S3 inputs | Minutes |
| `mind` | `scripts/03_mind/networks/01_prepare_mind_software.py`, `02_map_annotations.py`, `03_compute_mind_networks.py`, `04_compile_network_tables.py` | 66 matrices, masked edges and counts | 30 min–3 h |
| `primary` | `analysis/04_primary/01_prepare_primary_composites.R`, `02_run_primary_models.R`, `03_write_summary_tables.R` | Primary result estimates and permutation P values | 1–10 min |
| `nbs` | `analysis/05_nbs/01_prepare_nbs_inputs.py`, `run_targeted_nbs.m`, `run_thresholds.m`, `03_summarize_nbs_results.py` | S3; component edges/nodes for S4 and annotation | Minutes–hours |
| `reference` | `analysis/06_reference/01_prepare_specificity_composites.R`, `02_run_specificity_models.R`, `03_write_specificity_summary.R` | S2, reference/direct contrasts | 1–15 min |
| `nodes` | `analysis/07_nodes/01_prepare_node_metrics.R`, `02_run_node_models.R`, `03_write_node_summary.R` | S4 | 5–45 min |
| `robustness` | `analysis/08_sensitivity/primary_robustness.R` | Primary intervals, whole cortex, eTIV and S5A | 1–10 min |
| `clinical` | `analysis/08_sensitivity/qc_clinical.R` | Table 1, S1, CJV in S5A and S7 | 1–15 min |
| `surface` | `analysis/08_sensitivity/extract_surface_morphometry.py`, `analyze_surface_morphometry.R` | S5C, plus auxiliary regional diagnostics | 5–60 min |
| `parcellation` | `analysis/08_sensitivity/compute_alternative_parcellations.py`, `analyze_parcellations.R` | S5B | 30 min–3 h |
| `maps` | `analysis/10_spatial/maps/01_build_annotation_maps.py` | MIND effect and component maps | Minutes |
| `terms` | `analysis/10_spatial/neurosynth/01_prepare_neurosynth_term_matrix.py` | Neurosynth term matrix for S8–S9 | 15 min–3 h |
| `decode` | `analysis/10_spatial/neurosynth/02_run_neurosynth_decoding.py` | S8 and degree decoding in S9; saved Váša assignments | 5–60 min |
| `spatial` | `analysis/10_spatial/annotate_components.py` | Coverage, S9 background checks, component part of S10 | 1–10 min |
| `vbm-preprocess` | `analysis/09_vbm/preprocess.py`, `segment.m` | S6 normalized, smoothed gray-matter input | 30 min–4 h |
| `vbm` | `analysis/09_vbm/analyze.py` | S6 voxelwise FWE and atlas overlap | 15 min–3 h |

This repository writes analysis tables rather than journal-formatted tables. It omits publication figure generation; the outputs above provide the numerical inputs for network, component, regional, and annotation panels. The VBM script still writes a local registration QC mosaic. See [the output index](docs/outputs.md) for exact artifact names.

Selected downstream analyses were checked against saved results using existing derivatives. An end-to-end run from the raw T1-weighted images has not been validated.

### References

When referring to this workflow, cite the accompanying study, the OpenNeuro dataset snapshot (`10.18112/openneuro.ds004144.v1.0.2`), the dataset paper by Balducci et al. (`10.1038/s41597-022-01677-9`), and the MIND methods paper by Sebenius et al. (`10.1038/s41593-023-01376-7`).
