# Manuscript output index

`work/` refers to the configured `FM_WORK`. All artifacts below are generated locally.

| Manuscript item | File names (relative to the relevant stage in `work/`) |
|---|---|
| Table 1 | `qc_clinical/participant_characteristics.tsv` |
| Primary group differences and intervals | `004_schaefer200_primary_analysis/outputs/primary_composite_results.tsv`, `primary_robustness/model_results.tsv` |
| S1, T1w IQMs | `qc_clinical/mriqc_descriptive.tsv` |
| S2, edge counts and comparators | `002_ds004144_schaefer200_mind_networks/outputs/masks/network_pair_edge_counts.tsv`, `006_schaefer200_network_specificity/outputs/network_specificity_results.tsv` |
| S3, NBS | `nbs_thresholds/nbs_threshold_results.tsv`, `spatial_annotation/component_coverage.tsv`; component edges and nodes in `005_schaefer200_targeted_nbs/outputs/` |
| S4, parcel analyses | `007_schaefer200_node_contribution/outputs/node_contribution_results.tsv` |
| S5A | `primary_robustness/model_results.tsv`, `etiv_summary.tsv`, `low_control_influence.tsv`; `qc_clinical/mriqc_sensitivity.tsv` |
| S5B | `parcellation_sensitivity/parcellation_effects.tsv` |
| S5C | `surface_morphometry/mind_morphometry_sensitivity.tsv` |
| S6 | `voxel_morphometry/vbm_summary.tsv`, `significant_voxel_overlap.tsv` |
| S7 | `qc_clinical/clinical_associations.tsv` |
| S8 | `009_schaefer200_neurosynth_annotation/outputs/neurosynth_term_decoding_results.tsv`, full-cortex effect and standardized coefficient maps |
| S9 | `spatial_annotation/functional_background.tsv`, plus component degree rows from the S8 input table |
| S10 | `spatial_annotation/receptor_background.tsv`, component analyses |

The figure inputs are the corresponding primary, reference, component, parcel, and spatial tables, including the map outputs under `008_schaefer200_annotation_maps/maps/`.
