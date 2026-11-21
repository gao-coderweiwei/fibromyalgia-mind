#!/usr/bin/env bash
# Purpose: Run the selected reproduction stage using the repository configuration.
# Inputs: Repository configuration, public inputs, and upstream stage outputs.
# Outputs: Analysis artifacts in the configured work directory.
# Manuscript: Methods, Results, and Supplementary Methods.
set -euo pipefail
FM_REPO=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export FM_REPO
cd -- "$FM_REPO"
set -a
source config/defaults.sh
local_config=${FM_CONFIG:-config/local.sh}
if [[ -f $local_config ]]; then source "$local_config"; fi
for key in BIDS_DIR FM_WORK FS_DIR MRIQC_OUT FS_LOG_DIR TEMPLATEFLOW_HOME FS_CONTAINER MRIQC_CONTAINER FS_LICENSE MIND_HOME NBS_HOME SPM_HOME NEUROSYNTH_DATA PYTHON_BIN NEUROSYNTH_PYTHON SPIN_PYTHON MATLAB_BIN RSCRIPT_BIN APPTAINER_BIN; do
    value=${!key}
    case "$key" in
        *_BIN|NEUROSYNTH_PYTHON|SPIN_PYTHON) [[ $value == */* ]] || continue ;;
    esac
    if [[ $value != /* ]]; then printf -v "$key" '%s/%s' "$FM_REPO" "$value"; fi
done
for key in ATLAS MIND PRIMARY NBS REFERENCE NODE MAP NEUROSYNTH ROBUST THRESHOLD CLINICAL SURFACE PARCELLATION SPATIAL VBM; do
    name=${key}_STAGE_NAME
    printf -v "${key}_STAGE" '%s/%s' "$FM_WORK" "${!name}"
done
set +a
export PYTHONDONTWRITEBYTECODE=1
command_name=${1:-help}
shift $(( $# > 0 ? 1 : 0 ))
if [[ $command_name == help || $command_name == --help ]]; then
    printf '%s\n' 'Usage: bash scripts/run.sh STAGE' 'Stages: qc freesurfer atlas mind primary nbs reference nodes robustness clinical surface parcellation maps terms decode spatial vbm-preprocess vbm'
    exit 0
fi

init_stage() { mkdir -p "$1"/{outputs,config,software,logs}; }
run_models() {
    local stage=$1 folder=$2 prepare=$3 model=$4 summarize=$5 max_jobs=$6
    shift 6
    init_stage "$stage"
    if [[ -n $prepare ]]; then "$RSCRIPT_BIN" --vanilla "$folder/$prepare" "$stage"; fi
    local -a jobs=()
    local metric job
    for metric; do
        if [[ $metric == *'|'* ]]; then
            "$RSCRIPT_BIN" --vanilla "$folder/$model" "$stage" "${metric%%|*}" "${metric##*|}" &
        else
            "$RSCRIPT_BIN" --vanilla "$folder/$model" "$stage" "$metric" &
        fi
        jobs+=("$!")
        if (( ${#jobs[@]} >= max_jobs )); then
            for job in "${jobs[@]}"; do wait "$job"; done
            jobs=()
        fi
    done
    for job in "${jobs[@]}"; do wait "$job"; done
    "$RSCRIPT_BIN" --vanilla "$folder/$summarize" "$stage"
}

case "$command_name" in
    qc) bash scripts/01_qc/run_mriqc.sh ;;
    freesurfer) bash scripts/02_freesurfer/run_cohort.sh "$@" ;;
    atlas)
        init_stage "$ATLAS_STAGE"
        for script in 01_prepare_schaefer200_atlas.py 02_validate_cohort_inputs.py 03_write_analysis_specification.py; do
            args=(--stage-dir "$ATLAS_STAGE")
            if [[ $script == 02_validate_cohort_inputs.py ]]; then args+=(--ds004144 "$BIDS_DIR"); fi
            "$PYTHON_BIN" "scripts/03_mind/atlas/$script" "${args[@]}"
        done
        (cd "$ATLAS_STAGE/atlas" && sha256sum -c "$FM_REPO/config/schaefer200.sha256")
        ;;
    mind)
        init_stage "$MIND_STAGE"
        mkdir -p "$MIND_STAGE/outputs/"{matrices,tables,masks,checks}
        if [[ $(git -C "$MIND_HOME" rev-parse HEAD) != "$MIND_COMMIT" ]] || ! git -C "$MIND_HOME" diff --quiet HEAD --; then
            printf '%s\n' 'MIND must be checked out at the configured commit with clean tracked files.' >&2; exit 1
        fi
        "$PYTHON_BIN" scripts/03_mind/networks/01_prepare_mind_software.py
        "$PYTHON_BIN" scripts/03_mind/networks/02_map_annotations.py
        for script in 03_compute_mind_networks.py 04_compile_network_tables.py; do
            "$APPTAINER_BIN" exec --bind "$FM_WORK:$MIND_MOUNT_WORK,$FS_DIR:$FS_MOUNT_SUBJECTS:ro,$FM_REPO/scripts/03_mind/networks:$MIND_MOUNT_CODE:ro,$MIND_HOME:$MIND_MOUNT_DEP:ro,$FS_LICENSE:$FS_MOUNT_LICENSE:ro" \
                "$FS_CONTAINER" env "FM_WORK=$MIND_MOUNT_WORK" "MIND_STAGE=$MIND_MOUNT_WORK/$MIND_STAGE_NAME" \
                "MIND_HOME=$MIND_MOUNT_DEP" "FS_LICENSE=$FS_MOUNT_LICENSE" python3 "$MIND_MOUNT_CODE/$script"
        done
        ;;
    primary)
        run_models "$PRIMARY_STAGE" analysis/04_primary 01_prepare_primary_composites.R 02_run_primary_models.R 03_write_summary_tables.R 2 SalVentAttn_Default SalVentAttn_DorsAttn
        ;;
    nbs)
        init_stage "$NBS_STAGE"; init_stage "$THRESHOLD_STAGE"
        "$PYTHON_BIN" analysis/05_nbs/01_prepare_nbs_inputs.py
        jobs=()
        for pair in SalVentAttn-Default SalVentAttn-DorsAttn; do
            for direction in FM_gt_HC HC_gt_FM; do
                "$MATLAB_BIN" -batch "addpath('analysis/05_nbs'); run_targeted_nbs(getenv('NBS_STAGE'), '$pair', '$direction');" >"$NBS_STAGE/logs/${pair}_${direction}.log" 2>&1 &
                jobs+=("$!")
            done
        done
        for job in "${jobs[@]}"; do wait "$job"; done
        "$MATLAB_BIN" -singleCompThread -batch "addpath('analysis/05_nbs'); run_thresholds;"
        "$PYTHON_BIN" analysis/05_nbs/03_summarize_nbs_results.py
        ;;
    reference)
        run_models "$REFERENCE_STAGE" analysis/06_reference 01_prepare_specificity_composites.R 02_run_specificity_models.R 03_write_specificity_summary.R 10 \
            within_SalVentAttn within_Default within_DorsAttn Default_DorsAttn SalVentAttn_SomMot Default_SomMot DorsAttn_SomMot \
            SalVentAttn_Default_minus_Default_DorsAttn SalVentAttn_DorsAttn_minus_Default_DorsAttn SalVentAttn_DorsAttn_minus_SalVentAttn_Default
        ;;
    nodes)
        init_stage "$NODE_STAGE"
        folder=analysis/07_nodes
        "$RSCRIPT_BIN" --vanilla "$folder/01_prepare_node_metrics.R" "$NODE_STAGE"
        mapfile -t metrics < <(awk -F '\t' 'NR > 1 {print $7 "|" $8}' "$NODE_STAGE/outputs/node_contribution_values.tsv" | sort -u)
        run_models "$NODE_STAGE" "$folder" '' 02_run_node_models.R 03_write_node_summary.R "$NODE_WORKERS" "${metrics[@]}"
        ;;
    robustness) init_stage "$ROBUST_STAGE"; "$RSCRIPT_BIN" --vanilla analysis/08_sensitivity/primary_robustness.R ;;
    clinical) init_stage "$CLINICAL_STAGE"; "$RSCRIPT_BIN" --vanilla analysis/08_sensitivity/qc_clinical.R ;;
    surface)
        init_stage "$SURFACE_STAGE"
        "$PYTHON_BIN" analysis/08_sensitivity/extract_surface_morphometry.py
        "$RSCRIPT_BIN" --vanilla analysis/08_sensitivity/analyze_surface_morphometry.R
        ;;
    parcellation)
        init_stage "$PARCELLATION_STAGE"
        "$PYTHON_BIN" analysis/08_sensitivity/compute_alternative_parcellations.py
        "$RSCRIPT_BIN" --vanilla analysis/08_sensitivity/analyze_parcellations.R
        ;;
    maps) init_stage "$MAP_STAGE"; "$PYTHON_BIN" analysis/10_spatial/maps/01_build_annotation_maps.py ;;
    terms) init_stage "$NEUROSYNTH_STAGE"; "$NEUROSYNTH_PYTHON" analysis/10_spatial/neurosynth/01_prepare_neurosynth_term_matrix.py ;;
    decode) init_stage "$NEUROSYNTH_STAGE"; "$SPIN_PYTHON" analysis/10_spatial/neurosynth/02_run_neurosynth_decoding.py ;;
    spatial) init_stage "$SPATIAL_STAGE"; "$PYTHON_BIN" analysis/10_spatial/annotate_components.py ;;
    vbm-preprocess) init_stage "$VBM_STAGE"; "$PYTHON_BIN" analysis/09_vbm/preprocess.py "$VBM_STAGE/images" ;;
    vbm) init_stage "$VBM_STAGE"; "$PYTHON_BIN" analysis/09_vbm/analyze.py "$VBM_STAGE/images" ;;
    *) printf 'Unknown stage: %s\n' "$command_name" >&2; exit 2 ;;
esac
