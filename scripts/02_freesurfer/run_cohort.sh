#!/usr/bin/env bash
# Purpose: Run T1w FreeSurfer reconstructions for the BIDS sample.
# Inputs: BIDS participant table, T1w images, and optional BIDS participant IDs.
# Outputs: Participant reconstructions and shell logs.
# Manuscript: MRI Acquisition and Processing.
set -euo pipefail
if (( $# > 0 )); then
    subjects=("$@")
else
    mapfile -t subjects < <(awk -F '\t' 'NR > 1 {sub(/\r$/, "", $1); print $1}' "$BIDS_DIR/phenotype/participants.tsv")
fi
(( ${#subjects[@]} > 0 ))
jobs=()
for subject in "${subjects[@]}"; do
    if [[ -f $FS_DIR/$subject/scripts/recon-all.done ]]; then continue; fi
    bash "$FM_REPO/scripts/02_freesurfer/run_subject.sh" "$subject" &
    jobs+=("$!")
    if (( ${#jobs[@]} >= FS_MAX_JOBS )); then
        for job in "${jobs[@]}"; do wait "$job"; done
        jobs=()
    fi
done
for job in "${jobs[@]}"; do wait "$job"; done
