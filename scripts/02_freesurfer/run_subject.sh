#!/usr/bin/env bash
# Purpose: Reconstruct a single participant's T1w anatomy with FreeSurfer.
# Inputs: BIDS participant ID, T1w image, container, and local license.
# Outputs: FreeSurfer subject tree and execution log.
# Manuscript: MRI Acquisition and Processing.
set -euo pipefail
subject=${1:?Supply a BIDS participant ID}
[[ $subject =~ ^sub-[a-zA-Z0-9]+$ ]]
[[ -s $BIDS_DIR/$subject/anat/${subject}_T1w.nii.gz && -s $FS_LICENSE && -f $FS_CONTAINER ]]
mkdir -p "$FS_DIR" "$FS_LOG_DIR"
"$APPTAINER_BIN" exec --cleanenv \
    --bind "$BIDS_DIR:$FS_MOUNT_BIDS:ro" --bind "$FS_DIR:$FS_MOUNT_OUT" \
    --bind "$FS_LICENSE:$FS_NATIVE_LICENSE:ro" "$FS_CONTAINER" \
    recon-all -s "$subject" -sd "$FS_MOUNT_OUT" -i "$FS_MOUNT_BIDS/$subject/anat/${subject}_T1w.nii.gz" \
    -all -3T -parallel -openmp "$FS_OPENMP" -time >"$FS_LOG_DIR/${subject}_recon-all.log" 2>&1
