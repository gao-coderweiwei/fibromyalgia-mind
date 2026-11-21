#!/usr/bin/env bash
# Purpose: Run T1w MRIQC with the QI2 compatibility patch.
# Inputs: BIDS dataset, MRIQC container, TemplateFlow, and compatibility patch.
# Outputs: Participant reports and T1w IQM table in the configured output directory.
# Manuscript: MRI Acquisition and Processing; Table S1.
set -euo pipefail
mkdir -p "$MRIQC_OUT" "$FM_WORK/mriqc-work" "$TEMPLATEFLOW_HOME"
APPTAINERENV_PYTHONPATH="$MRIQC_MOUNT_PATCH" APPTAINERENV_TEMPLATEFLOW_HOME="$MRIQC_MOUNT_TEMPLATE" \
    "$APPTAINER_BIN" run --cleanenv \
        --bind "$BIDS_DIR:$MRIQC_MOUNT_BIDS:ro" \
        --bind "$MRIQC_OUT:$MRIQC_MOUNT_OUT" \
        --bind "$FM_WORK/mriqc-work:$MRIQC_MOUNT_WORK" \
        --bind "$TEMPLATEFLOW_HOME:$MRIQC_MOUNT_TEMPLATE" \
        --bind "$FM_REPO/scripts/01_qc/patches:$MRIQC_MOUNT_PATCH:ro" \
        "$MRIQC_CONTAINER" "$MRIQC_MOUNT_BIDS" "$MRIQC_MOUNT_OUT" participant group -m T1w \
        --nprocs "$MRIQC_NPROCS" --omp-nthreads "$MRIQC_OMP" --mem "$MRIQC_MEM" --no-sub --notrack -w "$MRIQC_MOUNT_WORK"
