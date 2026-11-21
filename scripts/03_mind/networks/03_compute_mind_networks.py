#!/usr/bin/env python3
# Purpose: 03 compute mind networks.
# Inputs: Native annotations and vertex feature files.
# Outputs: Subject MIND matrices.
# Manuscript: MIND construction.

"""Run the official compute_MIND function for each mapped subject."""

import csv
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import environ as cfg
from pathlib import Path

import pandas as pd

stage = Path(cfg["MIND_STAGE"])
root = Path(cfg["FM_WORK"])
sys.path.insert(0, cfg["MIND_HOME"])
from MIND import compute_MIND  # noqa: E402

parcellation = cfg["SCHAEFER_200"]
features = cfg["MIND_FEATURES"].split()
mapping_path = (
    Path(cfg["ATLAS_STAGE"]) / "outputs" / "ds004144_subject_freesurfer_mapping.tsv"
)
lookup = pd.read_csv(
    Path(cfg["ATLAS_STAGE"]) / "outputs" / "schaefer200_parcel_lookup.tsv", sep="\t"
)
expected = [
    f"{hemi}_{name}" for hemi, name in zip(lookup["hemisphere"], lookup["parcel_name"])
]
with mapping_path.open(newline="") as handle:
    subjects = [
        row
        for row in csv.DictReader(handle, delimiter="\t")
        if row["mapping_status"] == "complete"
    ]
pending_subjects = [
    row
    for row in subjects
    if not (
        stage / "outputs" / "matrices" / f"{row['subject_id']}_mind_matrix.tsv"
    ).is_file()
]


def compute(row):
    subject = row["subject_id"]
    matrix = compute_MIND(
        str(stage / "work" / "mind_subjects" / subject),
        features,
        parcellation,
        filter_vertices=cfg["MIND_FILTER_VERTICES"] == "true",
        resample=cfg["MIND_RESAMPLE"] == "true",
    )
    if set(matrix.index) != set(expected):
        missing = sorted(set(expected) - set(matrix.index))
        extra = sorted(set(matrix.index) - set(expected))
        raise RuntimeError(
            f"{subject}: unexpected parcel labels; missing={missing}; extra={extra}"
        )
    matrix = matrix.reindex(index=expected, columns=expected)
    matrix.index = lookup["parcel_name"]
    matrix.columns = lookup["parcel_name"]
    matrix.to_csv(
        stage / "outputs" / "matrices" / f"{subject}_mind_matrix.tsv",
        sep="\t",
        float_format="%.17g",
    )
    return subject


failures = []
if pending_subjects:
    with ProcessPoolExecutor(
        max_workers=min(int(cfg["MIND_WORKERS"]), len(pending_subjects))
    ) as executor:
        futures = {
            executor.submit(compute, row): row["subject_id"] for row in pending_subjects
        }
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as error:
                failures.append(f"{futures[future]}: {error}")
if failures:
    raise RuntimeError("\n".join(failures))
