#!/usr/bin/env python3
# Purpose: 02 validate cohort inputs.
# Inputs: BIDS phenotype and FreeSurfer outputs.
# Outputs: Subject mapping table.
# Manuscript: Participants.

"""Map participant records to FreeSurfer subjects and required surfaces."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from os import environ as cfg
from pathlib import Path

import pandas as pd

SURFACES = (
    "white",
    "pial",
    "sphere.reg",
    "thickness",
    "curv",
    "volume",
    "sulc",
    "area",
)


def inspect(row: dict[str, str], subjects_dir: Path) -> dict[str, str]:
    subject_id = row["subject_id"]
    subject_dir = subjects_dir / subject_id
    result = dict(row)
    result["freesurfer_subject_id"] = subject_id
    result["freesurfer_subject_directory"] = str(subject_dir)
    for hemisphere in ("lh", "rh"):
        for surface in SURFACES:
            result[f"{hemisphere}_{surface.replace('.', '_')}_exists"] = str(
                (subject_dir / "surf" / f"{hemisphere}.{surface}").is_file()
            ).lower()
    available = all(
        value == "true" for key, value in result.items() if key.endswith("_exists")
    )
    result["mapping_status"] = "complete" if available else "missing_required_file"
    return result


def cohort_rows(
    cohort: str, participants: Path, subjects_dir: Path
) -> list[dict[str, str]]:
    frame = pd.read_csv(participants, sep="\t", dtype=str, keep_default_na=False)
    required = {"participant_id", "age", "sex", "group"}
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"{cohort}: missing columns {sorted(missing)}")
    rows = []
    for record in frame.loc[:, ["participant_id", "group", "sex", "age"]].to_dict(
        "records"
    ):
        age = pd.to_numeric(record["age"], errors="coerce")
        if pd.isna(age):
            raise RuntimeError(
                f"{cohort}: unparseable age for {record['participant_id']}"
            )
        rows.append(
            {
                "cohort": cohort,
                "subject_id": record["participant_id"],
                "group": record["group"],
                "sex": record["sex"],
                "age": f"{float(age):g}",
            }
        )
    with ThreadPoolExecutor(max_workers=min(40, max(1, len(rows)))) as executor:
        return list(executor.map(lambda row: inspect(row, subjects_dir), rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-dir", type=Path, required=True)
    parser.add_argument("--ds004144", type=Path, required=True)
    args = parser.parse_args()
    stage = args.stage_dir.resolve()
    outputs = stage / "outputs"
    outputs.mkdir(parents=True, exist_ok=True)
    definitions = (
        (
            "ds004144",
            args.ds004144 / "phenotype" / "participants.tsv",
            Path(cfg["FS_DIR"]),
        ),
    )
    base_columns = [
        "cohort",
        "subject_id",
        "group",
        "sex",
        "age",
        "freesurfer_subject_id",
        "freesurfer_subject_directory",
    ]
    surface_columns = [
        f"{hemi}_{surface.replace('.', '_')}_exists"
        for hemi in ("lh", "rh")
        for surface in SURFACES
    ]
    fieldnames = base_columns + surface_columns + ["mapping_status"]
    for cohort, participants, subjects_dir in definitions:
        if not subjects_dir.is_dir():
            raise RuntimeError(f"{cohort}: FreeSurfer subjects directory not found")
        rows = cohort_rows(cohort, participants, subjects_dir)
        with (outputs / f"{cohort}_subject_freesurfer_mapping.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
