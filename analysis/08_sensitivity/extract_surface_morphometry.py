#!/usr/bin/env python3
# Purpose: Extract surface morphometry.
# Inputs: FreeSurfer surfaces and annotations.
# Outputs: Parcel morphometry table.
# Manuscript: S5C.

import csv
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from os import environ as cfg
from pathlib import Path
from tempfile import TemporaryDirectory

STAGE = Path(cfg["SURFACE_STAGE"])
PROJECT = Path(cfg["FM_WORK"])
SPEC = Path(cfg["ATLAS_STAGE"]) / "outputs"
ANNOTATIONS = Path(cfg["MIND_STAGE"]) / "outputs" / "annotations"
SUBJECTS = Path(cfg["FS_DIR"])
PARCELLATION = cfg["SCHAEFER_200"]
CONTAINER = cfg["FS_CONTAINER"]
LICENSE = cfg["FS_LICENSE"]


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract(task, temporary):
    subject, freesurfer_id, hemisphere, parcels = task
    output = f"{subject}_{hemisphere}.stats"
    annotation = f'{cfg["FS_MOUNT_ANNOTATIONS"]}/{subject}/label/{hemisphere}.{PARCELLATION}.annot'
    command = [
        cfg["APPTAINER_BIN"],
        "exec",
        "--bind",
        f'{SUBJECTS}:{cfg["FS_MOUNT_SUBJECTS"]}:ro,{ANNOTATIONS}:{cfg["FS_MOUNT_ANNOTATIONS"]}:ro,{temporary}:{cfg["FS_MOUNT_STATS"]},{LICENSE}:{cfg["FS_MOUNT_LICENSE"]}:ro',
        CONTAINER,
        "mris_anatomical_stats",
        "-th3",
        "-mgz",
        "-cortex",
        f'{cfg["FS_MOUNT_SUBJECTS"]}/{freesurfer_id}/label/{hemisphere}.cortex.label',
        "-f",
        f'{cfg["FS_MOUNT_STATS"]}/{output}',
        "-b",
        "-a",
        annotation,
        freesurfer_id,
        hemisphere,
        "white",
    ]
    env = os.environ | {
        "APPTAINERENV_SUBJECTS_DIR": cfg["FS_MOUNT_SUBJECTS"],
        "APPTAINERENV_FS_LICENSE": cfg["FS_MOUNT_LICENSE"],
        "APPTAINERENV_OMP_NUM_THREADS": "1",
    }
    result = subprocess.run(command, env=env, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            f"{subject} {hemisphere}: {result.stderr.strip() or result.stdout.strip()}"
        )

    by_parcel = {}
    with (temporary / output).open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parcel, vertices, area, volume, thickness, *_ = line.split()
            if parcel in parcels:
                by_parcel[parcel] = (
                    int(vertices),
                    float(area),
                    float(volume),
                    float(thickness),
                )
    if set(by_parcel) != set(parcels) or any(row[0] <= 0 for row in by_parcel.values()):
        raise ValueError(f"{subject} {hemisphere}: incomplete parcellation")
    return [
        (subject, hemisphere, parcel, network, *by_parcel[parcel])
        for parcel, network in parcels.items()
    ]


def main():
    lookup = read_tsv(SPEC / "schaefer200_parcel_lookup.tsv")
    parcels = {
        hemisphere: {
            row["parcel_name"]: row["network"]
            for row in lookup
            if row["hemisphere"] == hemisphere
        }
        for hemisphere in ("lh", "rh")
    }
    subjects = read_tsv(SPEC / "ds004144_subject_freesurfer_mapping.tsv")
    tasks = [
        (row["subject_id"], row["freesurfer_subject_id"], hemisphere, labels)
        for row in subjects
        if row["mapping_status"] == "complete"
        for hemisphere, labels in parcels.items()
    ]
    with TemporaryDirectory(prefix="fm_schaefer_stats_") as directory:
        with ThreadPoolExecutor(max_workers=8) as pool:
            batches = list(pool.map(lambda task: extract(task, Path(directory)), tasks))
    with (STAGE / "parcel_measures.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "subject_id",
                "hemisphere",
                "parcel",
                "network",
                "n_vertices",
                "area_mm2",
                "volume_mm3",
                "thickness_mm",
            ]
        )
        writer.writerows(row for batch in batches for row in batch)


if __name__ == "__main__":
    main()
