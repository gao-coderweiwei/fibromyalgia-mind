#!/usr/bin/env python3
# Purpose: 02 map annotations.
# Inputs: FreeSurfer surfaces and fsaverage6 annotations.
# Outputs: Native annotations and work links.
# Manuscript: MIND construction.

"""Map fsaverage6 annotations to stage-local subject work directories."""

import csv
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from os import environ as cfg
from pathlib import Path

stage = Path(cfg["MIND_STAGE"])
root = Path(cfg["FM_WORK"])
mapping = (
    Path(cfg["ATLAS_STAGE"]) / "outputs" / "ds004144_subject_freesurfer_mapping.tsv"
)
atlas = Path(cfg["ATLAS_STAGE"]) / "atlas"
derivatives = Path(cfg["FS_DIR"])
container = cfg["FS_CONTAINER"]
license_file = cfg["FS_LICENSE"]
parcellation = cfg["SCHAEFER_200"]

with mapping.open(newline="") as handle:
    subjects = [
        row
        for row in csv.DictReader(handle, delimiter="\t")
        if row["mapping_status"] == "complete"
    ]


def link(path, target):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink() or path.exists():
        path.unlink()
    path.symlink_to(target)


work_subjects = stage / "work" / "subjects"
work_subjects.mkdir(parents=True, exist_ok=True)
link(work_subjects / "fsaverage6", cfg["FS_AVG6"])
for row in subjects:
    subject = row["freesurfer_subject_id"]
    link(work_subjects / subject, f"{cfg['FS_MOUNT_SUBJECTS']}/{subject}")
    out_label = stage / "outputs" / "annotations" / row["subject_id"] / "label"
    out_label.mkdir(parents=True, exist_ok=True)
    subject_work = stage / "work" / "mind_subjects" / row["subject_id"]
    for hemi in ("lh", "rh"):
        for feature in ("thickness", "curv", "volume", "sulc", "area"):
            link(
                subject_work / "surf" / f"{hemi}.{feature}",
                f"{cfg['FS_MOUNT_SUBJECTS']}/{subject}/surf/{hemi}.{feature}",
            )
        link(
            subject_work / "label" / f"{hemi}.{parcellation}.annot",
            Path("../../../../outputs")
            / "annotations"
            / row["subject_id"]
            / "label"
            / f"{hemi}.{parcellation}.annot",
        )


def map_subject(row):
    subject = row["subject_id"]
    for hemi in ("lh", "rh"):
        command = [
            cfg["APPTAINER_BIN"],
            "exec",
            "--bind",
            f"{derivatives}:{cfg['FS_MOUNT_SUBJECTS']}:ro,{atlas}:{cfg['FS_MOUNT_ATLAS']}:ro,{stage}:{cfg['FS_MOUNT_STAGE']},{license_file}:{cfg['FS_MOUNT_LICENSE']}:ro",
            container,
            "bash",
            "-lc",
            f"export SUBJECTS_DIR={cfg['FS_MOUNT_STAGE']}/work/subjects FS_LICENSE={cfg['FS_MOUNT_LICENSE']}; mri_surf2surf --srcsubject fsaverage6 --trgsubject {row['freesurfer_subject_id']} --hemi {hemi} --sval-annot {cfg['FS_MOUNT_ATLAS']}/{hemi}.{parcellation}.annot --tval {cfg['FS_MOUNT_STAGE']}/outputs/annotations/{subject}/label/{hemi}.{parcellation}.annot",
        ]
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    return subject


failures = []
with ThreadPoolExecutor(
    max_workers=min(int(cfg["MIND_WORKERS"]), len(subjects))
) as executor:
    futures = {executor.submit(map_subject, row): row["subject_id"] for row in subjects}
    for future in as_completed(futures):
        try:
            future.result()
        except Exception as error:
            failures.append(f"{futures[future]}: {error}")
if failures:
    raise RuntimeError("\n".join(failures))
