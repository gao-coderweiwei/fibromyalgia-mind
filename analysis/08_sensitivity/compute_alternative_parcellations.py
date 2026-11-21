# Purpose: Compute alternative parcellations.
# Inputs: FreeSurfer surfaces and 100/300 atlases.
# Outputs: Alternative network means.
# Manuscript: S5B.

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from os import environ as cfg
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlretrieve

import numpy as np

STAGE = Path(cfg["PARCELLATION_STAGE"])
PROJECT = Path(cfg["FM_WORK"])
SPEC = Path(cfg["ATLAS_STAGE"])
NETWORKS = Path(cfg["MIND_STAGE"])
PRIMARY = Path(cfg["PRIMARY_STAGE"]) / "outputs"
DERIVATIVES = Path(cfg["FS_DIR"])
CONTAINER = cfg["FS_CONTAINER"]
LICENSE = cfg["FS_LICENSE"]
CBIG_COMMIT = cfg["CBIG_COMMIT"]
ATLAS_URL = (
    f"https://raw.githubusercontent.com/ThomasYeoLab/CBIG/{CBIG_COMMIT}/"
    "stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/"
    "Parcellations/FreeSurfer5.3/fsaverage6/label"
)
FEATURES = cfg["MIND_FEATURES"].split()
TARGETS = ("Default", "DorsAttn")

sys.dont_write_bytecode = True
sys.path.insert(0, cfg["MIND_HOME"])
from get_vertex_df import get_vertex_df  # noqa: E402
from MIND_helpers import get_KDTree, get_KL  # noqa: E402


def parcellation(resolution):
    return f"Schaefer2018_{resolution}Parcels_7Networks_order"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def download_atlases():
    atlas = STAGE / "atlas"
    atlas.mkdir(exist_ok=True)
    for resolution in (100, 200, 300):
        for hemisphere in ("lh", "rh"):
            name = f"{hemisphere}.{parcellation(resolution)}.annot"
            path = atlas / name
            if not path.exists():
                with TemporaryDirectory(dir=atlas) as directory:
                    temporary = Path(directory) / name
                    urlretrieve(f"{ATLAS_URL}/{name}", temporary)
                    if resolution == 200:
                        if (
                            temporary.read_bytes()
                            != (SPEC / "atlas" / name).read_bytes()
                        ):
                            raise ValueError(
                                f"CBIG reference differs from original: {name}"
                            )
                    else:
                        temporary.replace(path)
    return atlas


def map_annotations(subject, freesurfer_id, resolution, atlas, labels):
    env = os.environ | {
        "APPTAINERENV_FS_LICENSE": cfg["FS_MOUNT_LICENSE"],
        "APPTAINERENV_OMP_NUM_THREADS": "1",
    }
    mount = (
        f'{DERIVATIVES}:{cfg["FS_MOUNT_SUBJECTS"]}:ro,{atlas}:{cfg["FS_MOUNT_ATLAS"]}:ro,{labels}:{cfg["FS_MOUNT_LABELS"]},'
        f'{LICENSE}:{cfg["FS_MOUNT_LICENSE"]}:ro'
    )
    command = [
        cfg["APPTAINER_BIN"],
        "exec",
        "--bind",
        mount,
        CONTAINER,
        "bash",
        "-c",
        f'set -eu; SUBJECTS_DIR=$(mktemp -d {cfg["FS_MAP_TMP"]}); export SUBJECTS_DIR; '
        "trap 'rm -rf \"$SUBJECTS_DIR\"' EXIT; "
        f'ln -s {cfg["FS_AVG6"]} "$SUBJECTS_DIR/fsaverage6"; '
        f'ln -s "{cfg["FS_MOUNT_SUBJECTS"]}/$1" "$SUBJECTS_DIR/$1"; '
        "for hemi in lh rh; do "
        'mri_surf2surf --srcsubject fsaverage6 --trgsubject "$1" --hemi "$hemi" '
        f'--sval-annot "{cfg["FS_MOUNT_ATLAS"]}/$hemi.$2.annot" --tval "{cfg["FS_MOUNT_LABELS"]}/$hemi.$2.annot"; done',
        "map_annotations",
        freesurfer_id,
        parcellation(resolution),
    ]
    result = subprocess.run(command, env=env, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"{subject}, {resolution} parcels: {result.stderr or result.stdout}"
        )


def target_composites(work, resolution):
    with Path(os.devnull).open("w") as sink, redirect_stdout(sink):
        data, regions, columns = get_vertex_df(
            str(work), FEATURES, parcellation(resolution)
        )
    if len(regions) != resolution:
        raise ValueError(f"{work}: found {len(regions)} out of {resolution} parcels")
    for column in columns:
        data[column] = (data[column] - data[column].mean()) / data[column].std()
    cortex = data.loc[data["Label"].isin(regions)]
    if len(cortex[columns].drop_duplicates()) / len(cortex) < 0.8:
        raise ValueError(f"{work}: too many repeated feature vectors")
    groups = {name: [] for name in ("SalVentAttn", *TARGETS)}
    for parcel, vertices in cortex.groupby("Label"):
        network = next((name for name in groups if f"_{name}_" in parcel), None)
        if network:
            points = vertices[columns].to_numpy()
            groups[network].append((points, get_KDTree(points)))
    results = []
    for target in TARGETS:
        edges = [
            1
            / (
                1
                + get_KL(left, right, left_tree, right_tree)
                + get_KL(right, left, right_tree, left_tree)
            )
            for left, left_tree in groups["SalVentAttn"]
            for right, right_tree in groups[target]
        ]
        if not edges or not np.isfinite(edges).all():
            raise ValueError(f"{work}: invalid edges for {target}")
        results.append((target, len(edges), np.mean(edges)))
    return results


def process(task):
    subject, freesurfer_id, resolution, atlas = task
    with TemporaryDirectory(prefix=f"mind_{resolution}_{subject}_") as directory:
        work = Path(directory)
        (work / "surf").symlink_to(DERIVATIVES / freesurfer_id / "surf")
        if resolution == 200:
            (work / "label").symlink_to(
                NETWORKS / "outputs/annotations" / subject / "label"
            )
        else:
            (work / "label").mkdir()
            map_annotations(subject, freesurfer_id, resolution, atlas, work / "label")
        return [
            (subject, resolution, f"SalVentAttn_{target}", count, value)
            for target, count, value in target_composites(work, resolution)
        ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=int(cfg["PARCELLATION_WORKERS"]))
    args = parser.parse_args()
    atlas = download_atlases()
    cohort = read_tsv(PRIMARY / "primary_composite_values.tsv")
    mapping = {
        row["subject_id"]: row["freesurfer_subject_id"]
        for row in read_tsv(SPEC / "outputs/ds004144_subject_freesurfer_mapping.tsv")
    }
    for group in ("FM", "HC"):
        baseline = next(row for row in cohort if row["group"] == group)
        subject = baseline["subject_id"]
        for _, _, endpoint, _, value in process(
            (subject, mapping[subject], 200, atlas)
        ):
            if not np.isclose(value, float(baseline[endpoint]), rtol=0, atol=1e-10):
                raise ValueError(
                    f"Schaefer200 reference mismatch: {subject}, {endpoint}"
                )

    tasks = [
        (row["subject_id"], mapping[row["subject_id"]], resolution, atlas)
        for resolution in map(int, cfg["ALT_RESOLUTIONS"].split())
        for row in cohort
    ]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        batches = list(pool.map(process, tasks))
    with (STAGE / "alternative_composites.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["subject_id", "parcels", "endpoint", "n_edges", "mind"])
        writer.writerows(row for batch in batches for row in batch)


if __name__ == "__main__":
    main()
