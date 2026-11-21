# Purpose: Preprocess.
# Inputs: BIDS T1w and SPM.
# Outputs: Modulated and smoothed grey matter maps.
# Manuscript: S6.

import argparse
import csv
import gzip
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from os import environ as cfg
from pathlib import Path

STAGE = Path(cfg["VBM_STAGE"])
PROJECT = Path(cfg["FM_WORK"])
DATASET = Path(cfg["BIDS_DIR"])
PRIMARY = Path(cfg["PRIMARY_STAGE"]) / "outputs/primary_composite_values.tsv"
SPM = Path(cfg["SPM_HOME"])


def prepare(subject, work):
    folder = work / subject
    folder.mkdir(exist_ok=True)
    target = folder / "t1.nii"
    if not (folder / "gm_s8.nii").is_file():
        with gzip.open(
            DATASET / subject / "anat" / f"{subject}_T1w.nii.gz", "rb"
        ) as source:
            with target.open("wb") as handle:
                shutil.copyfileobj(source, handle)
    return target


def run_batch(manifest):
    command = (
        f"addpath('{Path(__file__).resolve().parent}'); segment('{manifest}', '{SPM}');"
    )
    result = subprocess.run(
        [cfg["MATLAB_BIN"], "-singleCompThread", "-batch", command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode:
        raise RuntimeError(result.stdout[-3500:])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("work", type=Path)
    parser.add_argument("--workers", type=int, default=int(cfg["VBM_WORKERS"]))
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    with PRIMARY.open(newline="") as handle:
        inputs = [
            prepare(row["subject_id"], work)
            for row in csv.DictReader(handle, delimiter="\t")
        ]
    manifests = []
    for index in range(args.workers):
        manifest = work / f"batch_{index:02}.txt"
        manifest.write_text("\n".join(map(str, inputs[index :: args.workers])) + "\n")
        manifests.append(manifest)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(run_batch, manifests))


if __name__ == "__main__":
    main()
