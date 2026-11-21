#!/usr/bin/env python3
# Purpose: 01 prepare schaefer200 atlas.
# Inputs: CBIG atlas resources.
# Outputs: Parcel lookup, checks, and provenance.
# Manuscript: Cortical parcellation; S2.

"""Download and validate the Schaefer 200/Yeo7 atlas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import tempfile
import urllib.request
from collections import Counter
from os import environ as cfg
from pathlib import Path

import nibabel as nib

FILES = (
    "lh.Schaefer2018_200Parcels_7Networks_order.annot",
    "rh.Schaefer2018_200Parcels_7Networks_order.annot",
    "Schaefer2018_200Parcels_7Networks_order.txt",
)
BASE_URL = (
    "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/"
    "stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/"
    "Parcellations/FreeSurfer5.3/fsaverage6/label"
)
TEXT_URL = (
    "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/"
    "stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/"
    "Parcellations/MNI/freeview_lut/Schaefer2018_200Parcels_7Networks_order.txt"
)
NETWORKS = ("Visual", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default")
NETWORK_TOKENS = {
    "Visual": "Vis",
    "SomMot": "SomMot",
    "DorsAttn": "DorsAttn",
    "SalVentAttn": "SalVentAttn",
    "Limbic": "Limbic",
    "Cont": "Cont",
    "Default": "Default",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(path: Path) -> str:
    url = TEXT_URL if path.name.endswith(".txt") else f"{BASE_URL}/{path.name}"
    if path.exists():
        return url
    with urllib.request.urlopen(url, timeout=60) as response, path.open("wb") as output:
        shutil.copyfileobj(response, output)
    return url


def parse_labels(path: Path, hemisphere: str) -> list[dict[str, str]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.strip().split()
            if len(fields) < 2 or not fields[0].isdigit():
                continue
            name = fields[1]
            if f"_{hemisphere.upper()}_" in name:
                rows.append(
                    {
                        "hemisphere": hemisphere,
                        "parcel_identifier": fields[0],
                        "parcel_name": name,
                    }
                )
    return rows


def validate_with_freesurfer(atlas: Path, work: Path) -> tuple[str, str]:
    container = Path(cfg["FS_CONTAINER"])
    license_file = Path(cfg["FS_LICENSE"])
    if not container.exists() or not license_file.exists():
        return "unavailable", "FreeSurfer container or license unavailable"
    with tempfile.TemporaryDirectory(dir=work) as temporary:
        outputs = Path(temporary)
        binds = [
            f"{atlas / FILES[0]}:{cfg['FS_AVG_LABEL']}/{FILES[0]}:ro",
            f"{atlas / FILES[1]}:{cfg['FS_AVG_LABEL']}/{FILES[1]}:ro",
            f"{outputs}:{outputs}",
            f"{license_file}:{cfg['FS_NATIVE_LICENSE']}:ro",
        ]
        commands = []
        for hemisphere in ("lh", "rh"):
            output = outputs / hemisphere
            output.mkdir(exist_ok=True)
            commands.append(
                [
                    cfg["APPTAINER_BIN"],
                    "exec",
                    "--cleanenv",
                    *sum((["--bind", item] for item in binds), []),
                    str(container),
                    "env",
                    f"FS_LICENSE={cfg['FS_NATIVE_LICENSE']}",
                    "mri_annotation2label",
                    "--subject",
                    "fsaverage",
                    "--hemi",
                    hemisphere,
                    "--annotation",
                    cfg["SCHAEFER_200"],
                    "--outdir",
                    str(output),
                ]
            )
        try:
            for command in commands:
                subprocess.run(command, check=True, capture_output=True, text=True)
        except (OSError, subprocess.CalledProcessError) as error:
            detail = (
                error.stderr.strip()
                if isinstance(error, subprocess.CalledProcessError)
                else str(error)
            )
            return "failed", detail[-500:]
    return "passed", "mri_annotation2label"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-dir", type=Path, required=True)
    args = parser.parse_args()
    stage = args.stage_dir.resolve()
    atlas = stage / "atlas"
    outputs = stage / "outputs"
    atlas.mkdir(parents=True, exist_ok=True)
    outputs.mkdir(parents=True, exist_ok=True)

    provenance = []
    for filename in FILES:
        path = atlas / filename
        source = download(path)
        provenance.append(
            {
                "file": filename,
                "source": source,
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "status": "available",
            }
        )

    labels = parse_labels(atlas / FILES[2], "lh") + parse_labels(atlas / FILES[2], "rh")
    if len(labels) != 200:
        raise RuntimeError(f"Expected 200 parcel labels, found {len(labels)}")
    names = [row["parcel_name"] for row in labels]
    if len(set(names)) != 200 or any(not name for name in names):
        raise RuntimeError("Parcel identifiers or names are invalid")
    for row in labels:
        row["network"] = next(
            (
                network
                for network, token in NETWORK_TOKENS.items()
                if f"_{token}_" in row["parcel_name"]
            ),
            "",
        )
    if any(not row["network"] for row in labels):
        raise RuntimeError("Missing Yeo7 network label")

    for hemisphere, filename in (("lh", FILES[0]), ("rh", FILES[1])):
        _, ctab, annot_names = nib.freesurfer.read_annot(
            str(atlas / filename), orig_ids=True
        )
        cortical_names = [
            name
            for name in annot_names
            if not name.decode("utf-8").startswith("Background+")
        ]
        if len(cortical_names) != 100 or ctab.shape[0] != 101:
            raise RuntimeError(f"{hemisphere} annotation does not contain 100 parcels")

    fs_status, fs_detail = validate_with_freesurfer(atlas, outputs)
    with (outputs / "atlas_provenance.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["file", "source", "size_bytes", "sha256", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(provenance)
    with (outputs / "schaefer200_parcel_lookup.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["hemisphere", "parcel_identifier", "parcel_name", "network"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(labels)
    counts = Counter(row["network"] for row in labels)
    with (outputs / "yeo7_network_parcel_counts.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["network", "parcel_count"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(
            {"network": network, "parcel_count": counts[network]}
            for network in NETWORKS
        )
    (outputs / "atlas_validation.tsv").write_text(
        "check\tvalue\n"
        f"freesurfer_annotation_readability\t{fs_status}\n"
        f"freesurfer_validation_detail\t{fs_detail}\n",
        encoding="utf-8",
    )
    if fs_status != "passed":
        raise RuntimeError(f"FreeSurfer annotation validation failed: {fs_detail}")


if __name__ == "__main__":
    main()
