#!/usr/bin/env python3
# Purpose: 03 write analysis specification.
# Inputs: Atlas and cohort validation tables.
# Outputs: Analysis specification and checks.
# Manuscript: Methods; S2–S3.

"""Write the fixed analysis specification and minimal quality checks."""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import subprocess
import sys
from os import environ as cfg
from pathlib import Path

import pandas as pd

SPECIFICATION = (
    ("atlas", cfg["SCHAEFER_200"]),
    ("network_labels", "Yeo7"),
    ("mind_features", "CT,MC,Vol,SD,SA"),
    ("model", "metric ~ group + age"),
    ("group_direction", "FM - HC"),
    ("primary_network_pair_1", "SalVentAttn-Default"),
    ("primary_network_pair_2", "SalVentAttn-DorsAttn"),
    ("within_network_reference", "SalVentAttn,Default,DorsAttn"),
    (
        "cross_network_reference",
        "Default-DorsAttn,SalVentAttn-SomMot,Default-SomMot,DorsAttn-SomMot",
    ),
    ("primary_permutations", "50000"),
    ("nbs_threshold", "t > 2.0"),
    ("nbs_component_statistic", "mass"),
)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def command_output(command: list[str]) -> str:
    try:
        return (
            subprocess.check_output(command, text=True, stderr=subprocess.STDOUT)
            .strip()
            .splitlines()[0]
        )
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-dir", type=Path, required=True)
    args = parser.parse_args()
    stage = args.stage_dir.resolve()
    outputs, config, software = stage / "outputs", stage / "config", stage / "software"
    config.mkdir(parents=True, exist_ok=True)
    software.mkdir(parents=True, exist_ok=True)
    with (config / "analysis_specification.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["parameter", "value"])
        writer.writerows(SPECIFICATION)

    atlas = read_tsv(outputs / "atlas_provenance.tsv")
    lookup = read_tsv(outputs / "schaefer200_parcel_lookup.tsv")
    network_counts = read_tsv(outputs / "yeo7_network_parcel_counts.tsv")
    atlas_validation = read_tsv(outputs / "atlas_validation.tsv")
    mappings = {
        "ds004144": read_tsv(outputs / "ds004144_subject_freesurfer_mapping.tsv")
    }
    expected = {"ds004144": {"FM": 33, "HC": 33}}
    required_networks = {
        "Visual",
        "SomMot",
        "DorsAttn",
        "SalVentAttn",
        "Limbic",
        "Cont",
        "Default",
    }
    surface_columns = [
        column for column in mappings["ds004144"].columns if column.endswith("_exists")
    ]
    atlas_status = atlas["status"].eq("available").all()
    readable = (
        atlas_validation.loc[
            atlas_validation["check"] == "freesurfer_annotation_readability", "value"
        ]
        .eq("passed")
        .all()
    )
    mapping_complete = {
        cohort: frame["mapping_status"].eq("complete").all()
        for cohort, frame in mappings.items()
    }
    group_counts = {
        cohort: frame["group"].value_counts().to_dict() == expected[cohort]
        for cohort, frame in mappings.items()
    }
    surfaces_available = all(
        frame[surface_columns].eq("true").all().all() for frame in mappings.values()
    )
    checks = (
        ("official_atlas_files_available", atlas_status),
        ("atlas_readable", readable),
        (
            "atlas_200_parcels_confirmed",
            len(lookup) == 200 and lookup["parcel_identifier"].nunique() == 200,
        ),
        (
            "yeo7_labels_confirmed",
            set(network_counts["network"]) == required_networks
            and network_counts["parcel_count"].astype(int).sum() == 200,
        ),
        ("ds004144_subject_mapping_complete", mapping_complete["ds004144"]),
        ("ds004144_group_counts_confirmed", group_counts["ds004144"]),
        ("required_surface_files_available", surfaces_available),
        (
            "analysis_specification_written",
            (config / "analysis_specification.tsv").is_file(),
        ),
    )
    with (outputs / "specification_checks.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "status"])
        writer.writerows(
            (name, "passed" if value else "failed") for name, value in checks
        )
    container = cfg["FS_CONTAINER"]
    fs_version = command_output(
        [cfg["APPTAINER_BIN"], "exec", "--cleanenv", container, "recon-all", "-version"]
    )
    versions = (
        ("python", sys.version.split()[0]),
        ("pandas", importlib.metadata.version("pandas")),
        ("numpy", importlib.metadata.version("numpy")),
        ("freesurfer", fs_version),
        ("container_or_executable_path", container),
    )
    with (software / "software_versions.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["software", "version_or_path"])
        writer.writerows(versions)


if __name__ == "__main__":
    main()
