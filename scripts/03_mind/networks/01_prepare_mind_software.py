#!/usr/bin/env python3
# Purpose: 01 prepare mind software.
# Inputs: Official MIND source and container.
# Outputs: Source and version records.
# Manuscript: MIND construction.

"""Record the fixed MIND source and runtime environment."""

import csv
import subprocess
from os import environ as cfg
from pathlib import Path

stage = Path(cfg["MIND_STAGE"])
container = cfg["FS_CONTAINER"]
commit = cfg["MIND_COMMIT"]
with (stage / "software" / "mind_software.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "source_url",
            "commit_hash",
            "local_path",
            "patch_status",
            "python_version",
            "freesurfer_executable_or_container",
            "freesurfer_version",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerow(
        {
            "source_url": "https://github.com/isebenius/MIND",
            "commit_hash": commit,
            "local_path": cfg["MIND_HOME"],
            "patch_status": "unmodified",
            "python_version": "container Python recorded in python_package_versions.tsv",
            "freesurfer_executable_or_container": container,
            "freesurfer_version": "7.3.2",
        }
    )
probe = "import sys,numpy,pandas,scipy,nibabel; print('python\\t'+sys.version.split()[0]); [print(k+'\\t'+v) for k,v in [('numpy',numpy.__version__),('pandas',pandas.__version__),('scipy',scipy.__version__),('nibabel',nibabel.__version__)]]"
result = subprocess.run(
    [cfg["APPTAINER_BIN"], "exec", container, "python3", "-c", probe],
    check=True,
    capture_output=True,
    text=True,
)
with (stage / "software" / "python_package_versions.tsv").open(
    "w", newline=""
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["package", "version"])
    writer.writerows(
        line.split("\t", 1) for line in result.stdout.splitlines() if "\t" in line
    )
