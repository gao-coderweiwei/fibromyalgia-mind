#!/usr/bin/env python
# Purpose: 01 prepare neurosynth term matrix.
# Inputs: Neurosynth v7 database, term panel, and atlas.
# Outputs: Term-by-parcel matrix.
# Manuscript: S8–S9.

from __future__ import annotations

import hashlib
import os
import sys
from os import environ as cfg
from pathlib import Path

import nibabel as nib
import nilearn
import nimare
import numpy as np
import pandas as pd
from nilearn.maskers import NiftiLabelsMasker
from nimare.decode.continuous import CorrelationDecoder
from nimare.io import convert_neurosynth_to_dataset
from nimare.meta.cbma.mkda import MKDAChi2
from scipy import sparse

stage = Path(cfg["NEUROSYNTH_STAGE"])
analysis_root = Path(cfg["FM_WORK"])
resources = Path(cfg["NEUROSYNTH_DATA"])
map_stage = Path(cfg["MAP_STAGE"])
atlas_img = (
    map_stage / "atlas" / "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
)
term_panel_file = Path(cfg["FM_REPO"]) / "config" / "term_panel.tsv"
term_matrix_file = (
    stage / "outputs" / "term_matrix" / "neurosynth_schaefer200_term_matrix.tsv"
)
seed = int(cfg["PERM_SEED"])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def resource_inventory() -> None:
    files = [
        resources / "data-neurosynth_version-7_coordinates.tsv.gz",
        resources / "data-neurosynth_version-7_metadata.tsv.gz",
        resources
        / "data-neurosynth_version-7_vocab-terms_source-abstract_type-tfidf_features.npz",
        resources / "data-neurosynth_version-7_vocab-terms_vocabulary.txt",
        atlas_img,
        term_panel_file,
    ]
    rows = []
    for path in files:
        rows.append(
            {
                "file": str(path),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else np.nan,
                "sha256": sha256(path) if path.exists() else "",
            }
        )
    write_tsv(pd.DataFrame(rows), stage / "outputs" / "resource_inventory.tsv")


def write_term_map_versions() -> None:
    rows = [
        {"component": "python", "version": sys.version.split()[0]},
        {"component": "nimare", "version": nimare.__version__},
        {"component": "nilearn", "version": nilearn.__version__},
        {"component": "nibabel", "version": nib.__version__},
        {"component": "numpy", "version": np.__version__},
        {"component": "pandas", "version": pd.__version__},
    ]
    write_tsv(pd.DataFrame(rows), stage / "software" / "term_map_software_versions.tsv")


def load_term_panel() -> pd.DataFrame:
    vocabulary_file = resources / "data-neurosynth_version-7_vocab-terms_vocabulary.txt"
    features_file = (
        resources
        / "data-neurosynth_version-7_vocab-terms_source-abstract_type-tfidf_features.npz"
    )
    vocabulary = vocabulary_file.read_text(encoding="utf-8").splitlines()
    counts = np.asarray((sparse.load_npz(features_file) > 0).sum(axis=0)).ravel()
    vocab = pd.DataFrame({"term": vocabulary, "n_studies": counts.astype(int)})
    panel = pd.read_csv(term_panel_file, sep="\t").drop_duplicates(["category", "term"])
    panel = panel.merge(vocab, on="term", how="left")
    panel["available"] = panel["n_studies"].notna()
    write_tsv(panel, stage / "outputs" / "neurosynth_term_panel.tsv")
    selected = panel.loc[panel["available"]].copy()
    selected["feature"] = "terms__" + selected["term"]
    selected["n_studies"] = selected["n_studies"].astype(int)
    if selected.empty:
        raise RuntimeError("No requested Neurosynth terms were available")
    return selected


def make_dataset():
    annotations = {
        "features": str(
            resources
            / "data-neurosynth_version-7_vocab-terms_source-abstract_type-tfidf_features.npz"
        ),
        "vocabulary": str(
            resources / "data-neurosynth_version-7_vocab-terms_vocabulary.txt"
        ),
    }
    return convert_neurosynth_to_dataset(
        coordinates_file=str(
            resources / "data-neurosynth_version-7_coordinates.tsv.gz"
        ),
        metadata_file=str(resources / "data-neurosynth_version-7_metadata.tsv.gz"),
        annotations_files=annotations,
        feature_groups=["terms"],
        target="mni152_2mm",
    )


def build_term_matrix(selected: pd.DataFrame) -> None:
    labels = np.unique(np.rint(nib.load(str(atlas_img)).get_fdata()).astype(int))
    labels = labels[labels > 0]
    if labels.tolist() != list(range(1, 201)):
        raise RuntimeError("Schaefer200 atlas labels are not consecutive")

    dataset = make_dataset()
    n_cores = min(16, max(1, (os.cpu_count() or 2) - 2))
    decoder = CorrelationDecoder(
        features=selected["feature"].tolist(),
        frequency_threshold=float(cfg["NEUROSYNTH_FEATURE_THRESHOLD"]),
        meta_estimator=MKDAChi2(),
        target_image="z_desc-association",
        n_cores=n_cores,
    )
    decoder.fit(dataset)

    masker = NiftiLabelsMasker(
        labels_img=str(atlas_img), strategy="mean", standardize=False
    )
    masker.fit()
    parcel_columns = [f"parcel_{index:03d}" for index in range(1, 201)]
    rows = []
    for record in selected.itertuples(index=False):
        image = decoder.results_.masker.inverse_transform(
            decoder.results_.maps[record.feature]
        )
        values = masker.transform(image).ravel()
        if len(values) != 200 or not np.isfinite(values).all():
            raise RuntimeError(f"Invalid parcellated term map for {record.term}")
        row = {
            "category": record.category,
            "term": record.term,
            "feature": record.feature,
            "n_studies": record.n_studies,
        }
        row.update(dict(zip(parcel_columns, values)))
        rows.append(row)
    write_tsv(pd.DataFrame(rows), term_matrix_file)


def main() -> None:
    np.random.seed(seed)
    resource_inventory()
    write_term_map_versions()
    selected = load_term_panel()
    build_term_matrix(selected)


if __name__ == "__main__":
    main()
