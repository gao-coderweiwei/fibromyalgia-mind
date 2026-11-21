#!/usr/bin/env python
# Purpose: 02 run neurosynth decoding.
# Inputs: Term matrix, atlas, and study maps.
# Outputs: Decoding statistics and spin indices.
# Manuscript: S8–S9.

from __future__ import annotations

import platform
import sys
from os import environ as cfg
from pathlib import Path

import neuromaps
import nibabel as nib
import numpy as np
import pandas as pd
import scipy
from neuromaps import nulls
from neuromaps.images import annot_to_gifti
from statsmodels.stats.multitest import multipletests

stage = Path(cfg["NEUROSYNTH_STAGE"])
analysis_root = Path(cfg["FM_WORK"])
map_stage = Path(cfg["MAP_STAGE"])
term_matrix_file = (
    stage / "outputs" / "term_matrix" / "neurosynth_schaefer200_term_matrix.tsv"
)
study_maps_file = Path(cfg["FM_REPO"]) / "config" / "study_maps.tsv"
seed = int(cfg["PERM_SEED"])
n_perm = int(cfg["N_SPINS"])


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def zscore_vector(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    centered = values - values.mean()
    norm = np.linalg.norm(centered)
    if norm == 0 or not np.isfinite(norm):
        raise RuntimeError("A map has zero or invalid variance")
    return centered / norm


def zscore_rows(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(centered, axis=1, keepdims=True)
    if np.any(norm == 0) or not np.isfinite(norm).all():
        raise RuntimeError("A term map has zero or invalid variance")
    return centered / norm


def zscore_columns(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=0, keepdims=True)
    norm = np.linalg.norm(centered, axis=0, keepdims=True)
    if np.any(norm == 0) or not np.isfinite(norm).all():
        raise RuntimeError("A null map has zero or invalid variance")
    return centered / norm


def load_study_maps() -> pd.DataFrame:
    lookup = pd.read_csv(map_stage / "atlas" / "schaefer200_parcel_lookup.tsv", sep="\t")
    definitions = pd.read_csv(study_maps_file, sep="\t")
    rows = []
    for record in definitions.itertuples(index=False):
        source = analysis_root / record.source_file
        frame = pd.read_csv(source, sep="\t")
        merged = lookup.merge(
            frame[["parcel_id", record.value_column]],
            on="parcel_id",
            how="left",
            validate="one_to_one",
        )
        values = merged[record.value_column].to_numpy(dtype=float)
        if len(values) != 200 or not np.isfinite(values).all():
            raise RuntimeError(f"Invalid study map: {record.map_name}")
        row = {"map_name": record.map_name, "role": record.role}
        row.update(
            {
                f"parcel_{index:03d}": value
                for index, value in enumerate(values, start=1)
            }
        )
        rows.append(row)
    output = pd.DataFrame(rows)
    write_tsv(output, stage / "outputs" / "prepared_study_maps.tsv")
    return output


def load_term_matrix() -> tuple[pd.DataFrame, np.ndarray]:
    terms = pd.read_csv(term_matrix_file, sep="\t")
    parcel_columns = [f"parcel_{index:03d}" for index in range(1, 201)]
    values = terms[parcel_columns].to_numpy(dtype=float)
    if len(terms) == 0 or values.shape[1] != 200 or not np.isfinite(values).all():
        raise RuntimeError("Invalid Neurosynth term matrix")
    return terms[["category", "term", "feature", "n_studies"]].copy(), values


def make_spin_indices() -> np.ndarray:
    atlas_dir = map_stage / "atlas"
    parcellation = annot_to_gifti(
        (
            str(atlas_dir / "lh.Schaefer2018_200Parcels_7Networks_order.annot"),
            str(atlas_dir / "rh.Schaefer2018_200Parcels_7Networks_order.annot"),
        )
    )
    spins = nulls.vasa(
        None,
        atlas="fsaverage",
        density="41k",
        parcellation=parcellation,
        n_perm=n_perm,
        seed=seed,
    ).astype(int)
    if spins.shape != (200, n_perm):
        raise RuntimeError("Spin index matrix has an unexpected shape")
    np.save(stage / "outputs" / f"spin_indices_vasa_seed{cfg['PERM_SEED']}.npy", spins)
    return spins


def decode(
    study_maps: pd.DataFrame,
    term_info: pd.DataFrame,
    term_values: np.ndarray,
    spins: np.ndarray,
) -> pd.DataFrame:
    term_z = zscore_rows(term_values)
    parcel_columns = [f"parcel_{index:03d}" for index in range(1, 201)]
    rows = []
    for record in study_maps.itertuples(index=False):
        values = np.asarray(
            [getattr(record, column) for column in parcel_columns], dtype=float
        )
        observed = term_z @ zscore_vector(values)
        null_values = values[spins]
        null_corr = term_z @ zscore_columns(null_values)
        p_values = (
            np.sum(np.abs(null_corr) >= np.abs(observed[:, None]), axis=1) + 1
        ) / (n_perm + 1)
        _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
        result = term_info.copy()
        result.insert(0, "map_name", record.map_name)
        result.insert(1, "map_role", record.role)
        result["pearson_r"] = observed
        result["p_spin"] = p_values
        result["p_fdr_bh"] = q_values
        result["direction"] = np.where(result["pearson_r"] >= 0, "positive", "negative")
        result["rank_abs_r"] = (
            result["pearson_r"].abs().rank(method="first", ascending=False).astype(int)
        )
        rows.append(result)
    return pd.concat(rows, ignore_index=True).sort_values(["map_name", "rank_abs_r"])


def summarize(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (map_name, category), group in results.groupby(
        ["map_name", "category"], sort=True
    ):
        positive = group.sort_values("pearson_r", ascending=False).iloc[0]
        negative = group.sort_values("pearson_r", ascending=True).iloc[0]
        rows.append(
            {
                "map_name": map_name,
                "category": category,
                "n_terms": len(group),
                "n_fdr_0_05": int((group["p_fdr_bh"] < 0.05).sum()),
                "top_positive_term": positive["term"],
                "top_positive_r": positive["pearson_r"],
                "top_positive_q": positive["p_fdr_bh"],
                "top_negative_term": negative["term"],
                "top_negative_r": negative["pearson_r"],
                "top_negative_q": negative["p_fdr_bh"],
            }
        )
    return pd.DataFrame(rows)


def write_versions() -> None:
    rows = [
        {"component": "python", "version": sys.version.split()[0]},
        {"component": "platform", "version": platform.platform()},
        {"component": "numpy", "version": np.__version__},
        {"component": "pandas", "version": pd.__version__},
        {"component": "scipy", "version": scipy.__version__},
        {"component": "nibabel", "version": nib.__version__},
        {"component": "neuromaps", "version": neuromaps.__version__},
    ]
    write_tsv(pd.DataFrame(rows), stage / "software" / "software_versions.tsv")


def write_checks(
    study_maps: pd.DataFrame,
    term_info: pd.DataFrame,
    term_values: np.ndarray,
    results: pd.DataFrame,
    spins: np.ndarray,
) -> None:
    expected_rows = len(study_maps) * len(term_info)
    checks = [
        ("term_matrix_exists", term_matrix_file.is_file()),
        ("term_matrix_terms_available", len(term_info) >= 60),
        ("term_matrix_200_parcels", term_values.shape[1] == 200),
        ("study_maps_available", len(study_maps) == 4),
        ("study_maps_200_parcels", study_maps.filter(like="parcel_").shape[1] == 200),
        ("spin_indices_shape", spins.shape == (200, n_perm)),
        ("decoding_results_complete", len(results) == expected_rows),
        (
            "decoding_values_finite",
            np.isfinite(
                results[["pearson_r", "p_spin", "p_fdr_bh"]].to_numpy(dtype=float)
            ).all(),
        ),
        ("fdr_range_valid", results["p_fdr_bh"].between(0, 1).all()),
    ]
    write_tsv(
        pd.DataFrame(checks, columns=["check", "passed"]),
        stage / "outputs" / "checks" / "neurosynth_decoding_checks.tsv",
    )


def main() -> None:
    np.random.seed(seed)
    study_maps = load_study_maps()
    term_info, term_values = load_term_matrix()
    spins = make_spin_indices()
    results = decode(study_maps, term_info, term_values, spins)
    write_tsv(results, stage / "outputs" / "neurosynth_term_decoding_results.tsv")
    write_tsv(
        summarize(results), stage / "outputs" / "neurosynth_term_category_summary.tsv"
    )
    top_terms = (
        results.groupby("map_name", group_keys=False)
        .apply(lambda frame: frame.nsmallest(20, "rank_abs_r"))
        .reset_index(drop=True)
    )
    write_tsv(top_terms, stage / "outputs" / "neurosynth_top_terms_by_map.tsv")
    write_checks(study_maps, term_info, term_values, results, spins)
    write_versions()


if __name__ == "__main__":
    main()
