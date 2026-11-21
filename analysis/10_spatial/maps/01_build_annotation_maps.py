#!/usr/bin/env python3
# Purpose: Build cortical effect maps and NBS component maps for spatial annotation.
# Inputs: MIND, atlas, node, and NBS tables.
# Outputs: Parcel maps and component edges.
# Manuscript: S8–S10.

from __future__ import annotations

import hashlib
import shutil
import urllib.request
from os import environ as cfg
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import statsmodels.api as sm

stage = Path(cfg["MAP_STAGE"])
specification = Path(cfg["ATLAS_STAGE"])
networks = Path(cfg["MIND_STAGE"])
primary = Path(cfg["PRIMARY_STAGE"])
nbs = Path(cfg["NBS_STAGE"])
specificity = Path(cfg["REFERENCE_STAGE"])
nodes = Path(cfg["NODE_STAGE"])
pairs = ("SalVentAttn-Default", "SalVentAttn-DorsAttn")
volume_name = "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
volume_url = (
    "https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/"
    "stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/"
    "MNI/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hedges_g(fm: np.ndarray, hc: np.ndarray) -> float:
    pooled = np.sqrt(
        ((len(fm) - 1) * np.var(fm, ddof=1) + (len(hc) - 1) * np.var(hc, ddof=1))
        / (len(fm) + len(hc) - 2)
    )
    return float(
        (1 - 3 / (4 * (len(fm) + len(hc)) - 9)) * (np.mean(fm) - np.mean(hc)) / pooled
    )


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def component_key(frame: pd.DataFrame) -> pd.Series:
    return (
        frame["network_pair"].astype(str)
        + "|"
        + frame["component_id"].astype(float).astype(int).astype(str)
    )


def build_lookup() -> pd.DataFrame:
    source = pd.read_csv(
        specification / "outputs" / "schaefer200_parcel_lookup.tsv", sep="\t"
    )
    edge_dictionary = pd.read_csv(
        networks / "outputs" / "tables" / "mind_edge_dictionary.tsv", sep="\t"
    )
    records = (
        pd.concat(
            [
                edge_dictionary[["node_i", "parcel_i", "network_i"]].rename(
                    columns={
                        "node_i": "parcel_id",
                        "parcel_i": "parcel_name",
                        "network_i": "network",
                    }
                ),
                edge_dictionary[["node_j", "parcel_j", "network_j"]].rename(
                    columns={
                        "node_j": "parcel_id",
                        "parcel_j": "parcel_name",
                        "network_j": "network",
                    }
                ),
            ],
            ignore_index=True,
        )
        .drop_duplicates("parcel_id")
        .sort_values("parcel_id")
    )
    source = source.rename(columns={"parcel_identifier": "parcel_id"})
    lookup = records.merge(
        source[["parcel_id", "hemisphere", "parcel_name", "network"]],
        on=["parcel_id", "parcel_name", "network"],
        how="left",
        validate="one_to_one",
    )
    lookup["hemisphere"] = lookup["hemisphere"].map({"lh": "L", "rh": "R"})
    lookup = (
        lookup[["parcel_id", "parcel_name", "hemisphere", "network"]]
        .sort_values("parcel_id")
        .reset_index(drop=True)
    )
    if (
        len(lookup) != 200
        or lookup["parcel_id"].tolist() != list(range(1, 201))
        or lookup.isna().any().any()
    ):
        raise RuntimeError(
            "Parcel lookup is inconsistent with the MIND edge dictionary"
        )
    return lookup


def prepare_atlas(lookup: pd.DataFrame) -> None:
    atlas_dir = stage / "atlas"
    atlas_dir.mkdir(parents=True, exist_ok=True)
    source_atlas = specification / "atlas"
    for hemisphere in ("lh", "rh"):
        name = f"{hemisphere}.Schaefer2018_200Parcels_7Networks_order.annot"
        shutil.copy2(source_atlas / name, atlas_dir / name)
    volume = atlas_dir / volume_name
    if not volume.exists() or volume.stat().st_size == 0:
        urllib.request.urlretrieve(volume_url, volume)
    if sha256(volume) != cfg["SCHAEFER_200_MNI_SHA256"]:
        raise ValueError("Schaefer 200 MNI atlas differs from the recorded input")
    labels = np.unique(np.rint(nib.load(volume).get_fdata()).astype(int))
    labels = labels[labels > 0]
    if labels.tolist() != list(range(1, 201)):
        raise RuntimeError("MNI volume atlas labels are not 1 through 200")
    write_tsv(
        pd.DataFrame(
            [
                {
                    "file": volume.name,
                    "source": volume_url,
                    "size_bytes": volume.stat().st_size,
                    "sha256": sha256(volume),
                },
                *[
                    {
                        "file": name,
                        "source": pd.read_csv(
                            specification / "outputs" / "atlas_provenance.tsv", sep="\t"
                        )
                        .set_index("file")
                        .loc[name, "source"],
                        "size_bytes": (atlas_dir / name).stat().st_size,
                        "sha256": sha256(atlas_dir / name),
                    }
                    for name in (
                        "lh.Schaefer2018_200Parcels_7Networks_order.annot",
                        "rh.Schaefer2018_200Parcels_7Networks_order.annot",
                    )
                ],
            ]
        ),
        atlas_dir / "atlas_file_provenance.tsv",
    )
    write_tsv(
        lookup.rename(columns={"parcel_id": "id", "parcel_name": "label"})[
            ["id", "label", "hemisphere"]
        ].assign(structure="cortex"),
        atlas_dir / "schaefer200_abagen_atlas_info.tsv",
    )
    image = nib.load(volume)
    data = np.rint(image.get_fdata()).astype(int)
    centroids = []
    for row in lookup.itertuples(index=False):
        coordinates = np.argwhere(data == row.parcel_id)
        if len(coordinates) == 0:
            raise RuntimeError(f"Volume atlas lacks parcel {row.parcel_id}")
        xyz = nib.affines.apply_affine(image.affine, coordinates).mean(axis=0)
        centroids.append(
            {
                "parcel_id": row.parcel_id,
                "parcel_name": row.parcel_name,
                "hemisphere": row.hemisphere,
                "network": row.network,
                "x": xyz[0],
                "y": xyz[1],
                "z": xyz[2],
            }
        )
    write_tsv(pd.DataFrame(centroids), atlas_dir / "schaefer200_mni_centroids.tsv")


def build_node_effect_maps(lookup: pd.DataFrame) -> None:
    results = pd.read_csv(nodes / "outputs" / "node_contribution_results.tsv", sep="\t")
    output = lookup.copy()
    for pair in pairs:
        prefix = pair.replace("-", "_")
        selected = results.loc[
            results["network_pair"] == pair,
            [
                "node",
                "beta_group",
                "t_group",
                "hedges_g",
                "p_two_sided_permutation",
                "p_fdr_module",
                "in_fm_gt_hc_nbs_component",
            ],
        ].rename(
            columns={
                "node": "parcel_id",
                "beta_group": f"{prefix}_beta",
                "t_group": f"{prefix}_t",
                "hedges_g": f"{prefix}_hedges_g",
                "p_two_sided_permutation": f"{prefix}_p",
                "p_fdr_module": f"{prefix}_q",
                "in_fm_gt_hc_nbs_component": f"{prefix}_in_nbs_component",
            }
        )
        output = output.merge(
            selected, on="parcel_id", how="left", validate="one_to_one"
        )
    write_tsv(output, stage / "maps" / "target_node_effect_maps.tsv")


def build_component_maps(lookup: pd.DataFrame) -> None:
    summaries = pd.read_csv(nbs / "outputs" / "nbs_results_summary.tsv", sep="\t")
    component_nodes = pd.read_csv(nbs / "outputs" / "nbs_component_nodes.tsv", sep="\t")
    component_edges = pd.read_csv(nbs / "outputs" / "nbs_component_edges.tsv", sep="\t")
    significant = summaries.loc[
        summaries["network_pair"].isin(pairs)
        & summaries["direction"].eq("FM_gt_HC")
        & summaries["p_fwe_across_masks_directions"].lt(0.05)
        & summaries["component_id"].notna(),
        ["network_pair", "component_id"],
    ].copy()
    significant["component_key"] = component_key(significant)
    component_nodes = component_nodes.copy()
    component_edges = component_edges.copy()
    component_nodes["component_key"] = component_key(component_nodes)
    component_edges["component_key"] = component_key(component_edges)
    selected_nodes = component_nodes.loc[
        component_nodes["component_key"].isin(significant["component_key"])
    ].copy()
    selected_edges = component_edges.loc[
        component_edges["component_key"].isin(significant["component_key"])
    ].copy()
    selected_nodes = selected_nodes.loc[
        selected_nodes["network_pair"].isin(pairs)
        & selected_nodes["direction"].eq("FM_gt_HC")
    ].drop(columns="component_key")
    selected_edges = selected_edges.loc[
        selected_edges["network_pair"].isin(pairs)
        & selected_edges["direction"].eq("FM_gt_HC")
    ].drop(columns="component_key")
    write_tsv(selected_edges, stage / "components" / "nbs_fm_gt_hc_component_edges.tsv")
    write_tsv(selected_nodes, stage / "components" / "nbs_fm_gt_hc_component_nodes.tsv")
    output = lookup.copy()
    for pair in pairs:
        prefix = pair.replace("-", "_")
        pair_nodes = selected_nodes.loc[
            selected_nodes["network_pair"] == pair, ["node", "degree_in_component"]
        ].rename(
            columns={
                "node": "parcel_id",
                "degree_in_component": f"{prefix}_component_degree",
            }
        )
        pair_edges = selected_edges.loc[selected_edges["network_pair"] == pair]
        strength = (
            pd.concat(
                [
                    pair_edges[["node_i", "t_value"]].rename(
                        columns={"node_i": "parcel_id"}
                    ),
                    pair_edges[["node_j", "t_value"]].rename(
                        columns={"node_j": "parcel_id"}
                    ),
                ],
                ignore_index=True,
            )
            .groupby("parcel_id", as_index=False)["t_value"]
            .sum()
            .rename(columns={"t_value": f"{prefix}_component_strength"})
        )
        pair_map = pair_nodes.merge(
            strength, on="parcel_id", how="outer", validate="one_to_one"
        )
        output = output.merge(
            pair_map, on="parcel_id", how="left", validate="one_to_one"
        )
        output[f"{prefix}_component_degree"] = (
            output[f"{prefix}_component_degree"].fillna(0).astype(int)
        )
        output[f"{prefix}_component_strength"] = output[
            f"{prefix}_component_strength"
        ].fillna(0.0)
    write_tsv(output, stage / "maps" / "nbs_component_node_maps.tsv")
    return significant


def build_full_cortex_map(lookup: pd.DataFrame) -> None:
    wide = pd.read_csv(
        networks / "outputs" / "tables" / "mind_edges_wide.tsv", sep="\t"
    )
    edges = pd.read_csv(
        networks / "outputs" / "tables" / "mind_edge_dictionary.tsv", sep="\t"
    ).sort_values("edge_id")
    edge_columns = [f"edge_{edge_id:05d}" for edge_id in edges["edge_id"]]
    if len(edge_columns) != 19900 or not set(edge_columns).issubset(wide.columns):
        raise RuntimeError("MIND edge table is inconsistent with the edge dictionary")
    edge_values = wide[edge_columns].to_numpy(dtype=float)
    strength = np.zeros((len(wide), 200), dtype=float)
    for index, edge in enumerate(edges.itertuples(index=False)):
        strength[:, edge.node_i - 1] += edge_values[:, index]
        strength[:, edge.node_j - 1] += edge_values[:, index]
    strength /= 199.0
    group = np.where(wide["group"].eq("FM"), 1.0, 0.0)
    age = wide["age"].to_numpy(dtype=float)
    design = sm.add_constant(np.column_stack((group, age)), has_constant="add")
    rows = []
    for parcel_index, row in enumerate(lookup.itertuples(index=False)):
        metric = strength[:, parcel_index]
        fit = sm.OLS(metric, design).fit()
        fm = metric[group == 1]
        hc = metric[group == 0]
        rows.append(
            {
                "parcel_id": row.parcel_id,
                "parcel_name": row.parcel_name,
                "hemisphere": row.hemisphere,
                "network": row.network,
                "beta_group": fit.params[1],
                "standardized_beta_group": fit.params[1]
                * np.std(group, ddof=1)
                / np.std(metric, ddof=1),
                "hedges_g": hedges_g(fm, hc),
                "mean_FM": np.mean(fm),
                "mean_HC": np.mean(hc),
            }
        )
    write_tsv(
        pd.DataFrame(rows),
        stage / "maps" / "full_cortex_node_strength_effects.tsv",
    )


def copy_tables() -> None:
    (stage / "tables").mkdir(parents=True, exist_ok=True)
    for source, destination in (
        (
            primary / "outputs" / "primary_composite_results.tsv",
            "primary_composite_results.tsv",
        ),
        (
            specificity / "outputs" / "network_specificity_results.tsv",
            "network_specificity_results.tsv",
        ),
        (
            nodes / "outputs" / "node_contribution_results.tsv",
            "node_contribution_results.tsv",
        ),
        (
            networks / "outputs" / "masks" / "network_pair_edge_counts.tsv",
            "network_pair_edge_counts.tsv",
        ),
    ):
        shutil.copy2(source, stage / "tables" / destination)


def write_manifest() -> None:
    entries = [
        (
            "target_node_effect_maps",
            "maps/target_node_effect_maps.tsv",
            [
                "SalVentAttn_Default_beta",
                "SalVentAttn_Default_t",
                "SalVentAttn_Default_hedges_g",
                "SalVentAttn_Default_p",
                "SalVentAttn_Default_q",
                "SalVentAttn_Default_in_nbs_component",
                "SalVentAttn_DorsAttn_beta",
                "SalVentAttn_DorsAttn_t",
                "SalVentAttn_DorsAttn_hedges_g",
                "SalVentAttn_DorsAttn_p",
                "SalVentAttn_DorsAttn_q",
                "SalVentAttn_DorsAttn_in_nbs_component",
            ],
            "target-pair parcel annotation",
        ),
        (
            "nbs_component_node_maps",
            "maps/nbs_component_node_maps.tsv",
            [
                "SalVentAttn_Default_component_degree",
                "SalVentAttn_Default_component_strength",
                "SalVentAttn_DorsAttn_component_degree",
                "SalVentAttn_DorsAttn_component_strength",
            ],
            "component topology annotation",
        ),
        (
            "full_cortex_node_strength_effects",
            "maps/full_cortex_node_strength_effects.tsv",
            ["beta_group", "standardized_beta_group", "hedges_g", "mean_FM", "mean_HC"],
            "whole-cortex spatial annotation",
        ),
    ]
    rows = [
        {
            "map_name": name,
            "file": file,
            "n_parcels": 200,
            "value_column": value,
            "recommended_use": use,
        }
        for name, file, values, use in entries
        for value in values
    ]
    write_tsv(pd.DataFrame(rows), stage / "maps" / "map_manifest.tsv")


def write_checks(lookup: pd.DataFrame, significant: pd.DataFrame) -> None:
    atlas_dir = stage / "atlas"
    checks = [
        ("parcel_lookup_200_rows", len(lookup) == 200),
        ("parcel_ids_unique", lookup["parcel_id"].is_unique),
        (
            "parcel_order_matches_mind_outputs",
            lookup["parcel_id"].tolist() == list(range(1, 201)),
        ),
        (
            "surface_annotation_files_present",
            all(
                (
                    atlas_dir
                    / f"{hemisphere}.Schaefer2018_200Parcels_7Networks_order.annot"
                ).is_file()
                for hemisphere in ("lh", "rh")
            ),
        ),
        ("mni_volume_atlas_present", (atlas_dir / volume_name).is_file()),
        (
            "abagen_atlas_info_200_rows",
            len(pd.read_csv(atlas_dir / "schaefer200_abagen_atlas_info.tsv", sep="\t"))
            == 200,
        ),
        (
            "centroid_table_200_rows",
            len(pd.read_csv(atlas_dir / "schaefer200_mni_centroids.tsv", sep="\t"))
            == 200,
        ),
        (
            "target_node_effect_map_200_rows",
            len(
                pd.read_csv(
                    stage / "maps" / "target_node_effect_maps.tsv", sep="\t"
                )
            )
            == 200,
        ),
        (
            "nbs_component_node_map_200_rows",
            len(
                pd.read_csv(
                    stage / "maps" / "nbs_component_node_maps.tsv", sep="\t"
                )
            )
            == 200,
        ),
        (
            "full_cortex_node_strength_map_200_rows",
            len(
                pd.read_csv(
                    stage / "maps" / "full_cortex_node_strength_effects.tsv",
                    sep="\t",
                )
            )
            == 200,
        ),
        (
            "required_primary_pairs_present",
            set(pairs).issubset(
                set(
                    pd.read_csv(
                        stage / "tables" / "primary_composite_results.tsv", sep="\t"
                    )["network_pair"]
                )
            ),
        ),
        (
            "significant_fm_gt_hc_components_present",
            set(pairs).issubset(set(significant["network_pair"])),
        ),
    ]
    write_tsv(
        pd.DataFrame(checks, columns=["check", "status"]).assign(
            status=lambda frame: np.where(frame["status"], "passed", "failed")
        ),
        stage / "checks" / "annotation_map_checks.tsv",
    )


def main() -> None:
    lookup = build_lookup()
    write_tsv(lookup, stage / "atlas" / "schaefer200_parcel_lookup.tsv")
    prepare_atlas(lookup)
    build_node_effect_maps(lookup)
    significant = build_component_maps(lookup)
    build_full_cortex_map(lookup)
    copy_tables()
    write_manifest()
    write_checks(lookup, significant)


if __name__ == "__main__":
    main()
