#!/usr/bin/env python3
# Purpose: 01 prepare nbs inputs.
# Inputs: MIND edges, masks, and age.
# Outputs: Network-pair design and MAT inputs.
# Manuscript: S3.

from __future__ import annotations

import csv
from os import environ as cfg
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import savemat

ROOT = Path(cfg["FM_WORK"])
STAGE = Path(cfg["NBS_STAGE"])
NETWORK = Path(cfg["MIND_STAGE"])
PAIRS = ("SalVentAttn-Default", "SalVentAttn-DorsAttn")


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def main() -> None:
    wide_path = NETWORK / "outputs" / "tables" / "mind_edges_wide.tsv"
    if not wide_path.exists():
        wide_path = wide_path.with_suffix(".tsv.gz")
    wide = pd.read_csv(wide_path, sep="\t")
    edge_table = read_tsv(NETWORK / "outputs" / "tables" / "mind_edge_dictionary.tsv")
    edge_table["edge_id"] = pd.to_numeric(edge_table["edge_id"], errors="raise").astype(
        int
    )
    for column in ("node_i", "node_j"):
        edge_table[column] = pd.to_numeric(edge_table[column], errors="raise").astype(
            int
        )
    masks = read_tsv(NETWORK / "outputs" / "masks" / "analysis_network_masks.tsv")
    primary_results = (
        Path(cfg["PRIMARY_STAGE"]) / "outputs" / "primary_composite_results.tsv"
    )
    if not primary_results.is_file():
        raise RuntimeError("Primary composite results are unavailable")
    if len(wide) != 66 or set(wide["group"]) != {"FM", "HC"}:
        raise RuntimeError("Cohort records are invalid")
    ages = pd.to_numeric(wide["age"], errors="raise").to_numpy(float)
    group = wide["group"].map({"HC": 0, "FM": 1})
    if group.isna().any() or not np.isfinite(ages).all():
        raise RuntimeError("Group or age values are invalid")
    design = np.column_stack((np.ones(len(wide)), group.to_numpy(float), ages))
    output_rows = []
    for pair in PAIRS:
        selected = masks.loc[masks["mask_name"] == pair]
        if len(selected) != 1:
            raise RuntimeError(f"Missing target mask: {pair}")
        edge_ids = np.array(
            [int(value) for value in selected.iloc[0]["edge_ids"].split(",")], dtype=int
        )
        edges = edge_table.loc[edge_table["edge_id"].isin(edge_ids)].sort_values(
            "edge_id"
        )
        if len(edges) != len(edge_ids):
            raise RuntimeError(f"Mask edge IDs are absent from the edge table: {pair}")
        nodes = np.sort(
            np.unique(
                np.concatenate((edges["node_i"].to_numpy(), edges["node_j"].to_numpy()))
            )
        )
        node_index = {node: index for index, node in enumerate(nodes)}
        matrices = np.zeros((len(nodes), len(nodes), len(wide)), dtype=np.float64)
        for edge in edges.itertuples(index=False):
            column = f"edge_{edge.edge_id:05d}"
            if column not in wide:
                raise RuntimeError(f"Missing wide edge column: {column}")
            values = wide[column].to_numpy(float)
            if not np.isfinite(values).all():
                raise RuntimeError(f"Non-finite edge values: {column}")
            i, j = node_index[edge.node_i], node_index[edge.node_j]
            matrices[i, j, :] = values
            matrices[j, i, :] = values
        input_dir = STAGE / "outputs" / "inputs" / pair
        input_dir.mkdir(parents=True, exist_ok=True)
        savemat(
            input_dir / "matrices.mat", {"connectivity": matrices}, do_compression=True
        )
        np.savetxt(input_dir / "design.txt", design, fmt="%.17g")
        lookup = (
            pd.concat(
                (
                    edges[["node_i", "parcel_i", "network_i"]].rename(
                        columns={
                            "node_i": "node",
                            "parcel_i": "parcel",
                            "network_i": "network",
                        }
                    ),
                    edges[["node_j", "parcel_j", "network_j"]].rename(
                        columns={
                            "node_j": "node",
                            "parcel_j": "parcel",
                            "network_j": "network",
                        }
                    ),
                )
            )
            .drop_duplicates()
            .sort_values("node")
        )
        lookup.to_csv(input_dir / "nodes.tsv", sep="\t", index=False)
        edges.to_csv(input_dir / "mask_edges.tsv", sep="\t", index=False)
        output_rows.append(
            {
                "network_pair": pair,
                "subject_count": len(wide),
                "node_count": len(nodes),
                "edge_count_in_mask": len(edges),
                "matrix_path": str(input_dir / "matrices.mat"),
                "design_path": str(input_dir / "design.txt"),
            }
        )
    pd.DataFrame(output_rows).to_csv(
        STAGE / "outputs" / "nbs_input_summary.tsv", sep="\t", index=False
    )
    with (STAGE / "config" / "nbs_parameters.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["parameter", "value"])
        writer.writerows(
            (
                ("model", "edge ~ group + age"),
                ("group_direction", "FM - HC"),
                ("primary_threshold", f"t > {cfg['NBS_BASE_THRESHOLD']}"),
                ("component_statistic", "Intensity"),
                ("permutations", cfg["NBS_PERMUTATIONS"]),
                ("seed", cfg["PERM_SEED"]),
            )
        )


if __name__ == "__main__":
    main()
