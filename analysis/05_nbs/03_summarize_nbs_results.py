#!/usr/bin/env python3
# Purpose: 03 summarize nbs results.
# Inputs: Baseline MAT files, edge lookup, and corrected threshold table.
# Outputs: Component edge, node, and summary tables.
# Manuscript: S3.

from __future__ import annotations

import subprocess
import sys
from os import environ as cfg
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(cfg["FM_WORK"])
STAGE = Path(cfg["NBS_STAGE"])
NETWORK = Path(cfg["MIND_STAGE"])
PAIRS = ("SalVentAttn-Default", "SalVentAttn-DorsAttn")
DIRECTIONS = ("FM_gt_HC", "HC_gt_FM")


def mat_cells(value: object) -> list[object]:
    if isinstance(value, np.ndarray) and value.dtype == object:
        return list(value.ravel())
    if value is None or (isinstance(value, np.ndarray) and value.size == 0):
        return []
    return [value]


def main() -> None:
    edge_table = pd.read_csv(
        NETWORK / "outputs" / "tables" / "mind_edge_dictionary.tsv", sep="\t"
    )
    edge_table["edge_id"] = pd.to_numeric(edge_table["edge_id"], errors="raise").astype(
        int
    )
    masks = pd.read_csv(
        NETWORK / "outputs" / "masks" / "analysis_network_masks.tsv",
        sep="\t",
        dtype=str,
    )
    input_summary = pd.read_csv(STAGE / "outputs" / "nbs_input_summary.tsv", sep="\t")
    summary_rows, edge_rows, node_rows = [], [], []
    raw_available = True
    for pair in PAIRS:
        pair_input = input_summary.loc[input_summary["network_pair"] == pair].iloc[0]
        nodes = pd.read_csv(STAGE / "outputs" / "inputs" / pair / "nodes.tsv", sep="\t")
        local_to_global = nodes["node"].to_numpy(int)
        local_lookup = nodes.set_index("node")
        for direction in DIRECTIONS:
            raw_path = STAGE / "outputs" / "raw" / f"{pair}_{direction}.mat"
            if not raw_path.is_file():
                raw_available = False
                continue
            result = loadmat(raw_path, squeeze_me=True)
            test_stat = np.asarray(result["test_stat"], dtype=float)
            component_cells = mat_cells(result.get("component_matrices"))
            pvalues = (
                np.atleast_1d(result.get("component_pvalues", np.array([])))
                .astype(float)
                .ravel()
            )
            if not component_cells:
                summary_rows.append(
                    {
                        "network_pair": pair,
                        "direction": direction,
                        "edge_count_in_mask": int(pair_input.edge_count_in_mask),
                        "primary_threshold": f"t > {cfg['NBS_BASE_THRESHOLD']}",
                        "component_statistic": "Intensity",
                        "permutations": int(cfg["NBS_PERMUTATIONS"]),
                        "seed": int(cfg["PERM_SEED"]),
                        "component_id": np.nan,
                        "component_edge_count": np.nan,
                        "component_node_count": np.nan,
                        "component_mass": np.nan,
                        "p_fwe": 1.0,
                        "p_fwe_across_masks_directions": np.nan,
                        "software": "NBS",
                        "software_version": "1.2",
                    }
                )
            for component_index, component in enumerate(component_cells, start=1):
                adjacency = np.asarray(
                    component.toarray() if hasattr(component, "toarray") else component,
                    dtype=float,
                )
                local_i, local_j = np.where(np.triu(adjacency, 1) > 0)
                component_edges = []
                for i, j in zip(local_i, local_j):
                    global_i, global_j = int(local_to_global[i]), int(
                        local_to_global[j]
                    )
                    match = edge_table.loc[
                        (
                            (edge_table["node_i"] == global_i)
                            & (edge_table["node_j"] == global_j)
                        )
                        | (
                            (edge_table["node_i"] == global_j)
                            & (edge_table["node_j"] == global_i)
                        )
                    ]
                    if len(match) != 1:
                        raise RuntimeError(
                            "Component edge does not map to the MIND edge dictionary"
                        )
                    edge = match.iloc[0]
                    t_value = float(test_stat[i, j])
                    component_edges.append((edge, t_value, i, j))
                    edge_rows.append(
                        {
                            "network_pair": pair,
                            "direction": direction,
                            "component_id": component_index,
                            "edge_id": int(edge.edge_id),
                            "node_i": int(edge.node_i),
                            "node_j": int(edge.node_j),
                            "parcel_i": edge.parcel_i,
                            "parcel_j": edge.parcel_j,
                            "network_i": edge.network_i,
                            "network_j": edge.network_j,
                            "t_value": t_value,
                        }
                    )
                component_nodes = sorted(
                    set(
                        [item[2] for item in component_edges]
                        + [item[3] for item in component_edges]
                    )
                )
                for node_index in component_nodes:
                    degree = sum(
                        node_index in (item[2], item[3]) for item in component_edges
                    )
                    global_node = int(local_to_global[node_index])
                    node_rows.append(
                        {
                            "network_pair": pair,
                            "direction": direction,
                            "component_id": component_index,
                            "node": global_node,
                            "parcel": local_lookup.loc[global_node, "parcel"],
                            "network": local_lookup.loc[global_node, "network"],
                            "degree_in_component": degree,
                        }
                    )
                mass = float(
                    sum(
                        value - float(cfg["NBS_BASE_THRESHOLD"])
                        for _, value, _, _ in component_edges
                    )
                )
                p_fwe = (
                    float(pvalues[component_index - 1])
                    if component_index <= len(pvalues)
                    else np.nan
                )
                summary_rows.append(
                    {
                        "network_pair": pair,
                        "direction": direction,
                        "edge_count_in_mask": int(pair_input.edge_count_in_mask),
                        "primary_threshold": f"t > {cfg['NBS_BASE_THRESHOLD']}",
                        "component_statistic": "Intensity",
                        "permutations": int(cfg["NBS_PERMUTATIONS"]),
                        "seed": int(cfg["PERM_SEED"]),
                        "component_id": component_index,
                        "component_edge_count": len(component_edges),
                        "component_node_count": len(component_nodes),
                        "component_mass": mass,
                        "p_fwe": p_fwe,
                        "p_fwe_across_masks_directions": np.nan,
                        "software": "NBS",
                        "software_version": "1.2",
                    }
                )
    if not raw_available:
        raise RuntimeError("One or more official NBS output files are missing")
    results = pd.DataFrame(
        summary_rows,
        columns=[
            "network_pair",
            "direction",
            "edge_count_in_mask",
            "primary_threshold",
            "component_statistic",
            "permutations",
            "seed",
            "component_id",
            "component_edge_count",
            "component_node_count",
            "component_mass",
            "p_fwe",
            "p_fwe_across_masks_directions",
            "software",
            "software_version",
        ],
    )
    # Read the fixed four-test Bonferroni correction from the threshold analysis.
    thresholds = pd.read_csv(
        Path(cfg["THRESHOLD_STAGE"]) / "nbs_threshold_results.tsv", sep="\t"
    )
    keys = ["network_pair", "direction", "component_id"]
    corrected = thresholds.loc[
        (thresholds["threshold"] == float(cfg["NBS_BASE_THRESHOLD"]))
        & thresholds["component_id"].notna()
    ]
    observed = results["component_id"].notna()
    matched = results.loc[observed, keys + ["p_fwe"]].merge(
        corrected,
        on=keys,
        how="left",
        validate="one_to_one",
        suffixes=("", "_threshold"),
    )
    if len(matched) != len(corrected) or not np.allclose(
        matched["p_fwe"], matched["p_fwe_threshold"], rtol=0, atol=1e-12
    ):
        raise ValueError("Threshold table does not match the baseline NBS components")
    results.loc[observed, "p_fwe_across_masks_directions"] = matched[
        "p_bonferroni_within_threshold"
    ].to_numpy()
    results.to_csv(STAGE / "outputs" / "nbs_results_summary.tsv", sep="\t", index=False)
    pd.DataFrame(
        edge_rows,
        columns=[
            "network_pair",
            "direction",
            "component_id",
            "edge_id",
            "node_i",
            "node_j",
            "parcel_i",
            "parcel_j",
            "network_i",
            "network_j",
            "t_value",
        ],
    ).to_csv(STAGE / "outputs" / "nbs_component_edges.tsv", sep="\t", index=False)
    pd.DataFrame(
        node_rows,
        columns=[
            "network_pair",
            "direction",
            "component_id",
            "node",
            "parcel",
            "network",
            "degree_in_component",
        ],
    ).to_csv(STAGE / "outputs" / "nbs_component_nodes.tsv", sep="\t", index=False)
    target_counts = masks.loc[masks["mask_name"].isin(PAIRS), "mask_name"].value_counts()
    input_paths = [STAGE / "outputs" / "inputs" / pair / "matrices.mat" for pair in PAIRS]
    checks = pd.DataFrame(
        [
            ("target_masks_present", target_counts.to_dict() == dict.fromkeys(PAIRS, 1)),
            (
                "nbs_input_matrices_finite",
                all(
                    path.is_file() and np.isfinite(loadmat(path)["connectivity"]).all()
                    for path in input_paths
                ),
            ),
        ],
        columns=["check", "status"],
    )
    checks["status"] = np.where(checks["status"], "passed", "failed")
    checks.to_csv(STAGE / "outputs" / "analysis_checks.tsv", sep="\t", index=False)
    matlab_version = (
        subprocess.check_output(
            [cfg["MATLAB_BIN"], "-batch", "disp(version)"],
            text=True,
            stderr=subprocess.STDOUT,
        )
        .splitlines()[-1]
        .strip()
    )
    pd.DataFrame(
        {
            "software": ["NBS", "MATLAB", "nbs_toolbox_path"],
            "version_or_path": ["1.2", matlab_version, cfg["NBS_HOME"]],
        }
    ).to_csv(STAGE / "software" / "nbs_software_versions.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "software": ["python", "numpy", "pandas", "scipy", "MATLAB executable"],
            "version_or_path": [
                sys.version.split()[0],
                np.__version__,
                pd.__version__,
                __import__("scipy").__version__,
                cfg["MATLAB_BIN"],
            ],
        }
    ).to_csv(STAGE / "software" / "software_versions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
