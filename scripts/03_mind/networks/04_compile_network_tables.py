#!/usr/bin/env python3
# Purpose: 04 compile network tables.
# Inputs: Subject matrices and parcel lookup.
# Outputs: Edge tables, masks, and checks.
# Manuscript: MIND construction; S2.

"""Compile subject matrices into edge tables and predefined network masks."""

from os import environ as cfg
from pathlib import Path

import numpy as np
import pandas as pd

stage = Path(cfg["MIND_STAGE"])
root = Path(cfg["FM_WORK"])
lookup = pd.read_csv(
    Path(cfg["ATLAS_STAGE"]) / "outputs" / "schaefer200_parcel_lookup.tsv", sep="\t"
)
mapping = pd.read_csv(
    Path(cfg["ATLAS_STAGE"]) / "outputs" / "ds004144_subject_freesurfer_mapping.tsv",
    sep="\t",
)
mapping = mapping.loc[mapping["mapping_status"] == "complete"].copy()
node_names = lookup["parcel_name"].tolist()
upper_i, upper_j = np.triu_indices(len(lookup), k=1)
mask_pairs = [
    "SalVentAttn-Default",
    "SalVentAttn-DorsAttn",
    "Default-DorsAttn",
    "SalVentAttn-SomMot",
    "Default-SomMot",
    "DorsAttn-SomMot",
    "within SalVentAttn",
    "within Default",
    "within DorsAttn",
]


def pair_label(first, second):
    for label in mask_pairs[:6]:
        left, right = label.split("-")
        if {first, second} == {left, right}:
            return label
    return f"{first}-{second}"


edge_rows = []
for edge_id, (i, j) in enumerate(zip(upper_i, upper_j), start=1):
    network_i, network_j = lookup.iloc[i]["network"], lookup.iloc[j]["network"]
    edge_rows.append(
        {
            "edge_id": edge_id,
            "node_i": i + 1,
            "node_j": j + 1,
            "parcel_i": lookup.iloc[i]["parcel_name"],
            "parcel_j": lookup.iloc[j]["parcel_name"],
            "network_i": network_i,
            "network_j": network_j,
            "network_pair": pair_label(network_i, network_j),
        }
    )
edges = pd.DataFrame(edge_rows)
edges.to_csv(
    stage / "outputs" / "tables" / "mind_edge_dictionary.tsv", sep="\t", index=False
)

wide_rows = []
summary_rows = []
for _, row in mapping.iterrows():
    subject = row["subject_id"]
    matrix = pd.read_csv(
        stage / "outputs" / "matrices" / f"{subject}_mind_matrix.tsv",
        sep="\t",
        index_col=0,
    )
    matrix = matrix.reindex(index=node_names, columns=node_names)
    values = matrix.to_numpy()[upper_i, upper_j]
    wide = {
        "cohort": row["cohort"],
        "subject_id": subject,
        "group": row["group"],
        "sex": row["sex"],
        "age": row["age"],
    }
    wide.update(
        {f"edge_{edge_id:05d}": value for edge_id, value in enumerate(values, start=1)}
    )
    wide_rows.append(wide)
    summary_rows.append(
        {
            "cohort": row["cohort"],
            "subject_id": subject,
            "group": row["group"],
            "sex": row["sex"],
            "age": row["age"],
            "edge_count": len(values),
            "mind_min": np.min(values),
            "mind_max": np.max(values),
            "mind_mean": np.mean(values),
            "finite_edge_values": bool(np.isfinite(values).all()),
        }
    )
pd.DataFrame(wide_rows).to_csv(
    stage / "outputs" / "tables" / "mind_edges_wide.tsv",
    sep="\t",
    index=False,
    float_format="%.17g",
)
pd.DataFrame(summary_rows).to_csv(
    stage / "outputs" / "tables" / "subject_mind_summary.tsv",
    sep="\t",
    index=False,
    float_format="%.17g",
)

mask_rows = []
for label in mask_pairs:
    if label.startswith("within "):
        network = label.removeprefix("within ")
        selected = (edges["network_i"] == network) & (edges["network_j"] == network)
    else:
        left, right = label.split("-")
        selected = ((edges["network_i"] == left) & (edges["network_j"] == right)) | (
            (edges["network_i"] == right) & (edges["network_j"] == left)
        )
    selected_edges = edges.loc[selected, "edge_id"].tolist()
    mask_rows.append(
        {
            "mask_name": label,
            "edge_count": len(selected_edges),
            "edge_ids": ",".join(map(str, selected_edges)),
        }
    )
pd.DataFrame(mask_rows).to_csv(
    stage / "outputs" / "masks" / "analysis_network_masks.tsv", sep="\t", index=False
)
pd.DataFrame(mask_rows)[["mask_name", "edge_count"]].to_csv(
    stage / "outputs" / "masks" / "network_pair_edge_counts.tsv", sep="\t", index=False
)

checks = []


def add(name, passed, details):
    checks.append(
        {"check": name, "status": "pass" if passed else "fail", "details": details}
    )


matrices = []
for subject in mapping["subject_id"]:
    matrix = (
        pd.read_csv(
            stage / "outputs" / "matrices" / f"{subject}_mind_matrix.tsv",
            sep="\t",
            index_col=0,
        )
        .reindex(index=node_names, columns=node_names)
        .to_numpy()
    )
    matrices.append(matrix)
add("subject_count", len(mapping) == 66, str(len(mapping)))
add(
    "group_counts",
    mapping["group"].value_counts().to_dict() == {"FM": 33, "HC": 33},
    str(mapping["group"].value_counts().to_dict()),
)
add(
    "annotation_files_created",
    all(
        (
            stage
            / "outputs"
            / "annotations"
            / subject
            / "label"
            / f"{hemi}.Schaefer2018_200Parcels_7Networks_order.annot"
        ).is_file()
        for subject in mapping["subject_id"]
        for hemi in ("lh", "rh")
    ),
    "two annotations per subject",
)
add("mind_matrices_created", len(matrices) == len(mapping), str(len(matrices)))
add(
    "matrix_dimension_200_by_200",
    all(matrix.shape == (200, 200) for matrix in matrices),
    "200 by 200",
)
add(
    "matrix_symmetric_atol_1e-10",
    all(np.allclose(matrix, matrix.T, atol=1e-10, rtol=0) for matrix in matrices),
    "atol=1e-10",
)
add(
    "finite_mind_values",
    all(np.isfinite(matrix).all() for matrix in matrices),
    "all matrices",
)
add("edge_count_19900", len(edges) == 19900, str(len(edges)))
add(
    "parcel_order_matches_lookup",
    all(
        list(
            pd.read_csv(
                stage / "outputs" / "matrices" / f"{subject}_mind_matrix.tsv",
                sep="\t",
                index_col=0,
            ).index
        )
        == node_names
        for subject in mapping["subject_id"]
    ),
    "lookup order",
)
add(
    "network_labels_present",
    set(
        ["Visual", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"]
    ).issubset(set(lookup["network"])),
    "Yeo7",
)
add(
    "primary_masks_nonempty",
    all(row["edge_count"] > 0 for row in mask_rows[:2]),
    "SalVentAttn pairs",
)
add(
    "primary_mask_edge_counts",
    all(row["edge_count"] > 0 for row in mask_rows[:2]),
    "; ".join(f"{row['mask_name']}={row['edge_count']}" for row in mask_rows[:2]),
)
pd.DataFrame(checks).to_csv(
    stage / "outputs" / "checks" / "mind_network_checks.tsv", sep="\t", index=False
)
