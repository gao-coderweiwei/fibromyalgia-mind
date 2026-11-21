# Purpose: Annotate components.
# Inputs: Component maps, term matrix, PET maps, and spins.
# Outputs: Background contrasts and coverage tables.
# Manuscript: S9 and component portion of S10.

from os import environ as cfg
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

STAGE = Path(cfg["SPATIAL_STAGE"])
PROJECT = Path(cfg["FM_WORK"])
MAP_STAGE = Path(cfg["MAP_STAGE"])
NEUROSYNTH = Path(cfg["NEUROSYNTH_STAGE"]) / "outputs"
PET_REVISION = cfg["PET_COMMIT"]
PET_ROOT = (
    f"https://raw.githubusercontent.com/netneurolab/hansen_receptors/{PET_REVISION}"
)
SYSTEMS = {
    "serotonin": ["5HT1a", "5HT1b", "5HT2a", "5HT4", "5HT6", "5HTT"],
    "acetylcholine": ["A4B2", "M1", "VAChT"],
    "cannabinoid": ["CB1"],
    "dopamine": ["D1", "D2", "DAT"],
    "GABA": ["GABAa"],
    "histamine": ["H3"],
    "glutamate": ["mGluR5", "NMDA"],
    "opioid": ["MOR"],
    "norepinephrine": ["NET"],
}


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def save(table, filename):
    table.to_csv(STAGE / filename, sep="\t", index=False, float_format="%.10g")


def unit_center(matrix):
    centered = matrix - matrix.mean(axis=-1, keepdims=True)
    return centered / np.linalg.norm(centered, axis=-1, keepdims=True)


def zscore(matrix):
    return unit_center(matrix) * np.sqrt(matrix.shape[-1])


def p_spin(observed, simulated):
    return (
        1 + (np.abs(simulated) >= np.abs(observed[:, None]) - 1e-12).sum(axis=1)
    ) / (simulated.shape[1] + 1)


def residualize(values, strata):
    result = values.astype(float, copy=True)
    for stratum in np.unique(strata):
        selection = strata == stratum
        result[..., selection] -= result[..., selection].mean(axis=-1, keepdims=True)
    return result


def load_pet():
    path = STAGE / "hansen_scale200.tsv"
    if not path.exists():
        with TemporaryDirectory() as directory:
            data = Path(directory) / "receptor_data.csv"
            names = Path(directory) / "receptor_names.npy"
            urlretrieve(f"{PET_ROOT}/results/receptor_data_scale200.csv", data)
            urlretrieve(f"{PET_ROOT}/data/receptor_names_pet.npy", names)
            table = pd.DataFrame(
                np.loadtxt(data, delimiter=","),
                columns=np.load(names, allow_pickle=False),
            )
            table.insert(0, "parcel_id", np.arange(1, len(table) + 1))
            save(table, path.name)
    table = read_tsv(path).set_index("parcel_id").loc[np.arange(1, 201)]
    receptors = pd.DataFrame(zscore(table.to_numpy().T), index=table.columns)
    systems = pd.DataFrame(
        {
            system: receptors.loc[members].mean(axis=0)
            for system, members in SYSTEMS.items()
        }
    ).T
    systems.loc[:, :] = zscore(systems.to_numpy())
    return {"receptor": receptors, "system": systems}


def functional_terms(pair, terminal, degree, strength, networks, terms, labels, spins):
    eligible = np.isin(networks, ["SalVentAttn", terminal])
    strata = networks[eligible]
    profiles = labels.copy()
    profiles.insert(0, "network_pair", pair)
    profiles["r_complete_network_union"] = unit_center(terms) @ unit_center(
        eligible.astype(float)
    )
    for network in ["SalVentAttn", terminal]:
        profiles[f"term_mean_{network}"] = terms[:, networks == network].mean(axis=1)

    # Center the term and component values within each candidate network.
    within_terms = unit_center(residualize(terms[:, eligible], strata))
    within_maps = {
        name: unit_center(residualize(values[eligible], strata))
        for name, values in {"degree": degree, "strength": strength}.items()
    }
    within_null = {name: [] for name in within_maps}
    for spin_block in np.array_split(spins[eligible, :].T, 20):
        spun = unit_center(residualize(terms[:, spin_block], strata))
        for name, values in within_maps.items():
            within_null[name].append(spun @ values)

    for name, values in {"degree": degree, "strength": strength}.items():
        full_unit = unit_center(values)
        full_terms = unit_center(terms)
        full_r = full_terms @ full_unit
        full_null = full_terms @ unit_center(values[spins].T).T
        within_r = within_terms @ within_maps[name]
        profiles[f"r_whole_{name}"] = full_r
        profiles[f"p_whole_{name}"] = p_spin(full_r, full_null)
        profiles[f"q_whole_{name}"] = false_discovery_control(
            profiles[f"p_whole_{name}"]
        )
        profiles[f"r_within_{name}"] = within_r
        profiles[f"p_within_{name}"] = p_spin(
            within_r, np.concatenate(within_null[name], axis=1)
        )
        profiles[f"q_within_{name}"] = false_discovery_control(
            profiles[f"p_within_{name}"]
        )
    return profiles


def receptor_contrasts(pair, terminal, support, networks, pet, spins):
    component = support.astype(float) / support.sum()
    baseline = np.zeros(len(networks))
    for network in ["SalVentAttn", terminal]:
        selection = networks == network
        baseline[selection] = component[selection].sum() / selection.sum()
    contrast = component - baseline
    result = []
    for family, table in pet.items():
        reference = table.to_numpy()
        observed = reference @ contrast
        null = np.concatenate(
            [
                reference[:, block.T] @ contrast
                for block in np.array_split(spins, 20, axis=1)
            ],
            axis=1,
        )
        original = reference @ component
        whole_p = p_spin(original, reference @ component[spins])
        local_p = p_spin(observed, null)
        result.append(
            pd.DataFrame(
                {
                    "network_pair": pair,
                    "family": family,
                    "reference": table.index,
                    "component_mean_z": original,
                    "matched_complete_network_mean_z": reference @ baseline,
                    "difference_z": observed,
                    "p_whole": whole_p,
                    "q_whole": false_discovery_control(whole_p),
                    "p_network_contrast": local_p,
                    "q_network_contrast": false_discovery_control(local_p),
                }
            )
        )
    return pd.concat(result, ignore_index=True)


def main():
    parcels = read_tsv(MAP_STAGE / "atlas/schaefer200_parcel_lookup.tsv").sort_values(
        "parcel_id"
    )
    assert np.array_equal(parcels.parcel_id, np.arange(1, 201))
    networks = parcels.network.to_numpy()
    edges = read_tsv(MAP_STAGE / "components/nbs_fm_gt_hc_component_edges.tsv")
    maps = read_tsv(
        MAP_STAGE / "maps/nbs_component_node_maps.tsv"
    ).sort_values("parcel_id")
    terms = read_tsv(NEUROSYNTH / "term_matrix/neurosynth_schaefer200_term_matrix.tsv")
    term_values = terms[
        [f"parcel_{index:03}" for index in parcels.parcel_id]
    ].to_numpy()
    spins = np.load(NEUROSYNTH / f"spin_indices_vasa_seed{cfg['PERM_SEED']}.npy")
    assert spins.shape == (200, int(cfg["N_SPINS"]))
    pet = load_pet()
    functional, receptors, coverage, map_rows = [], [], [], []
    for terminal in ["Default", "DorsAttn"]:
        pair = f"SalVentAttn-{terminal}"
        selected = edges[edges.network_pair == pair]
        nodes = selected[["node_i", "node_j"]].to_numpy().ravel() - 1
        degree = np.bincount(nodes, minlength=200)
        strength = np.bincount(
            nodes, weights=np.repeat(np.abs(selected.t_value), 2), minlength=200
        )
        assert np.array_equal(degree, maps[f"SalVentAttn_{terminal}_component_degree"])
        assert np.allclose(strength, maps[f"SalVentAttn_{terminal}_component_strength"])
        support = degree > 0
        eligible = np.isin(networks, ["SalVentAttn", terminal])
        eligible_edges = int(
            (networks == "SalVentAttn").sum() * (networks == terminal).sum()
        )
        coverage.append(
            {
                "network_pair": pair,
                "edges": len(selected),
                "candidate_edges": eligible_edges,
                "edge_fraction": len(selected) / eligible_edges,
                "nodes": support.sum(),
                "candidate_nodes": eligible.sum(),
                "node_fraction": support.sum() / eligible.sum(),
                "salventattn_nodes": (support & (networks == "SalVentAttn")).sum(),
                "partner_nodes": (support & (networks == terminal)).sum(),
            }
        )
        map_rows.append(
            parcels.assign(
                network_pair=pair,
                in_candidate=eligible,
                degree=degree,
                strength=strength,
            )
        )
        functional.append(
            functional_terms(
                pair,
                terminal,
                degree,
                strength,
                networks,
                term_values,
                terms[["term", "category"]],
                spins,
            )
        )
        receptors.append(
            receptor_contrasts(pair, terminal, support, networks, pet, spins)
        )
    save(pd.DataFrame(coverage), "component_coverage.tsv")
    save(pd.concat(map_rows), "component_maps.tsv")
    save(pd.concat(functional), "functional_background.tsv")
    save(pd.concat(receptors), "receptor_background.tsv")


if __name__ == "__main__":
    main()
