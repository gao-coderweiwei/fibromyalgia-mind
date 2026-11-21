# Purpose: Analyze.
# Inputs: Smoothed maps, covariates, and atlas maps.
# Outputs: Voxel inference, peaks, and overlap table.
# Manuscript: S6.

import argparse
from os import environ as cfg
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.processing import resample_from_to
from preprocess import PRIMARY, SPM, STAGE


def intracranial_volume(subject):
    path = Path(cfg["FS_DIR"]) / subject / "stats/aseg.stats"
    line = next(
        line
        for line in path.read_text().splitlines()
        if line.startswith("# Measure EstimatedTotalIntraCranialVol,")
    )
    return float(line.split(",")[3]) / 1e6


def voxel_glm(outcomes, nuisance, group, repetitions):
    reduced, _ = np.linalg.qr(nuisance)
    adjusted_group = group - reduced @ (reduced.T @ group)
    magnitude = np.linalg.norm(adjusted_group)
    basis = np.column_stack((adjusted_group / magnitude, reduced))
    residual = outcomes - reduced @ (reduced.T @ outcomes)
    residual_df = len(group) - basis.shape[1]
    rng = np.random.default_rng(int(cfg["PERM_SEED"]))
    orders = np.stack(
        [
            np.arange(len(group)),
            *(rng.permutation(len(group)) for _ in range(repetitions)),
        ]
    )
    weights = basis[np.argsort(orders)].transpose(0, 2, 1)
    t_values = np.empty(outcomes.shape[1])
    betas = np.empty_like(t_values)
    maxima = np.zeros(repetitions)
    for start in range(0, outcomes.shape[1], 4096):
        section = slice(start, start + 4096)
        data = residual[:, section]
        total = np.einsum("ij,ij->j", data, data)
        if np.any(total <= 0):
            raise ValueError(
                "Constant smoothed gray matter values inside inference mask"
            )
        for offset in range(0, len(weights), 256):
            subset = weights[offset : offset + 256]
            projected = (subset.reshape(-1, len(group)) @ data).reshape(
                len(subset), basis.shape[1], -1
            )
            rss = total - (projected**2).sum(axis=1)
            statistics = projected[:, 0, :] / np.sqrt(
                np.maximum(rss, 1e-20) / residual_df
            )
            if offset == 0:
                t_values[section] = statistics[0]
                betas[section] = projected[0, 0] / magnitude
                null = slice(0, len(subset) - 1)
                values = statistics[1:]
            else:
                null = slice(offset - 1, offset + len(subset) - 1)
                values = statistics
            if len(values):
                maxima[null] = np.maximum(maxima[null], np.abs(values).max(axis=1))
    p_fwe = (
        1
        + repetitions
        - np.searchsorted(np.sort(maxima), np.abs(t_values), side="left")
    ) / (1 + repetitions)
    return betas, t_values, p_fwe


def save_map(name, values, mask, image):
    volume = np.zeros(mask.shape, dtype=np.float32)
    volume[mask] = values
    nib.save(nib.Nifti1Image(volume, image.affine), STAGE / f"{name}.nii.gz")


def qc_mosaic(images, subjects):
    figure, axes = plt.subplots(11, 6, figsize=(12, 15.4), layout="constrained")
    for axis, data, subject in zip(axes.flat, images, subjects):
        canonical = nib.as_closest_canonical(data).get_fdata(dtype=np.float32)
        world = nib.affines.apply_affine(
            np.linalg.inv(nib.as_closest_canonical(data).affine), [0, 0, 24]
        )
        axis.imshow(
            canonical[:, :, round(world[2])].T,
            origin="lower",
            cmap="gray",
            vmin=0,
            vmax=0.8,
        )
        axis.set_title(subject, fontsize=6, pad=0.4)
        axis.axis("off")
    figure.savefig(STAGE / "gm_registration_qc.png", dpi=120)
    plt.close(figure)


def report_significant_voxels():
    atlas_dir = Path(cfg["MAP_STAGE"]) / "atlas"
    atlas = nib.load(
        atlas_dir / "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"
    )
    labels = pd.read_csv(
        atlas_dir / "schaefer200_parcel_lookup.tsv", sep="\t"
    ).set_index("parcel_id")
    components = pd.read_csv(
        Path(cfg["SPATIAL_STAGE"]) / "component_maps.tsv", sep="\t"
    )
    statistic = nib.load(STAGE / "vbm_t.nii.gz")
    probabilities = nib.load(STAGE / "vbm_p_fwe.nii.gz").get_fdata()
    mask = nib.load(STAGE / "vbm_mask.nii.gz").get_fdata() > 0
    selected = mask & (probabilities < 0.05)
    coordinates = nib.affines.apply_affine(statistic.affine, np.argwhere(selected))
    parcel_ids = resample_from_to(
        atlas, (statistic.shape, statistic.affine), order=0
    ).get_fdata()[selected]
    table = pd.DataFrame(coordinates, columns=["x_mni_mm", "y_mni_mm", "z_mni_mm"])
    table["t_fm_minus_hc"] = statistic.get_fdata()[selected]
    table["p_fwe"] = probabilities[selected]
    table["parcel_id"] = parcel_ids.astype(int)
    table["parcel"] = table.parcel_id.map(labels.parcel_name)
    table["network"] = table.parcel_id.map(labels.network)
    for pair, maps in components.groupby("network_pair"):
        table[f"in_{pair}_component"] = table.parcel_id.isin(
            maps.loc[maps.degree > 0, "parcel_id"]
        )
    table = table.sort_values("p_fwe")
    table.to_csv(STAGE / "significant_voxel_overlap.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("work", type=Path)
    parser.add_argument(
        "--permutations", type=int, default=int(cfg["VBM_PERMUTATIONS"])
    )
    args = parser.parse_args()
    table = pd.read_csv(PRIMARY, sep="\t")
    images = [
        nib.load(args.work / subject / "gm_s8.nii") for subject in table.subject_id
    ]
    image = images[0]
    if not all(
        other.shape == image.shape and np.allclose(other.affine, image.affine)
        for other in images
    ):
        raise ValueError("Inconsistent normalized image grid")
    tpm = nib.funcs.four_to_three(nib.load(SPM / "tpm/TPM.nii"))[0]
    mask = resample_from_to(
        tpm, (image.shape, image.affine), order=1
    ).get_fdata() >= float(cfg["VBM_TPM_MIN"])
    data = np.stack(
        [other.get_fdata(dtype=np.float32)[mask] for other in images]
    ).astype(np.float64)
    if not np.isfinite(data).all():
        raise ValueError("Non-finite gray matter data inside mask")
    template = data.mean(axis=0)
    similarity = [np.corrcoef(subject, template)[0, 1] for subject in data]
    diagnostics = pd.DataFrame(
        {"subject_id": table.subject_id, "gm_profile_r_to_sample_mean": similarity}
    )
    diagnostics.to_csv(STAGE / "registration_qc.tsv", sep="\t", index=False)
    qc_mosaic(images, table.subject_id)

    covariates = np.column_stack(
        [table.age, [intracranial_volume(subject) for subject in table.subject_id]]
    )
    nuisance = np.column_stack(
        [np.ones(len(table)), covariates - covariates.mean(axis=0)]
    )
    beta, t_value, p_fwe = voxel_glm(
        data, nuisance, table.group_code.to_numpy(), args.permutations
    )
    locations = np.argwhere(mask)
    peaks = []
    for direction, pick in [("FM_gt_HC", np.argmax), ("HC_gt_FM", np.argmin)]:
        index = pick(t_value)
        x_mm, y_mm, z_mm = nib.affines.apply_affine(image.affine, locations[index])
        sign = 1 if direction == "FM_gt_HC" else -1
        peaks.append(
            {
                "direction": direction,
                "n": len(table),
                "mask_voxels": mask.sum(),
                "permutations": args.permutations,
                "peak_x_mni_mm": x_mm,
                "peak_y_mni_mm": y_mm,
                "peak_z_mni_mm": z_mm,
                "peak_beta_fm_minus_hc": beta[index],
                "peak_t_fm_minus_hc": t_value[index],
                "peak_p_fwe": p_fwe[index],
                "significant_voxels": ((p_fwe < 0.05) & (sign * t_value > 0)).sum(),
            }
        )
    pd.DataFrame(peaks).to_csv(STAGE / "vbm_summary.tsv", sep="\t", index=False)
    for name, values in [
        ("vbm_mask", np.ones(mask.sum())),
        ("vbm_beta", beta),
        ("vbm_t", t_value),
        ("vbm_p_fwe", p_fwe),
    ]:
        save_map(name, values, mask, image)
    report_significant_voxels()


if __name__ == "__main__":
    main()
