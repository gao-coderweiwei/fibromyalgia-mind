# Purpose: Handle degenerate MRIQC QI2 fits.
# Inputs: Anatomical image and air mask.
# Outputs: QI2 value and diagnostic status.
# Manuscript: MRI acquisition and processing; S1.

"""Runtime MRIQC patch for degenerate anatomical QI2 fits.

Some anatomical scans can yield an all-zero/non-finite background distribution
inside MRIQC's air mask. MRIQC 24.1.0.dev0 then raises from scikit-learn's
KernelDensity and aborts the whole dataset. Returning NaN keeps the rest of the
IQMs and reports available while making the failed metric explicit.
"""

from __future__ import annotations

import os

import numpy as np


def _patch_mriqc_art_qi2() -> None:
    try:
        import mriqc.interfaces.anatomical as anatomical_ifaces
        import mriqc.qc.anatomical as anatomical_qc
    except Exception:
        return

    original_art_qi2 = anatomical_qc.art_qi2

    def safe_art_qi2(
        img,
        airmask,
        min_voxels=int(1e3),
        max_voxels=int(3e5),
        save_plot=True,
        coil_elements=32,
    ):
        data = np.nan_to_num(img[airmask > 0], posinf=0.0, neginf=0.0)
        data[data < 0] = 0
        out_file = os.path.abspath("error.svg")

        if (data > 0).sum() < min_voxels:
            return 0.0, out_file

        p99 = np.percentile(data, 99)
        if not np.isfinite(p99) or p99 <= 0:
            with open(out_file, "w", encoding="utf-8") as ofh:
                ofh.write(
                    "<p>Background noise fitting was skipped: degenerate air-mask intensities.</p>"
                )
            return float("nan"), out_file

        try:
            return original_art_qi2(
                img,
                airmask,
                min_voxels=min_voxels,
                max_voxels=max_voxels,
                save_plot=save_plot,
                coil_elements=coil_elements,
            )
        except ValueError as exc:
            if "Input X contains NaN" not in str(exc):
                raise
            with open(out_file, "w", encoding="utf-8") as ofh:
                ofh.write(
                    "<p>Background noise fitting was skipped: non-finite KDE input.</p>"
                )
            return float("nan"), out_file

    anatomical_qc.art_qi2 = safe_art_qi2
    anatomical_ifaces.art_qi2 = safe_art_qi2


_patch_mriqc_art_qi2()
