"""Marker observability (G1) and the two donor-level tables the text quotes.

Writes
    data/detection_astrocyte.csv   per-marker detection, donors clearing the
                                   floor, and the G1 verdict
    data/hmox1_donors.csv          per-donor HMOX1 detection with the donor
                                   ambient index
    data/bin1_donors.csv           the first bin's donors with their metadata

These three were previously carried in the repository as files with no script
behind them. They back numbers in the text, so they are generated here from the
deposited inputs rather than shipped as results somebody has to take on trust.

The only inputs are inputs/detection_by_donor_Astrocyte.csv (per-donor, per-gene
detection rate), inputs/donor_mean_Astrocyte.csv and inputs/donor_meta.csv. The
44-marker panel and its categories are fixed in advance and written out below,
because the panel is an input of the study, not a result of it.
"""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from paths import DATA, INPUT

FLOOR = 0.02            # the floor the engine applies to background genes
MIN_DONORS = 60         # of 84
AMBIENT = ["P2RY12", "PLP1", "RBFOX3"]   # cell-type-exclusive markers, read in astrocytes

PANEL = {
    "Reactive astrocyte": ["GFAP", "SERPINA3", "VIM", "FBLN5", "AQP4", "MERTK",
                           "MFGE8", "SOX9", "FKBP5"],
    "A1 signaling": ["C3", "STAT1", "STAT3", "SOCS3", "NFKB1"],
    "Inflammation": ["IL1B", "IL6", "TNF", "CCL2", "CXCL12", "CXCL16"],
    "Vascular": ["ICAM1", "VCAM1", "VWF", "SELE", "OCLN", "PECAM1"],
    "AD risk": ["APOE", "CLU", "TREM2", "CR1", "CD33"],
    "Iron/oxidative": ["FTH1", "FTL", "HMOX1", "NFE2L2", "TFRC"],
    "Injury/protective": ["NEFL", "BDNF", "GDNF", "NGF"],
    "AD core/cell death": ["APP", "CASP3", "BCL2", "TP53"],
}
GENE2CAT = {g: c for c, gs in PANEL.items() for g in gs}
PANEL44 = [g for gs in PANEL.values() for g in gs]
assert len(PANEL44) == 44 and len(set(PANEL44)) == 44

if __name__ == "__main__":
    dd = pd.read_csv(INPUT + "detection_by_donor_Astrocyte.csv")
    meta = pd.read_csv(INPUT + "donor_meta.csv")
    mn = pd.read_csv(INPUT + "donor_mean_Astrocyte.csv")
    mn = mn[~mn.cps.isna()].reset_index(drop=True)
    DONOR = "donor" if "donor" in mn.columns else "donor_id"

    nc = dd.n_cells.to_numpy(float)
    genes = [c for c in dd.columns if c not in ("donor_id", "cps", "n_cells")]

    # donor-level ambient index: cell-type-exclusive markers read inside astrocytes
    amb = dd[AMBIENT].mean(axis=1)
    wmean = lambda g: float((dd[g].to_numpy(float) * nc).sum() / nc.sum())
    # The band is the nucleus-weighted detection of those markers, on the same
    # scale as every other detection rate here, not an unweighted donor mean.
    AMB_LO, AMB_HI = min(map(wmean, AMBIENT)), max(map(wmean, AMBIENT))

    rows = []
    for g in genes:
        v = dd[g].to_numpy(float)
        det = wmean(g)
        npass = int((v >= FLOOR).sum())
        # Order matters and is the order the gates are applied in: a marker that
        # too few donors carry is reported as such even when it also sits in the
        # ambient band. The band's UPPER bound is what the verdict uses.
        if det < FLOOR:
            verdict = "below floor"
        elif npass < MIN_DONORS:
            verdict = "too few donors"
        elif det <= AMB_HI:
            verdict = "at/below ambient"
        else:
            verdict = "observable"
        rows.append({"gene": g, "detection": round(det, 5),
                     "donors_passing_floor": npass, "verdict": verdict,
                     "in_panel": g in GENE2CAT,
                     "category": GENE2CAT.get(g, "control")})
    det = pd.DataFrame(rows)
    det.to_csv(DATA + "detection_astrocyte.csv", index=False)
    print(f"detection_astrocyte.csv: {len(det)} genes | ambient band "
          f"{AMB_LO:.5f}-{AMB_HI:.5f}")
    print("  " + " | ".join(f"{k} {v}" for k, v in det.verdict.value_counts().items()))

    # ---- HMOX1 by donor -----------------------------------------------------
    z = mn["HMOX1"].to_numpy(float)
    h = pd.DataFrame({"donor": mn[DONOR], "cps": mn.cps, "n_astro": mn.n_cells.astype(int),
                      "z": (z - z.mean()) / z.std(ddof=0)})
    h = h.merge(dd[["donor_id", "HMOX1"]].rename(columns={"donor_id": "donor",
                                                          "HMOX1": "detection"}),
                on="donor", how="left")
    h = h.merge(pd.DataFrame({"donor": dd.donor_id, "ambient": amb}), on="donor", how="left")
    h.to_csv(DATA + "hmox1_donors.csv", index=False)
    print(f"hmox1_donors.csv: {len(h)} donors | cohort median detection "
          f"{h.detection.median() * 100:.3f}% | max {h.detection.max() * 100:.2f}%")

    # ---- the first bin ------------------------------------------------------
    b1 = mn[mn.cps < 0.2][[DONOR, "cps", "n_cells"]].rename(
        columns={DONOR: "donor_id", "n_cells": "n_Astrocyte"})
    b1 = b1.merge(meta, on="donor_id", how="left")
    b1 = b1.merge(h[["donor", "detection", "ambient"]].rename(
        columns={"donor": "donor_id", "detection": "hmox1_det"}), on="donor_id", how="left")
    b1.to_csv(DATA + "bin1_donors.csv", index=False)
    print(f"bin1_donors.csv: {len(b1)} donors | HMOX1 "
          + ", ".join(f"{x * 100:.2f}%" for x in b1.hmox1_det))
