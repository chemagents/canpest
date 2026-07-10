#!/usr/bin/env python3
"""Reproduce Entry 44 RMT residue selection for the Entry 53 veto rerun.

The script intentionally mirrors LOG/44_rmt_residue_selection/rmt_select_full.py:
RMT prior s_i is computed on the full labelled dataset, point-biserial
correlation r_pb is computed against pesticide activity, residues are ranked by
|r_pb| * s_i, and m is chosen by the same inner PR-AUC scan.

Run from the repository root:
    python3 select_rmt_residues.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pointbiserialr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedShuffleSplit


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
BUNDLE_DATA_DIR = SCRIPT_DIR / "data"
DATA_DIR = BUNDLE_DATA_DIR if (BUNDLE_DATA_DIR / "data.csv").exists() else PROJECT_ROOT / "data"
OUT_DIR = SCRIPT_DIR / "rmt_selection"
RESIDUE_DIR = SCRIPT_DIR / "source_inputs" / "residue_matrices"
if not RESIDUE_DIR.exists():
    RESIDUE_DIR = DATA_DIR / "residue_matrices"

SOURCE_CODE_DIR = SCRIPT_DIR / "source_code"
if (SOURCE_CODE_DIR / "rmt_filter.py").exists():
    sys.path.insert(0, str(SOURCE_CODE_DIR))
else:
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from rmt_filter import load_residue_matrix, normalize_ligand_id, rmt_prior_from_train  # noqa: E402


PROTEINS = [
    {"gene": "GLC1", "key": "glc1", "top_m_expected": 45},
    {"gene": "ACHE", "key": "D8V7J0", "top_m_expected": 50},
]
MAX_M = 50
INNER_REPEATS = 3
INNER_TRAIN_FRAC = 0.7
SEED = 42
ALPHA = 1.0


def point_biserial_by_feature(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros(x.shape[1], dtype=float)
    for i in range(x.shape[1]):
        xi = x[:, i]
        if np.std(xi) == 0:
            out[i] = 0.0
            continue
        r, _ = pointbiserialr(xi, y)
        out[i] = r if not np.isnan(r) else 0.0
    return out


def select_m_inner_cv_pr_auc(
    x_full: np.ndarray,
    y_full: np.ndarray,
    alpha: float = ALPHA,
    max_m: int = MAX_M,
    repeats: int = INNER_REPEATS,
    inner_train_frac: float = INNER_TRAIN_FRAC,
    seed: int = SEED,
) -> tuple[int, pd.DataFrame]:
    _, p = x_full.shape
    max_m_eff = min(max_m, p)
    pr_sum = np.zeros(max_m_eff, dtype=float)
    pr_sq_sum = np.zeros(max_m_eff, dtype=float)
    valid_n = np.zeros(max_m_eff, dtype=int)

    for rep in range(repeats):
        sss = StratifiedShuffleSplit(n_splits=1, test_size=1 - inner_train_frac, random_state=seed + rep)
        for tr_idx, va_idx in sss.split(x_full, y_full):
            x_tr = x_full[tr_idx]
            x_va = x_full[va_idx]
            y_tr = y_full[tr_idx]
            y_va = y_full[va_idx]

            s_tr, _, _, _ = rmt_prior_from_train(x_tr)
            pb_tr = point_biserial_by_feature(x_tr, y_tr)
            s_norm = (s_tr - s_tr.min()) / (s_tr.max() - s_tr.min() + 1e-12)
            rank_score = np.abs(pb_tr) * np.power(s_norm, alpha)
            ranking = np.argsort(rank_score)[::-1]

            score_mat = np.cumsum(x_va[:, ranking[:max_m_eff]], axis=1)
            for m_idx in range(max_m_eff):
                score = score_mat[:, m_idx]
                if np.std(score) == 0:
                    continue
                lr = LogisticRegression(max_iter=10000, solver="lbfgs", random_state=seed)
                try:
                    # This reproduces Entry 44 exactly, including fitting the
                    # 1D diagnostic LR on the inner validation score vector.
                    lr.fit(score.reshape(-1, 1), y_va)
                    pred = lr.predict_proba(score.reshape(-1, 1))[:, 1]
                    pr = average_precision_score(y_va, pred)
                    pr_sum[m_idx] += pr
                    pr_sq_sum[m_idx] += pr * pr
                    valid_n[m_idx] += 1
                except Exception:
                    continue

    mean_pr = np.divide(pr_sum, valid_n, out=np.full(max_m_eff, np.nan), where=valid_n > 0)
    mean_sq = np.divide(pr_sq_sum, valid_n, out=np.full(max_m_eff, np.nan), where=valid_n > 0)
    std_pr = np.sqrt(np.clip(mean_sq - mean_pr * mean_pr, a_min=0.0, a_max=None))
    m_opt = 1 if np.all(np.isnan(mean_pr)) else int(np.nanargmax(mean_pr) + 1)

    scan = pd.DataFrame(
        {
            "m": np.arange(1, max_m_eff + 1),
            "mean_PR_AUC": mean_pr,
            "std_PR_AUC": std_pr,
            "n_valid": valid_n,
        }
    )
    return m_opt, scan


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA_DIR / "data.csv")
    labelled = data[data["activity"].notna()].copy()
    labelled["ligand_id"] = labelled["ligand_id"].astype(int)
    labelled["activity"] = labelled["activity"].astype(int)

    all_diag: dict[str, dict[str, float | int | str]] = {}
    selected_rows: list[dict[str, float | int | str]] = []

    for protein in PROTEINS:
        gene = protein["gene"]
        t0 = time.time()
        matrix = load_residue_matrix(RESIDUE_DIR, protein["key"])
        matrix["ligand_id"] = normalize_ligand_id(matrix["ligand_id"]).astype(int)
        merged = labelled[["ligand_id", "activity"]].merge(matrix, on="ligand_id", how="inner", validate="one_to_one")
        residue_cols = [c for c in merged.columns if c not in ("ligand_id", "activity")]
        y = merged["activity"].to_numpy(dtype=float)
        x = merged[residue_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)

        s_i, lambda_plus, n_signal, q = rmt_prior_from_train(x)
        r_pb = point_biserial_by_feature(x, y)
        s_norm = (s_i - s_i.min()) / (s_i.max() - s_i.min() + 1e-12)
        rank_score = np.abs(r_pb) * np.power(s_norm, ALPHA)
        ranking = np.argsort(rank_score)[::-1]

        ranked = pd.DataFrame(
            {
                "gene": gene,
                "residue": [residue_cols[i] for i in range(len(residue_cols))],
                "feature_name": [f"{gene}_{residue_cols[i]}" for i in range(len(residue_cols))],
                "r_pb": r_pb,
                "abs_r_pb": np.abs(r_pb),
                "s_score": s_i,
                "rank_score": rank_score,
                "rank": 0,
            }
        )
        ranked.loc[ranking, "rank"] = np.arange(1, len(ranking) + 1)
        ranked = ranked.sort_values("rank").reset_index(drop=True)
        ranked["sign"] = ranked["r_pb"].apply(lambda v: "active > inactive" if v < 0 else "inactive > active")

        m_opt, scan = select_m_inner_cv_pr_auc(x, y)
        best_pr = float(scan.loc[scan["m"] == m_opt, "mean_PR_AUC"].iloc[0])

        ranked.to_csv(OUT_DIR / f"{gene}_ranked_residues.csv", index=False)
        scan.to_csv(OUT_DIR / f"{gene}_m_scan.csv", index=False)

        top = ranked.head(m_opt).copy()
        top["selected_m_opt"] = m_opt
        selected_rows.extend(top.to_dict(orient="records"))

        all_diag[gene] = {
            "gene": gene,
            "n_residues": int(len(residue_cols)),
            "n_labelled": int(len(labelled)),
            "lambda_plus": float(lambda_plus),
            "n_signal": int(n_signal),
            "q": float(q),
            "m_opt": int(m_opt),
            "m_opt_expected_from_entry44": int(protein["top_m_expected"]),
            "best_PR_AUC": best_pr,
            "top1": str(ranked.iloc[0]["residue"]),
            "time_seconds": round(time.time() - t0, 2),
        }

    selected = pd.DataFrame(selected_rows)
    selected.to_csv(OUT_DIR / "selected_rmt95_features.csv", index=False)
    (OUT_DIR / "diag.json").write_text(json.dumps(all_diag, indent=2), encoding="utf-8")

    print("Saved RMT selection to", OUT_DIR)
    print(selected.groupby("gene").size().to_string())


if __name__ == "__main__":
    main()
