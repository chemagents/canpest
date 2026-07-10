#!/usr/bin/env python3
"""Run one-split DMPNN x CatBoost(RMT95) veto reproduction.

Default is split_00. RTE390 is intentionally not used in this rerun.

Run from the repository root:
    python3 run_split_veto.py --split-id split_00
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier
from sklearn.pipeline import make_pipeline
from sklearn.metrics import average_precision_score, f1_score, precision_score, recall_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUT_DIR = SCRIPT_DIR
SELECTION_CSV = OUT_DIR / "rmt_selection" / "selected_rmt95_features.csv"
RMT_FEATURE_DIR = OUT_DIR / "data" / "rmt95_features"
MODEL_DIR = OUT_DIR / "models"
PRED_DIR = OUT_DIR / "predictions"
METRIC_DIR = OUT_DIR / "metrics"
REPORT_DIR = OUT_DIR / "reports"

DATA_CSV = OUT_DIR / "data" / "data.csv"
if not DATA_CSV.exists():
    DATA_CSV = PROJECT_ROOT / "data" / "data.csv"
SPLIT_CSV = OUT_DIR / "data" / "splits" / "split_registry.csv"
if not SPLIT_CSV.exists():
    SPLIT_CSV = PROJECT_ROOT / "pipeline_runs" / "staged_pipeline" / "split_registry.csv"
DMPNN_DIR = MODEL_DIR / "dmpnn_rdkit2d"
if not DMPNN_DIR.exists():
    DMPNN_DIR = PROJECT_ROOT / "pipeline_runs" / "dmpnn_rdkit2d" / "rdkit2d"
RESIDUE_DIR = OUT_DIR / "source_inputs" / "residue_matrices"
if not RESIDUE_DIR.exists():
    RESIDUE_DIR = PROJECT_ROOT / "data" / "residue_matrices"

SOURCE_CODE_DIR = OUT_DIR / "source_code"
if (SOURCE_CODE_DIR / "rmt_filter.py").exists():
    sys.path.insert(0, str(SOURCE_CODE_DIR))
else:
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from rmt_filter import load_residue_matrix, normalize_ligand_id  # noqa: E402


NON_FEATURE = {"split_id", "set", "ligand_id", "activity"}
PROTEIN_KEYS = {"GLC1": "glc1", "ACHE": "D8V7J0"}
EPS = 1e-8


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reproduce DMPNN x CatBoost(RMT95) veto for one split")
    p.add_argument("--split-id", default="split_00")
    p.add_argument("--threshold", type=float, default=0.5)
    return p.parse_args()


def ensure_dirs() -> None:
    for path in [RMT_FEATURE_DIR, MODEL_DIR, PRED_DIR, METRIC_DIR, REPORT_DIR, OUT_DIR / "data" / "splits", OUT_DIR / "data" / "dmpnn_predictions"]:
        path.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        if src.resolve() == dst.resolve():
            return
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def display_path(path: Path) -> str:
    """Return a path relative to the repository root when possible."""
    try:
        return str(path.resolve().relative_to(OUT_DIR.resolve()))
    except ValueError:
        return str(path)


def catboost_model(seed: int) -> CatBoostClassifier:
    return CatBoostClassifier(
        iterations=500,
        learning_rate=0.05,
        depth=6,
        loss_function="Logloss",
        eval_metric="AUC",
        early_stopping_rounds=30,
        random_seed=seed,
        verbose=0,
        thread_count=1,
    )


def svm_model(seed: int):
    return make_pipeline(
        StandardScaler(),
        SVC(
            C=1.0,
            kernel="rbf",
            gamma="scale",
            class_weight="balanced",
            probability=True,
            random_state=seed,
        ),
    )


def confusion_counts(y: np.ndarray, p: np.ndarray, threshold: float) -> dict[str, int]:
    pred = p >= threshold
    return {
        "tp": int((pred & (y == 1)).sum()),
        "fp": int((pred & (y == 0)).sum()),
        "tn": int(((~pred) & (y == 0)).sum()),
        "fn": int(((~pred) & (y == 1)).sum()),
        "n_pred_pos": int(pred.sum()),
    }


def metrics_row(split_id: str, method: str, y: np.ndarray, p: np.ndarray, dmpnn: np.ndarray, threshold: float) -> dict[str, float | int | str]:
    c = confusion_counts(y, p, threshold)
    b = confusion_counts(y, dmpnn, threshold)
    pred = p >= threshold
    base_pred = dmpnn >= threshold
    fp_removed = int((base_pred & (y == 0) & ~pred).sum())
    tp_lost = int((base_pred & (y == 1) & ~pred).sum())
    return {
        "split_id": split_id,
        "method": method,
        "threshold": threshold,
        "n_test": int(len(y)),
        "n_pos_true": int((y == 1).sum()),
        "n_neg_true": int((y == 0).sum()),
        "roc_auc": float(roc_auc_score(y, p)),
        "pr_auc": float(average_precision_score(y, p)),
        "fpr": c["fp"] / max(c["fp"] + c["tn"], 1),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "recall": float(recall_score(y, pred, zero_division=0)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        **c,
        "dmpnn_tp": b["tp"],
        "dmpnn_fp": b["fp"],
        "dmpnn_tn": b["tn"],
        "dmpnn_fn": b["fn"],
        "fp_removed": fp_removed,
        "tp_lost": tp_lost,
        "fp_removed_frac": fp_removed / max(b["fp"], 1),
        "tp_lost_frac": tp_lost / max(b["tp"], 1),
    }


def load_selected_feature_frames() -> pd.DataFrame:
    selected = pd.read_csv(SELECTION_CSV)
    all_frames = []
    for gene, key in PROTEIN_KEYS.items():
        residues = selected[selected["gene"] == gene]["residue"].tolist()
        matrix = load_residue_matrix(RESIDUE_DIR, key)
        matrix["ligand_id"] = normalize_ligand_id(matrix["ligand_id"]).astype(int)
        missing = [r for r in residues if r not in matrix.columns]
        if missing:
            raise ValueError(f"{gene}: {len(missing)} selected residues missing from residue matrix: {missing[:5]}")
        frame = matrix[["ligand_id"] + residues].copy()
        frame.columns = ["ligand_id"] + [f"{gene}_{r}" for r in residues]
        all_frames.append(frame)

    merged = all_frames[0]
    for frame in all_frames[1:]:
        merged = merged.merge(frame, on="ligand_id", how="outer")
    return merged.fillna(0.0)


def write_split_feature_csv(split_id: str, features: pd.DataFrame, reg: pd.DataFrame, data: pd.DataFrame) -> tuple[Path, Path]:
    feature_cols = [c for c in features.columns if c != "ligand_id"]
    reg_s = reg[reg["split_id"] == split_id].copy()
    split = reg_s[["split_id", "set", "ligand_id", "activity"]].merge(features, on="ligand_id", how="left")
    split[feature_cols] = split[feature_cols].fillna(0.0)

    split_dir = RMT_FEATURE_DIR / split_id
    split_dir.mkdir(parents=True, exist_ok=True)
    train_csv = split_dir / "rmt95_train.csv"
    test_csv = split_dir / "rmt95_test.csv"
    split[split["set"] == "train"].to_csv(train_csv, index=False)
    split[split["set"] == "test"].to_csv(test_csv, index=False)

    selected_features = pd.DataFrame({"feature_name": feature_cols})
    selected_features.to_csv(split_dir / "feature_names.csv", index=False)
    return train_csv, test_csv


def main() -> None:
    args = parse_args()
    ensure_dirs()
    split_id = args.split_id
    split_num = int(split_id.split("_")[1])
    threshold = args.threshold

    copy_if_exists(DATA_CSV, OUT_DIR / "data" / "data.csv")
    copy_if_exists(SPLIT_CSV, OUT_DIR / "data" / "splits" / "split_registry.csv")
    copy_if_exists(PROJECT_ROOT / "pipeline_runs" / "staged_pipeline" / "split_summary.csv", OUT_DIR / "data" / "splits" / "split_summary.csv")
    copy_if_exists(PROJECT_ROOT / "pipeline_runs" / "staged_pipeline" / "metadata.json", OUT_DIR / "data" / "splits" / "metadata.json")

    data = pd.read_csv(DATA_CSV)
    data["ligand_id"] = data["ligand_id"].astype(int)
    reg = pd.read_csv(SPLIT_CSV)
    reg["ligand_id"] = reg["ligand_id"].astype(int)

    features = load_selected_feature_frames()
    train_csv, test_csv = write_split_feature_csv(split_id, features, reg, data)

    split_out = MODEL_DIR / "dmpnn_rdkit2d" / split_id
    for filename in [
        "dmpnn_model.pt",
        "hgb_model.pkl",
        "feature_scaler.pkl",
        "global_feature_names.npy",
        "global_feature_names_after_filter.npy",
        "selected_rmt_feature_names.npy",
        "metrics.json",
        "val_predictions.csv",
        "test_predictions.csv",
    ]:
        copy_if_exists(DMPNN_DIR / split_id / filename, split_out / filename)

    dmpnn_val = pd.read_csv(DMPNN_DIR / split_id / "val_predictions.csv")
    dmpnn_test = pd.read_csv(DMPNN_DIR / split_id / "test_predictions.csv")
    dmpnn_val["ligand_id"] = dmpnn_val["ligand_id"].astype(int)
    dmpnn_test["ligand_id"] = dmpnn_test["ligand_id"].astype(int)
    copy_if_exists(DMPNN_DIR / split_id / "val_predictions.csv", OUT_DIR / "data" / "dmpnn_predictions" / f"{split_id}_val_predictions.csv")
    copy_if_exists(DMPNN_DIR / split_id / "test_predictions.csv", OUT_DIR / "data" / "dmpnn_predictions" / f"{split_id}_test_predictions.csv")

    train_df = pd.read_csv(train_csv)
    test_df = pd.read_csv(test_csv)
    feature_cols = [c for c in train_df.columns if c not in NON_FEATURE]

    val_ids = dmpnn_val["ligand_id"].to_numpy(dtype=int)
    test_ids = dmpnn_test["ligand_id"].to_numpy(dtype=int)
    train_ids_all = train_df["ligand_id"].to_numpy(dtype=int)
    train_ids = np.array([x for x in train_ids_all if x not in set(val_ids)], dtype=int)

    train_index = train_df.set_index("ligand_id")
    test_index = test_df.set_index("ligand_id")
    x_train = train_index.loc[train_ids, feature_cols].to_numpy(dtype=np.float64)
    y_train = train_index.loc[train_ids, "activity"].to_numpy(dtype=int)
    x_val = train_index.loc[val_ids, feature_cols].to_numpy(dtype=np.float64)
    y_val = train_index.loc[val_ids, "activity"].to_numpy(dtype=int)
    x_test = test_index.loc[test_ids, feature_cols].to_numpy(dtype=np.float64)
    y_test = dmpnn_test.set_index("ligand_id").loc[test_ids, "activity"].to_numpy(dtype=int)
    p_dmpnn = dmpnn_test.set_index("ligand_id").loc[test_ids, "blend_proba"].to_numpy(dtype=float)

    model = catboost_model(seed=42 + split_num)
    model.fit(x_train, y_train, eval_set=(x_val, y_val), use_best_model=True)
    p_cb_val = model.predict_proba(x_val)[:, 1]
    p_cb_test = model.predict_proba(x_test)[:, 1]

    svm = svm_model(seed=42 + split_num)
    svm.fit(x_train, y_train)
    p_svm_val = svm.predict_proba(x_val)[:, 1]
    p_svm_test = svm.predict_proba(x_test)[:, 1]
    p_combo = p_dmpnn * np.clip(p_cb_test, EPS, 1.0)

    cb_dir = MODEL_DIR / "catboost_rmt95" / split_id
    cb_dir.mkdir(parents=True, exist_ok=True)
    model.save_model(cb_dir / "model.cbm")
    pd.DataFrame({"feature_name": feature_cols}).to_csv(cb_dir / "feature_names.csv", index=False)

    svm_dir = MODEL_DIR / "svm_rmt95" / split_id
    svm_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(svm, svm_dir / "model.joblib")
    pd.DataFrame({"feature_name": feature_cols}).to_csv(svm_dir / "feature_names.csv", index=False)

    val_pred = pd.DataFrame({"ligand_id": val_ids, "activity": y_val, "p_catboost_rmt95": p_cb_val})
    val_pred.to_csv(PRED_DIR / f"{split_id}_catboost_rmt95_val.csv", index=False)
    svm_val_pred = pd.DataFrame({"ligand_id": val_ids, "activity": y_val, "p_svm_rmt95": p_svm_val})
    svm_val_pred.to_csv(PRED_DIR / f"{split_id}_svm_rmt95_val.csv", index=False)

    pred = pd.DataFrame(
        {
            "split_id": split_id,
            "ligand_id": test_ids,
            "activity": y_test,
            "p_dmpnn": p_dmpnn,
            "p_catboost_rmt95": p_cb_test,
            "p_svm_rmt95": p_svm_test,
            "p_final_product_rmt95": p_combo,
        }
    )
    pred["dmpnn_pred"] = (pred["p_dmpnn"] >= threshold).astype(int)
    pred["veto_pred"] = (pred["p_final_product_rmt95"] >= threshold).astype(int)
    pred["fp_before"] = ((pred["activity"] == 0) & (pred["dmpnn_pred"] == 1)).astype(int)
    pred["fp_after"] = ((pred["activity"] == 0) & (pred["veto_pred"] == 1)).astype(int)
    pred["tp_lost_by_veto"] = ((pred["activity"] == 1) & (pred["dmpnn_pred"] == 1) & (pred["veto_pred"] == 0)).astype(int)
    pred["vetoed"] = ((pred["dmpnn_pred"] == 1) & (pred["veto_pred"] == 0)).astype(int)
    pred.to_csv(PRED_DIR / f"{split_id}_test_predictions_with_veto.csv", index=False)

    rows = [
        metrics_row(split_id, "dmpnn", y_test, p_dmpnn, p_dmpnn, threshold),
        metrics_row(split_id, "catboost_rmt95_standalone", y_test, p_cb_test, p_dmpnn, threshold),
        metrics_row(split_id, "svm_rmt95_standalone", y_test, p_svm_test, p_dmpnn, threshold),
        metrics_row(split_id, "product_rmt95_veto", y_test, p_combo, p_dmpnn, threshold),
    ]
    metrics = pd.DataFrame(rows)
    metrics.to_csv(METRIC_DIR / f"{split_id}_metrics.csv", index=False)

    baseline = metrics[metrics["method"] == "dmpnn"].iloc[0]
    veto = metrics[metrics["method"] == "product_rmt95_veto"].iloc[0]
    support = metrics[metrics["method"] == "catboost_rmt95_standalone"].iloc[0]
    svm_support = metrics[metrics["method"] == "svm_rmt95_standalone"].iloc[0]
    summary = {
        "split_id": split_id,
        "threshold": threshold,
        "n_features_rmt95": len(feature_cols),
        "catboost": {
            "iterations": 500,
            "learning_rate": 0.05,
            "depth": 6,
            "early_stopping_rounds": 30,
            "random_seed": 42 + split_num,
            "best_iteration": int(model.get_best_iteration() or 0),
        },
        "validation_protocol": "DMPNN val_predictions ligand_ids excluded from CatBoost train and used as eval_set.",
        "dmpnn_fpr": float(baseline["fpr"]),
        "veto_fpr": float(veto["fpr"]),
        "fpr_absolute_reduction_pp": float((baseline["fpr"] - veto["fpr"]) * 100),
        "fpr_relative_reduction_frac": float((baseline["fpr"] - veto["fpr"]) / max(float(baseline["fpr"]), EPS)),
        "dmpnn_fp": int(baseline["fp"]),
        "veto_fp": int(veto["fp"]),
        "fp_removed": int(veto["fp_removed"]),
        "tp_lost": int(veto["tp_lost"]),
        "dmpnn_roc_auc": float(baseline["roc_auc"]),
        "veto_roc_auc": float(veto["roc_auc"]),
        "catboost_rmt95_roc_auc": float(support["roc_auc"]),
        "svm_rmt95_roc_auc": float(svm_support["roc_auc"]),
        "outputs": {
            "catboost_model": display_path(cb_dir / "model.cbm"),
            "svm_model": display_path(svm_dir / "model.joblib"),
            "predictions": display_path(PRED_DIR / f"{split_id}_test_predictions_with_veto.csv"),
            "metrics": display_path(METRIC_DIR / f"{split_id}_metrics.csv"),
        },
    }
    (METRIC_DIR / f"{split_id}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    manifest_rows = []
    manifest_paths = [
        OUT_DIR / "data" / "data.csv",
        OUT_DIR / "data" / "splits" / "split_registry.csv",
        OUT_DIR / "data" / "splits" / "split_summary.csv",
        OUT_DIR / "data" / "splits" / "metadata.json",
        SELECTION_CSV,
        train_csv,
        test_csv,
        OUT_DIR / "data" / "dmpnn_predictions" / f"{split_id}_val_predictions.csv",
        OUT_DIR / "data" / "dmpnn_predictions" / f"{split_id}_test_predictions.csv",
        split_out / "dmpnn_model.pt",
        split_out / "feature_scaler.pkl",
        split_out / "global_feature_names_after_filter.npy",
        cb_dir / "model.cbm",
        svm_dir / "model.joblib",
        PRED_DIR / f"{split_id}_test_predictions_with_veto.csv",
        METRIC_DIR / f"{split_id}_metrics.csv",
    ]
    for path in manifest_paths:
        manifest_rows.append({"path": str(path.relative_to(OUT_DIR)), "sha256": sha256(path), "bytes": path.stat().st_size})
    pd.DataFrame(manifest_rows).to_csv(OUT_DIR / "sha256_manifest.csv", index=False)

    report = f"""# Split {split_id} RMT95 Veto Reproduction\n\n""" \
        f"Dataset split: `{display_path(SPLIT_CSV)}`. DMPNN baseline: `{display_path(DMPNN_DIR / split_id)}`.\n\n" \
        f"RMT95 features: {len(feature_cols)} residues (45 GLC1 + 50 ACHE). RTE390 was not used.\n\n" \
        f"Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.\n\n" \
        f"| Method | FPR | FP | ROC AUC | PR AUC | Recall |\n" \
        f"|---|---:|---:|---:|---:|---:|\n" \
        f"| DMPNN | {baseline['fpr']*100:.2f}% | {int(baseline['fp'])} | {baseline['roc_auc']:.4f} | {baseline['pr_auc']:.4f} | {baseline['recall']:.4f} |\n" \
        f"| CatBoost RMT95 standalone | {support['fpr']*100:.2f}% | {int(support['fp'])} | {support['roc_auc']:.4f} | {support['pr_auc']:.4f} | {support['recall']:.4f} |\n" \
        f"| SVM RMT95 standalone | {svm_support['fpr']*100:.2f}% | {int(svm_support['fp'])} | {svm_support['roc_auc']:.4f} | {svm_support['pr_auc']:.4f} | {svm_support['recall']:.4f} |\n" \
        f"| DMPNN x CatBoost RMT95 | {veto['fpr']*100:.2f}% | {int(veto['fp'])} | {veto['roc_auc']:.4f} | {veto['pr_auc']:.4f} | {veto['recall']:.4f} |\n\n" \
        f"False positives removed by veto: {int(veto['fp_removed'])}/{int(baseline['fp'])} " \
        f"({veto['fp_removed_frac']*100:.1f}%). True positives lost: {int(veto['tp_lost'])}/{int(baseline['tp'])} " \
        f"({veto['tp_lost_frac']*100:.1f}%).\n\n" \
        f"FPR reduction: {summary['fpr_absolute_reduction_pp']:.2f} percentage points, " \
        f"{summary['fpr_relative_reduction_frac']*100:.1f}% relative to DMPNN on this split.\n"
    (REPORT_DIR / f"{split_id}_report.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
