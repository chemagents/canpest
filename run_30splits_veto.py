#!/usr/bin/env python3
"""Run DMPNN x CatBoost(RMT95) veto reproduction across split registry.

This is the all-split extension of run_split_veto.py. It intentionally does not
use RTE390. The support model is CatBoost trained only on the 95 RMT-selected
residue e_total features reproduced by select_rmt_residues.py.

Run from the repository root:
    python3 run_30splits_veto.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_split_veto as one


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run RMT95 product-veto over multiple split IDs")
    p.add_argument("--split-ids", default="all", help="all or comma-list, e.g. split_00,split_01")
    p.add_argument("--threshold", type=float, default=0.5)
    return p.parse_args()


def split_list(value: str, reg: pd.DataFrame) -> list[str]:
    all_ids = sorted(reg["split_id"].unique())
    if value == "all":
        return all_ids
    wanted = [x.strip() for x in value.split(",") if x.strip()]
    missing = sorted(set(wanted) - set(all_ids))
    if missing:
        raise ValueError(f"Unknown split IDs: {missing}")
    return wanted


def summarize(per_split: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, g in per_split.groupby("method", sort=False):
        rows.append(
            {
                "method": method,
                "n_splits": int(len(g)),
                "fpr_mean": float(g["fpr"].mean()),
                "fpr_std": float(g["fpr"].std()),
                "roc_auc_mean": float(g["roc_auc"].mean()),
                "roc_auc_std": float(g["roc_auc"].std()),
                "pr_auc_mean": float(g["pr_auc"].mean()),
                "pr_auc_std": float(g["pr_auc"].std()),
                "precision_mean": float(g["precision"].mean()),
                "precision_std": float(g["precision"].std()),
                "recall_mean": float(g["recall"].mean()),
                "recall_std": float(g["recall"].std()),
                "f1_mean": float(g["f1"].mean()),
                "f1_std": float(g["f1"].std()),
                "fp_total": int(g["fp"].sum()),
                "tp_total": int(g["tp"].sum()),
                "tn_total": int(g["tn"].sum()),
                "fn_total": int(g["fn"].sum()),
                "fp_removed_total": int(g["fp_removed"].sum()),
                "tp_lost_total": int(g["tp_lost"].sum()),
                "fp_removed_frac_total": float(g["fp_removed"].sum() / max(g["dmpnn_fp"].sum(), 1)),
                "tp_lost_frac_total": float(g["tp_lost"].sum() / max(g["dmpnn_tp"].sum(), 1)),
            }
        )
    return pd.DataFrame(rows)


def compare_to_entry52(per_split: pd.DataFrame) -> pd.DataFrame:
    hist_path = one.PROJECT_ROOT / "LOG" / "52_full_docking_mix_model" / "per_split_methods.csv"
    if not hist_path.exists():
        return pd.DataFrame()

    current = per_split[per_split["method"].isin(["dmpnn", "product_rmt95_veto"])].copy()
    current["entry52_method"] = current["method"].replace({"product_rmt95_veto": "product_rmt95_prev"})
    hist = pd.read_csv(hist_path)
    hist = hist[hist["method"].isin(["dmpnn", "product_rmt95_prev"])]

    merged = current.merge(
        hist,
        left_on=["split_id", "entry52_method"],
        right_on=["split_id", "method"],
        suffixes=("_re", "_entry52"),
    )
    rows = []
    for _, row in merged.iterrows():
        rows.append(
            {
                "split_id": row["split_id"],
                "method_re": row["method_re"],
                "method_entry52": row["method_entry52"],
                "fpr_diff": float(row["fpr_re"] - row["fpr_entry52"]),
                "roc_auc_diff": float(row["roc_auc_re"] - row["roc_auc_entry52"]),
                "pr_auc_diff": float(row["pr_auc_re"] - row["pr_auc_entry52"]),
                "fp_diff": int(row["fp_re"] - row["fp_entry52"]),
                "tp_diff": int(row["tp_re"] - row["tp_entry52"]),
                "fp_removed_diff": int(row["fp_removed_re"] - row["fp_removed_entry52"]),
                "tp_lost_diff": int(row["tp_lost_re"] - row["tp_lost_entry52"]),
            }
        )
    return pd.DataFrame(rows)


def process_split(
    split_id: str,
    features: pd.DataFrame,
    reg: pd.DataFrame,
    data: pd.DataFrame,
    threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    split_num = int(split_id.split("_")[1])
    train_csv, test_csv = one.write_split_feature_csv(split_id, features, reg, data)

    split_out = one.MODEL_DIR / "dmpnn_rdkit2d" / split_id
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
        one.copy_if_exists(one.DMPNN_DIR / split_id / filename, split_out / filename)

    dmpnn_val = pd.read_csv(one.DMPNN_DIR / split_id / "val_predictions.csv")
    dmpnn_test = pd.read_csv(one.DMPNN_DIR / split_id / "test_predictions.csv")
    dmpnn_val["ligand_id"] = dmpnn_val["ligand_id"].astype(int)
    dmpnn_test["ligand_id"] = dmpnn_test["ligand_id"].astype(int)
    one.copy_if_exists(one.DMPNN_DIR / split_id / "val_predictions.csv", one.OUT_DIR / "data" / "dmpnn_predictions" / f"{split_id}_val_predictions.csv")
    one.copy_if_exists(one.DMPNN_DIR / split_id / "test_predictions.csv", one.OUT_DIR / "data" / "dmpnn_predictions" / f"{split_id}_test_predictions.csv")

    train_df = pd.read_csv(train_csv)
    test_df = pd.read_csv(test_csv)
    feature_cols = [c for c in train_df.columns if c not in one.NON_FEATURE]

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

    model = one.catboost_model(seed=42 + split_num)
    model.fit(x_train, y_train, eval_set=(x_val, y_val), use_best_model=True)
    p_cb_val = model.predict_proba(x_val)[:, 1]
    p_cb_test = model.predict_proba(x_test)[:, 1]

    svm = one.svm_model(seed=42 + split_num)
    svm.fit(x_train, y_train)
    p_svm_val = svm.predict_proba(x_val)[:, 1]
    p_svm_test = svm.predict_proba(x_test)[:, 1]
    p_combo = p_dmpnn * np.clip(p_cb_test, one.EPS, 1.0)

    cb_dir = one.MODEL_DIR / "catboost_rmt95" / split_id
    cb_dir.mkdir(parents=True, exist_ok=True)
    model.save_model(cb_dir / "model.cbm")
    pd.DataFrame({"feature_name": feature_cols}).to_csv(cb_dir / "feature_names.csv", index=False)

    svm_dir = one.MODEL_DIR / "svm_rmt95" / split_id
    svm_dir.mkdir(parents=True, exist_ok=True)
    one.joblib.dump(svm, svm_dir / "model.joblib")
    pd.DataFrame({"feature_name": feature_cols}).to_csv(svm_dir / "feature_names.csv", index=False)

    val_pred = pd.DataFrame({"split_id": split_id, "ligand_id": val_ids, "activity": y_val, "p_catboost_rmt95": p_cb_val})
    val_pred.to_csv(one.PRED_DIR / f"{split_id}_catboost_rmt95_val.csv", index=False)
    svm_val_pred = pd.DataFrame({"split_id": split_id, "ligand_id": val_ids, "activity": y_val, "p_svm_rmt95": p_svm_val})
    svm_val_pred.to_csv(one.PRED_DIR / f"{split_id}_svm_rmt95_val.csv", index=False)

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
    pred.to_csv(one.PRED_DIR / f"{split_id}_test_predictions_with_veto.csv", index=False)

    rows = [
        one.metrics_row(split_id, "dmpnn", y_test, p_dmpnn, p_dmpnn, threshold),
        one.metrics_row(split_id, "catboost_rmt95_standalone", y_test, p_cb_test, p_dmpnn, threshold),
        one.metrics_row(split_id, "svm_rmt95_standalone", y_test, p_svm_test, p_dmpnn, threshold),
        one.metrics_row(split_id, "product_rmt95_veto", y_test, p_combo, p_dmpnn, threshold),
    ]
    metrics = pd.DataFrame(rows)
    metrics.to_csv(one.METRIC_DIR / f"{split_id}_metrics.csv", index=False)

    baseline = metrics[metrics["method"] == "dmpnn"].iloc[0]
    veto = metrics[metrics["method"] == "product_rmt95_veto"].iloc[0]
    support = metrics[metrics["method"] == "catboost_rmt95_standalone"].iloc[0]
    svm_support = metrics[metrics["method"] == "svm_rmt95_standalone"].iloc[0]
    split_summary = {
        "split_id": split_id,
        "threshold": threshold,
        "n_features_rmt95": len(feature_cols),
        "catboost_best_iteration": int(model.get_best_iteration() or 0),
        "dmpnn_fpr": float(baseline["fpr"]),
        "veto_fpr": float(veto["fpr"]),
        "fpr_absolute_reduction_pp": float((baseline["fpr"] - veto["fpr"]) * 100),
        "fpr_relative_reduction_frac": float((baseline["fpr"] - veto["fpr"]) / max(float(baseline["fpr"]), one.EPS)),
        "dmpnn_fp": int(baseline["fp"]),
        "veto_fp": int(veto["fp"]),
        "fp_removed": int(veto["fp_removed"]),
        "tp_lost": int(veto["tp_lost"]),
        "dmpnn_roc_auc": float(baseline["roc_auc"]),
        "veto_roc_auc": float(veto["roc_auc"]),
        "catboost_rmt95_roc_auc": float(support["roc_auc"]),
        "svm_rmt95_roc_auc": float(svm_support["roc_auc"]),
    }
    (one.METRIC_DIR / f"{split_id}_summary.json").write_text(json.dumps(split_summary, indent=2), encoding="utf-8")

    report = f"""# Split {split_id} RMT95 Veto Reproduction\n\n""" \
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
        f"({veto['tp_lost_frac']*100:.1f}%).\n"
    (one.REPORT_DIR / f"{split_id}_report.md").write_text(report, encoding="utf-8")

    print(
        f"{split_id}: DMPNN FPR={baseline['fpr']*100:.2f}% FP={int(baseline['fp'])}; "
        f"veto FPR={veto['fpr']*100:.2f}% FP={int(veto['fp'])}; "
        f"removed={int(veto['fp_removed'])}/{int(baseline['fp'])}; "
        f"CB_ROC={support['roc_auc']:.4f}; SVM_ROC={svm_support['roc_auc']:.4f}; "
        f"best_iter={int(model.get_best_iteration() or 0)}"
    )
    return metrics, pred, split_summary


def main() -> None:
    args = parse_args()
    one.ensure_dirs()
    one.copy_if_exists(one.DATA_CSV, one.OUT_DIR / "data" / "data.csv")
    one.copy_if_exists(one.SPLIT_CSV, one.OUT_DIR / "data" / "splits" / "split_registry.csv")
    one.copy_if_exists(one.PROJECT_ROOT / "pipeline_runs" / "staged_pipeline" / "split_summary.csv", one.OUT_DIR / "data" / "splits" / "split_summary.csv")
    one.copy_if_exists(one.PROJECT_ROOT / "pipeline_runs" / "staged_pipeline" / "metadata.json", one.OUT_DIR / "data" / "splits" / "metadata.json")

    data = pd.read_csv(one.DATA_CSV)
    data["ligand_id"] = data["ligand_id"].astype(int)
    reg = pd.read_csv(one.SPLIT_CSV)
    reg["ligand_id"] = reg["ligand_id"].astype(int)
    ids = split_list(args.split_ids, reg)
    features = one.load_selected_feature_frames()

    all_metrics = []
    all_predictions = []
    split_summaries = []
    for split_id in ids:
        metrics, pred, split_summary = process_split(split_id, features, reg, data, args.threshold)
        all_metrics.append(metrics)
        all_predictions.append(pred)
        split_summaries.append(split_summary)

    per_split = pd.concat(all_metrics, ignore_index=True)
    pooled = pd.concat(all_predictions, ignore_index=True)
    summary = summarize(per_split)
    comparison = compare_to_entry52(per_split)

    per_split.to_csv(one.METRIC_DIR / "per_split_methods_rmt95.csv", index=False)
    pooled.to_csv(one.PRED_DIR / "pooled_30splits_predictions_with_veto.csv", index=False)
    summary.to_csv(one.METRIC_DIR / "method_summary_30splits.csv", index=False)
    pd.DataFrame(split_summaries).to_csv(one.METRIC_DIR / "split_summary_30splits.csv", index=False)
    if not comparison.empty:
        comparison.to_csv(one.METRIC_DIR / "compare_to_entry52_product_rmt95.csv", index=False)

    dmpnn = summary[summary["method"] == "dmpnn"].iloc[0]
    veto = summary[summary["method"] == "product_rmt95_veto"].iloc[0]
    result = {
        "split_ids": ids,
        "n_splits": len(ids),
        "threshold": args.threshold,
        "n_features_rmt95": int(len([c for c in features.columns if c != "ligand_id"])),
        "formula": "p_final = p_DMPNN * p_CatBoost_RMT95",
        "rte390_used": False,
        "method_summary": summary.to_dict(orient="records"),
        "headline": {
            "dmpnn_fpr_mean": float(dmpnn["fpr_mean"]),
            "veto_fpr_mean": float(veto["fpr_mean"]),
            "fpr_absolute_reduction_pp": float((dmpnn["fpr_mean"] - veto["fpr_mean"]) * 100),
            "fpr_relative_reduction_frac": float((dmpnn["fpr_mean"] - veto["fpr_mean"]) / max(float(dmpnn["fpr_mean"]), one.EPS)),
            "fp_removed_total": int(veto["fp_removed_total"]),
            "dmpnn_fp_total": int(dmpnn["fp_total"]),
            "tp_lost_total": int(veto["tp_lost_total"]),
            "dmpnn_tp_total": int(dmpnn["tp_total"]),
        },
    }
    if not comparison.empty:
        result["entry52_comparison"] = {
            "rows": int(len(comparison)),
            "max_abs_fpr_diff": float(comparison["fpr_diff"].abs().max()),
            "max_abs_roc_auc_diff": float(comparison["roc_auc_diff"].abs().max()),
            "max_abs_pr_auc_diff": float(comparison["pr_auc_diff"].abs().max()),
            "max_abs_fp_diff": int(comparison["fp_diff"].abs().max()),
            "max_abs_tp_diff": int(comparison["tp_diff"].abs().max()),
        }
    (one.METRIC_DIR / "summary_30splits.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    manifest_rows = []
    manifest_paths = [
        one.OUT_DIR / "data" / "data.csv",
        one.OUT_DIR / "data" / "splits" / "split_registry.csv",
        one.SELECTION_CSV,
        one.METRIC_DIR / "per_split_methods_rmt95.csv",
        one.METRIC_DIR / "method_summary_30splits.csv",
        one.METRIC_DIR / "summary_30splits.json",
        one.PRED_DIR / "pooled_30splits_predictions_with_veto.csv",
    ]
    for split_id in ids:
        manifest_paths.extend(
            [
                one.MODEL_DIR / "catboost_rmt95" / split_id / "model.cbm",
                one.MODEL_DIR / "svm_rmt95" / split_id / "model.joblib",
                one.MODEL_DIR / "dmpnn_rdkit2d" / split_id / "dmpnn_model.pt",
                one.OUT_DIR / "data" / "rmt95_features" / split_id / "rmt95_train.csv",
                one.OUT_DIR / "data" / "rmt95_features" / split_id / "rmt95_test.csv",
            ]
        )
    for path in manifest_paths:
        if path.exists():
            manifest_rows.append({"path": str(path.relative_to(one.OUT_DIR)), "sha256": one.sha256(path), "bytes": path.stat().st_size})
    pd.DataFrame(manifest_rows).to_csv(one.OUT_DIR / "sha256_manifest_30splits.csv", index=False)

    catboost = summary[summary["method"] == "catboost_rmt95_standalone"].iloc[0]
    svm = summary[summary["method"] == "svm_rmt95_standalone"].iloc[0]
    report = f"""# 30-Split DMPNN x CatBoost(RMT95) Veto Reproduction\n\n""" \
        f"Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`. RTE390 was not used.\n\n" \
        f"| Method | FPR mean ± std | ROC AUC mean ± std | PR AUC mean | Recall mean | FP total | FP removed vs DMPNN | TP lost vs DMPNN |\n" \
        f"|---|---:|---:|---:|---:|---:|---:|---:|\n" \
        f"| DMPNN | {dmpnn['fpr_mean']*100:.2f}% ± {dmpnn['fpr_std']*100:.2f} | {dmpnn['roc_auc_mean']:.4f} ± {dmpnn['roc_auc_std']:.4f} | {dmpnn['pr_auc_mean']:.4f} | {dmpnn['recall_mean']:.4f} | {int(dmpnn['fp_total'])} | 0 | 0 |\n" \
        f"| CatBoost(RMT95) standalone | {catboost['fpr_mean']*100:.2f}% ± {catboost['fpr_std']*100:.2f} | {catboost['roc_auc_mean']:.4f} ± {catboost['roc_auc_std']:.4f} | {catboost['pr_auc_mean']:.4f} | {catboost['recall_mean']:.4f} | {int(catboost['fp_total'])} | {int(catboost['fp_removed_total'])} ({catboost['fp_removed_frac_total']*100:.1f}%) | {int(catboost['tp_lost_total'])} ({catboost['tp_lost_frac_total']*100:.1f}%) |\n" \
        f"| SVM(RMT95) standalone | {svm['fpr_mean']*100:.2f}% ± {svm['fpr_std']*100:.2f} | {svm['roc_auc_mean']:.4f} ± {svm['roc_auc_std']:.4f} | {svm['pr_auc_mean']:.4f} | {svm['recall_mean']:.4f} | {int(svm['fp_total'])} | {int(svm['fp_removed_total'])} ({svm['fp_removed_frac_total']*100:.1f}%) | {int(svm['tp_lost_total'])} ({svm['tp_lost_frac_total']*100:.1f}%) |\n" \
        f"| DMPNN x CatBoost(RMT95) | {veto['fpr_mean']*100:.2f}% ± {veto['fpr_std']*100:.2f} | {veto['roc_auc_mean']:.4f} ± {veto['roc_auc_std']:.4f} | {veto['pr_auc_mean']:.4f} | {veto['recall_mean']:.4f} | {int(veto['fp_total'])} | {int(veto['fp_removed_total'])} ({veto['fp_removed_frac_total']*100:.1f}%) | {int(veto['tp_lost_total'])} ({veto['tp_lost_frac_total']*100:.1f}%) |\n\n" \
        f"Mean FPR reduction: {(dmpnn['fpr_mean'] - veto['fpr_mean'])*100:.2f} percentage points, " \
        f"{((dmpnn['fpr_mean'] - veto['fpr_mean']) / max(dmpnn['fpr_mean'], one.EPS))*100:.1f}% relative to DMPNN.\n\n"
    if not comparison.empty:
        report += (
            "Entry 52 comparison: max absolute differences are "
            f"FPR={comparison['fpr_diff'].abs().max():.3g}, "
            f"ROC={comparison['roc_auc_diff'].abs().max():.3g}, "
            f"PR={comparison['pr_auc_diff'].abs().max():.3g}, "
            f"FP={int(comparison['fp_diff'].abs().max())}, "
            f"TP={int(comparison['tp_diff'].abs().max())}.\n"
        )
    (one.REPORT_DIR / "30splits_report.md").write_text(report, encoding="utf-8")

    print("\n=== 30-split summary ===")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nSaved to {one.OUT_DIR}")


if __name__ == "__main__":
    main()
