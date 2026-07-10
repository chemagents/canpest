# Reproduction Process

## Commands

Run from the repository root.

```bash
python3 select_rmt_residues.py
python3 run_split_veto.py --split-id split_00
python3 run_30splits_veto.py
```

## What Was Reproduced

The current run reproduces only the RMT95 product-veto branch, not the later RTE390 geometric mix.

The reproduced classifier chain is:

```text
SMILES + RDKit2D -> DMPNN probability p_DMPNN
RMT95 residue energies -> CatBoost probability p_CatBoost_RMT95
RMT95 residue energies -> SVM probability p_SVM_RMT95 (diagnostic only)
p_final = p_DMPNN * p_CatBoost_RMT95
```

If CatBoost sees weak residue-level docking support, the product becomes smaller and the candidate may fall below the 0.5 pesticide threshold. This is the mechanistic veto.

## Why It Reduces FPR

Many DMPNN false positives sit near the decision boundary. Multiplying by the docking-support probability pushes many of those weakly supported predictions below 0.5.

For `split_00`:

| Quantity | DMPNN | With RMT95 veto |
|---|---:|---:|
| False positives | 33 | 14 |
| FPR | 11.04% | 4.68% |
| ROC AUC | 0.9366 | 0.9265 |

Thus the veto removed `19` false positives from `33`, a `57.6%` reduction on this split.

The historical all-split RMT95 product-veto result is DMPNN FPR `12.20%` to product-veto FPR `5.14%`, approximately `58%` relative FPR reduction. The later Entry 53 manuscript text used the slightly stronger RMT95+RTE390 geometric mix (`4.92%`), but RTE390 is intentionally excluded from this rerun.

## 30-Split Reproduction Result

| Method | FPR mean ± std | ROC AUC mean ± std | FP total | FP removed | TP lost |
|---|---:|---:|---:|---:|---:|
| DMPNN | 12.20% ± 1.86 | 0.9283 ± 0.0092 | 1094 | 0 | 0 |
| CatBoost(RMT95) standalone | 35.69% ± 3.71 | 0.7886 ± 0.0240 | 3201 | 331 | 1079 |
| SVM(RMT95) standalone | 33.48% ± 2.58 | 0.7851 ± 0.0164 | 3003 | 366 | 1192 |
| DMPNN x CatBoost(RMT95) | 5.14% ± 1.03 | 0.9167 ± 0.0088 | 461 | 633 | 1727 |

Relative FPR reduction: `57.9%`.

The new outputs match Entry 52 `product_rmt95_prev` exactly in counts and up to floating-point roundoff in probabilities-derived metrics.

## Completeness Boundary

This folder is sufficient to reproduce the reported RMT95 veto experiment:

- Rerun RMT selection from raw residue matrices.
- Rebuild the 95 selected residue features.
- Retrain CatBoost(RMT95) and SVM(RMT95) for each split.
- Combine saved DMPNN probabilities with CatBoost probabilities.
- Recompute FPR reduction across 30 splits.

The folder includes saved DMPNN models and predictions. It also includes the original DMPNN training script and RDKit2D descriptor cache for traceability, but it does not yet include a one-command wrapper that retrains all 30 DMPNN models from scratch. Therefore, if “full reproduction” is defined as retraining DMPNN weights from raw SMILES, one more wrapper should be added. If “full reproduction” is defined as reproducing the logged FPR-veto experiment and verifying the result, the folder is sufficient.
