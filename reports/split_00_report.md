# Split split_00 RMT95 Veto Reproduction

Dataset split: `data/splits/split_registry.csv`. DMPNN baseline: `models/dmpnn_rdkit2d/split_00`.

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.04% | 33 | 0.9366 | 0.9488 | 0.8631 |
| CatBoost RMT95 standalone | 35.79% | 107 | 0.7815 | 0.7844 | 0.7798 |
| SVM RMT95 standalone | 31.77% | 95 | 0.7876 | 0.7797 | 0.7619 |
| DMPNN x CatBoost RMT95 | 4.68% | 14 | 0.9265 | 0.9265 | 0.6964 |

False positives removed by veto: 19/33 (57.6%). True positives lost: 56/290 (19.3%).

FPR reduction: 6.35 percentage points, 57.6% relative to DMPNN on this split.
