# Split split_08 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.71% | 35 | 0.9306 | 0.9455 | 0.8333 |
| CatBoost RMT95 standalone | 34.11% | 102 | 0.7725 | 0.7513 | 0.7679 |
| SVM RMT95 standalone | 34.78% | 104 | 0.7601 | 0.7388 | 0.8006 |
| DMPNN x CatBoost RMT95 | 4.68% | 14 | 0.9162 | 0.9249 | 0.6429 |

False positives removed by veto: 21/35 (60.0%). True positives lost: 64/280 (22.9%).
