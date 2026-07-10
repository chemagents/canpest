# Split split_10 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.71% | 35 | 0.9210 | 0.9405 | 0.8423 |
| CatBoost RMT95 standalone | 29.77% | 89 | 0.8001 | 0.8105 | 0.7679 |
| SVM RMT95 standalone | 33.44% | 100 | 0.7813 | 0.7864 | 0.7679 |
| DMPNN x CatBoost RMT95 | 4.35% | 13 | 0.9105 | 0.9314 | 0.6786 |

False positives removed by veto: 22/35 (62.9%). True positives lost: 55/283 (19.4%).
