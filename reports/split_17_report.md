# Split split_17 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 9.70% | 29 | 0.9332 | 0.9470 | 0.8601 |
| CatBoost RMT95 standalone | 34.45% | 103 | 0.8159 | 0.8264 | 0.8185 |
| SVM RMT95 standalone | 32.44% | 97 | 0.7973 | 0.8018 | 0.8006 |
| DMPNN x CatBoost RMT95 | 3.68% | 11 | 0.9213 | 0.9356 | 0.7381 |

False positives removed by veto: 18/29 (62.1%). True positives lost: 41/289 (14.2%).
