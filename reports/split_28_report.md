# Split split_28 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 15.05% | 45 | 0.9336 | 0.9468 | 0.8512 |
| CatBoost RMT95 standalone | 33.11% | 99 | 0.7922 | 0.7994 | 0.7679 |
| SVM RMT95 standalone | 31.44% | 94 | 0.7879 | 0.7783 | 0.7976 |
| DMPNN x CatBoost RMT95 | 6.69% | 20 | 0.9178 | 0.9187 | 0.6845 |

False positives removed by veto: 25/45 (55.6%). True positives lost: 56/286 (19.6%).
