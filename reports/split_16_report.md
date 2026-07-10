# Split split_16 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 15.38% | 46 | 0.9092 | 0.9257 | 0.8333 |
| CatBoost RMT95 standalone | 41.47% | 124 | 0.7287 | 0.7172 | 0.7887 |
| SVM RMT95 standalone | 37.46% | 112 | 0.7490 | 0.7320 | 0.7857 |
| DMPNN x CatBoost RMT95 | 8.03% | 24 | 0.8914 | 0.8957 | 0.6756 |

False positives removed by veto: 22/46 (47.8%). True positives lost: 53/280 (18.9%).
