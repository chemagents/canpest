# Split split_11 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 10.37% | 31 | 0.9365 | 0.9489 | 0.8214 |
| CatBoost RMT95 standalone | 38.46% | 115 | 0.7937 | 0.8064 | 0.8214 |
| SVM RMT95 standalone | 38.46% | 115 | 0.7892 | 0.7911 | 0.7917 |
| DMPNN x CatBoost RMT95 | 6.02% | 18 | 0.9239 | 0.9363 | 0.6786 |

False positives removed by veto: 13/31 (41.9%). True positives lost: 48/276 (17.4%).
