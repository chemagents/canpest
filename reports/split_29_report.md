# Split split_29 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.37% | 34 | 0.9321 | 0.9469 | 0.8512 |
| CatBoost RMT95 standalone | 31.77% | 95 | 0.8012 | 0.8002 | 0.7619 |
| SVM RMT95 standalone | 28.76% | 86 | 0.7863 | 0.7771 | 0.7530 |
| DMPNN x CatBoost RMT95 | 5.02% | 15 | 0.9141 | 0.9221 | 0.6667 |

False positives removed by veto: 19/34 (55.9%). True positives lost: 62/286 (21.7%).
