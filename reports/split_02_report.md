# Split split_02 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 15.38% | 46 | 0.9333 | 0.9469 | 0.8810 |
| CatBoost RMT95 standalone | 35.45% | 106 | 0.7959 | 0.8086 | 0.7946 |
| SVM RMT95 standalone | 31.10% | 93 | 0.7929 | 0.8061 | 0.7857 |
| DMPNN x CatBoost RMT95 | 5.02% | 15 | 0.9176 | 0.9307 | 0.6429 |

False positives removed by veto: 31/46 (67.4%). True positives lost: 80/296 (27.0%).
