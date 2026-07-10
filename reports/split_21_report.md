# Split split_21 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.04% | 33 | 0.9328 | 0.9464 | 0.8601 |
| CatBoost RMT95 standalone | 32.78% | 98 | 0.7942 | 0.7934 | 0.8095 |
| SVM RMT95 standalone | 33.78% | 101 | 0.7848 | 0.7739 | 0.7917 |
| DMPNN x CatBoost RMT95 | 5.35% | 16 | 0.9203 | 0.9277 | 0.6845 |

False positives removed by veto: 17/33 (51.5%). True positives lost: 59/289 (20.4%).
