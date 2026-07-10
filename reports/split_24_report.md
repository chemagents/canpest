# Split split_24 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.04% | 33 | 0.9327 | 0.9409 | 0.8452 |
| CatBoost RMT95 standalone | 36.79% | 110 | 0.7741 | 0.7872 | 0.7708 |
| SVM RMT95 standalone | 32.11% | 96 | 0.7809 | 0.7824 | 0.7619 |
| DMPNN x CatBoost RMT95 | 5.35% | 16 | 0.9174 | 0.9243 | 0.6696 |

False positives removed by veto: 17/33 (51.5%). True positives lost: 59/284 (20.8%).
