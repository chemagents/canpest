# Split split_20 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 12.04% | 36 | 0.9451 | 0.9542 | 0.8810 |
| CatBoost RMT95 standalone | 33.11% | 99 | 0.8014 | 0.8043 | 0.8185 |
| SVM RMT95 standalone | 32.44% | 97 | 0.7883 | 0.7740 | 0.8036 |
| DMPNN x CatBoost RMT95 | 5.35% | 16 | 0.9289 | 0.9405 | 0.6786 |

False positives removed by veto: 20/36 (55.6%). True positives lost: 68/296 (23.0%).
