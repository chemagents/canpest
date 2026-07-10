# Split split_13 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 12.71% | 38 | 0.9194 | 0.9378 | 0.8155 |
| CatBoost RMT95 standalone | 35.79% | 107 | 0.7753 | 0.7830 | 0.7976 |
| SVM RMT95 standalone | 33.78% | 101 | 0.7847 | 0.7813 | 0.8036 |
| DMPNN x CatBoost RMT95 | 3.68% | 11 | 0.9126 | 0.9215 | 0.6905 |

False positives removed by veto: 27/38 (71.1%). True positives lost: 42/274 (15.3%).
