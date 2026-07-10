# Split split_05 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 9.70% | 29 | 0.9187 | 0.9371 | 0.8036 |
| CatBoost RMT95 standalone | 35.45% | 106 | 0.7908 | 0.7949 | 0.8185 |
| SVM RMT95 standalone | 34.78% | 104 | 0.7699 | 0.7619 | 0.7798 |
| DMPNN x CatBoost RMT95 | 6.35% | 19 | 0.9105 | 0.9275 | 0.6786 |

False positives removed by veto: 10/29 (34.5%). True positives lost: 42/270 (15.6%).
