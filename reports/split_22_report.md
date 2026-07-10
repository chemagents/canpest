# Split split_22 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 9.36% | 28 | 0.9344 | 0.9454 | 0.8542 |
| CatBoost RMT95 standalone | 34.11% | 102 | 0.8145 | 0.8144 | 0.8095 |
| SVM RMT95 standalone | 33.11% | 99 | 0.8052 | 0.8054 | 0.7887 |
| DMPNN x CatBoost RMT95 | 4.35% | 13 | 0.9335 | 0.9445 | 0.7054 |

False positives removed by veto: 15/28 (53.6%). True positives lost: 50/287 (17.4%).
