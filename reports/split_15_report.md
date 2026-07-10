# Split split_15 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 12.37% | 37 | 0.9288 | 0.9418 | 0.8571 |
| CatBoost RMT95 standalone | 31.77% | 95 | 0.8038 | 0.7931 | 0.8393 |
| SVM RMT95 standalone | 32.11% | 96 | 0.7849 | 0.7566 | 0.8006 |
| DMPNN x CatBoost RMT95 | 5.35% | 16 | 0.9147 | 0.9182 | 0.7113 |

False positives removed by veto: 21/37 (56.8%). True positives lost: 49/288 (17.0%).
