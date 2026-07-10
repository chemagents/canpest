# Split split_23 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 13.04% | 39 | 0.9343 | 0.9474 | 0.8690 |
| CatBoost RMT95 standalone | 34.78% | 104 | 0.7947 | 0.8017 | 0.8036 |
| SVM RMT95 standalone | 34.11% | 102 | 0.7790 | 0.7617 | 0.7798 |
| DMPNN x CatBoost RMT95 | 5.02% | 15 | 0.9167 | 0.9284 | 0.6696 |

False positives removed by veto: 24/39 (61.5%). True positives lost: 67/292 (22.9%).
