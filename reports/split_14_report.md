# Split split_14 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 10.37% | 31 | 0.9300 | 0.9474 | 0.8393 |
| CatBoost RMT95 standalone | 38.13% | 114 | 0.7756 | 0.7835 | 0.7857 |
| SVM RMT95 standalone | 36.79% | 110 | 0.7684 | 0.7675 | 0.7798 |
| DMPNN x CatBoost RMT95 | 4.68% | 14 | 0.9174 | 0.9350 | 0.6845 |

False positives removed by veto: 17/31 (54.8%). True positives lost: 52/282 (18.4%).
