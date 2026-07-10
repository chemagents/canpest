# Split split_26 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 13.04% | 39 | 0.9218 | 0.9310 | 0.8542 |
| CatBoost RMT95 standalone | 42.47% | 127 | 0.7477 | 0.7638 | 0.7976 |
| SVM RMT95 standalone | 36.12% | 108 | 0.7647 | 0.7816 | 0.7589 |
| DMPNN x CatBoost RMT95 | 4.01% | 12 | 0.9151 | 0.9301 | 0.6607 |

False positives removed by veto: 27/39 (69.2%). True positives lost: 65/287 (22.6%).
