# Split split_07 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 15.05% | 45 | 0.9204 | 0.9405 | 0.8631 |
| CatBoost RMT95 standalone | 39.46% | 118 | 0.7804 | 0.7708 | 0.8393 |
| SVM RMT95 standalone | 31.44% | 94 | 0.8097 | 0.8123 | 0.8006 |
| DMPNN x CatBoost RMT95 | 5.69% | 17 | 0.9104 | 0.9126 | 0.7024 |

False positives removed by veto: 28/45 (62.2%). True positives lost: 54/290 (18.6%).
