# Split split_01 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.04% | 33 | 0.9219 | 0.9397 | 0.8542 |
| CatBoost RMT95 standalone | 37.79% | 113 | 0.7889 | 0.8019 | 0.7917 |
| SVM RMT95 standalone | 36.45% | 109 | 0.7846 | 0.7834 | 0.7917 |
| DMPNN x CatBoost RMT95 | 6.69% | 20 | 0.9115 | 0.9294 | 0.6696 |

False positives removed by veto: 13/33 (39.4%). True positives lost: 62/287 (21.6%).
