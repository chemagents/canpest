# Split split_19 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 15.72% | 47 | 0.9207 | 0.9392 | 0.8482 |
| CatBoost RMT95 standalone | 36.45% | 109 | 0.7900 | 0.8036 | 0.8125 |
| SVM RMT95 standalone | 32.44% | 97 | 0.7887 | 0.7893 | 0.7917 |
| DMPNN x CatBoost RMT95 | 5.69% | 17 | 0.9141 | 0.9318 | 0.7113 |

False positives removed by veto: 30/47 (63.8%). True positives lost: 46/285 (16.1%).
