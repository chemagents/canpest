# Split split_04 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 14.38% | 43 | 0.9253 | 0.9343 | 0.8690 |
| CatBoost RMT95 standalone | 30.10% | 90 | 0.8238 | 0.8375 | 0.7946 |
| SVM RMT95 standalone | 30.77% | 92 | 0.8115 | 0.8180 | 0.7946 |
| DMPNN x CatBoost RMT95 | 6.02% | 18 | 0.9188 | 0.9236 | 0.6667 |

False positives removed by veto: 25/43 (58.1%). True positives lost: 68/292 (23.3%).
