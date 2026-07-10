# Split split_27 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 13.71% | 41 | 0.9324 | 0.9469 | 0.8690 |
| CatBoost RMT95 standalone | 35.12% | 105 | 0.8061 | 0.8251 | 0.7946 |
| SVM RMT95 standalone | 32.44% | 97 | 0.7934 | 0.7822 | 0.7798 |
| DMPNN x CatBoost RMT95 | 5.35% | 16 | 0.9189 | 0.9344 | 0.6935 |

False positives removed by veto: 25/41 (61.0%). True positives lost: 59/292 (20.2%).
