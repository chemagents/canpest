# Split split_03 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 12.37% | 37 | 0.9060 | 0.9216 | 0.8274 |
| CatBoost RMT95 standalone | 45.15% | 135 | 0.7288 | 0.7225 | 0.7946 |
| SVM RMT95 standalone | 38.80% | 116 | 0.7445 | 0.7310 | 0.7738 |
| DMPNN x CatBoost RMT95 | 5.35% | 16 | 0.9011 | 0.9086 | 0.5893 |

False positives removed by veto: 21/37 (56.8%). True positives lost: 80/278 (28.8%).
