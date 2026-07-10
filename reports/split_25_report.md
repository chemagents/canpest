# Split split_25 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 10.03% | 30 | 0.9343 | 0.9517 | 0.8482 |
| CatBoost RMT95 standalone | 29.77% | 89 | 0.8259 | 0.8311 | 0.7946 |
| SVM RMT95 standalone | 28.43% | 85 | 0.8109 | 0.8136 | 0.7619 |
| DMPNN x CatBoost RMT95 | 4.01% | 12 | 0.9304 | 0.9400 | 0.7054 |

False positives removed by veto: 18/30 (60.0%). True positives lost: 48/285 (16.8%).
