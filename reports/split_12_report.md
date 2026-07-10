# Split split_12 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.71% | 35 | 0.9137 | 0.9325 | 0.8333 |
| CatBoost RMT95 standalone | 38.13% | 114 | 0.7624 | 0.7602 | 0.8065 |
| SVM RMT95 standalone | 34.11% | 102 | 0.7766 | 0.7698 | 0.7768 |
| DMPNN x CatBoost RMT95 | 3.68% | 11 | 0.9011 | 0.9101 | 0.6399 |

False positives removed by veto: 24/35 (68.6%). True positives lost: 65/280 (23.2%).
