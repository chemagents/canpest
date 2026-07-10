# Split split_18 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 10.70% | 32 | 0.9419 | 0.9540 | 0.8780 |
| CatBoost RMT95 standalone | 32.11% | 96 | 0.8176 | 0.8309 | 0.7798 |
| SVM RMT95 standalone | 32.11% | 96 | 0.7934 | 0.7861 | 0.7917 |
| DMPNN x CatBoost RMT95 | 3.68% | 11 | 0.9257 | 0.9393 | 0.6905 |

False positives removed by veto: 21/32 (65.6%). True positives lost: 63/295 (21.4%).
