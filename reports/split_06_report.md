# Split split_06 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.37% | 34 | 0.9335 | 0.9492 | 0.8363 |
| CatBoost RMT95 standalone | 38.46% | 115 | 0.7816 | 0.7706 | 0.8214 |
| SVM RMT95 standalone | 32.11% | 96 | 0.7992 | 0.7986 | 0.7946 |
| DMPNN x CatBoost RMT95 | 4.68% | 14 | 0.9203 | 0.9309 | 0.6696 |

False positives removed by veto: 20/34 (58.8%). True positives lost: 56/281 (19.9%).
