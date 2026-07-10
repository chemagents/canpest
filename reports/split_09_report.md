# Split split_09 RMT95 Veto Reproduction

RMT95 features: 95 residues (45 GLC1 + 50 ACHE). RTE390 was not used.

Formula: `p_final = p_DMPNN * p_CatBoost_RMT95`.

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 13.38% | 40 | 0.9334 | 0.9422 | 0.8988 |
| CatBoost RMT95 standalone | 38.46% | 115 | 0.7987 | 0.8179 | 0.8363 |
| SVM RMT95 standalone | 36.45% | 109 | 0.7972 | 0.8003 | 0.8036 |
| DMPNN x CatBoost RMT95 | 5.69% | 17 | 0.9214 | 0.9320 | 0.7262 |

False positives removed by veto: 23/40 (57.5%). True positives lost: 58/302 (19.2%).
