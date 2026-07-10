# RMT95 Veto Reproduction: подробное воспроизведение DMPNN x CatBoost(RMT95) veto

## Краткий смысл

Эта папка содержит воспроизводимый набор для эксперимента по снижению доли ложноположительных предсказаний пестицидов (`FPR`, false-positive rate) с помощью механистического docking-veto.

Главная идея:

1. `DMPNN(rdkit2d)` смотрит на химическую структуру молекулы и предсказывает вероятность класса `pesticide`.
2. `CatBoost(RMT95)` смотрит только на 95 RMT-отобранных docking residue-energy признаков по двум белкам, `GLC1` и `ACHE`, и оценивает, есть ли механистическая docking-поддержка активности.
3. Итоговая вероятность считается асимметрично:

```text
p_final = p_DMPNN * p_CatBoost_RMT_RTE_SEL
```

Такой product-veto может только уменьшить предсказание DMPNN. Он не может повысить вероятность, если DMPNN сам не видит структурный pesticide-сигнал.

Дополнительно обучается `SVM(RMT_RTE_SEL)` на тех же 95 признаках (RMT_RTE_SEL с отобранными 95 признаками). SVM нужен для диагностического CatBoost/SVM benchmark-а residue-level docking signal. SVM не используется в финальной veto-формуле.

## Что именно воспроизводится

Воспроизводится RMT95-only ветка эксперимента:

```text
DMPNN(rdkit2d) x CatBoost(RMT95)
```

Здесь воспроизводится формула для одновременного использования модели DMPNN и CatBoost(RMT_RTE_SEL):

```text
p_final = p_DMPNN * p_CatBoost_RMT_RTE_SEL
```

Она даёт `FPR = 5.14%` и показывает почти такое же снижение ложноположительных срабатываний, но только на 95 RMT-отобранных residue features.

## Основной результат

30 split-ов, threshold `0.5`, test set каждого split-а.

| Method | FPR mean ± std | ROC AUC mean ± std | PR AUC mean | Recall mean | FP total | FP removed vs DMPNN | TP lost vs DMPNN |
|---|---:|---:|---:|---:|---:|---:|---:|
| DMPNN | 12.20% ± 1.86 | 0.9283 ± 0.0092 | 0.9426 | 0.8514 | 1094 | 0 | 0 |
| CatBoost(RMT95) standalone | 35.69% ± 3.71 | 0.7886 ± 0.0240 | 0.7932 | 0.7995 | 3201 | 331 (30.3%) | 1079 (12.6%) |
| SVM(RMT95) standalone | 33.48% ± 2.58 | 0.7851 ± 0.0164 | 0.7807 | 0.7851 | 3003 | 366 (33.5%) | 1192 (13.9%) |
| DMPNN x CatBoost(RMT95) | 5.14% ± 1.03 | 0.9167 ± 0.0088 | 0.9271 | 0.6801 | 461 | 633 (57.9%) | 1727 (20.1%) |

Главное число:

```text
FPR: 12.20% -> 5.14%
absolute reduction: 7.06 percentage points
relative reduction: 57.9%
```

Иными словами, veto удаляет `633` из `1094` ложноположительных DMPNN-pesticide предсказаний.

Цена фильтра: теряется `1727` из `8582` true positives, то есть `20.1%` истинных положительных DMPNN-предсказаний.

## Биологическая и ML-логика эксперимента

В этой задаче есть две разные информации о молекуле.

Первая информация: химическая структура.

DMPNN получает SMILES и RDKit2D descriptors. Он хорошо классифицирует pesticide/non-pesticide по общей химической структуре. Его ROC AUC высокий: `0.9283`.

Вторая информация: механистическая docking-поддержка.

Если DMPNN говорит, что молекула похожа на pesticide, но docking-профиль по важным остаткам `GLC1` и `ACHE` слабый или не похож на активные соединения, такую молекулу можно понизить. Это и есть veto.

Здесь veto мягкий, не бинарный:

```text
p_final = p_DMPNN * p_CatBoost_RMT95
```

Пример:

| p_DMPNN | p_CatBoost_RMT95 | p_final | Итог при threshold 0.5 |
|---:|---:|---:|---|
| 0.90 | 0.90 | 0.81 | pesticide remains pesticide |
| 0.90 | 0.40 | 0.36 | veto -> not pesticide |
| 0.60 | 0.70 | 0.42 | veto -> not pesticide |
| 0.40 | 0.90 | 0.36 | not pesticide already |

Поэтому CatBoost(RMT95) не заменяет DMPNN. Он работает как физико-химический фильтр согласованности.

## Термины

| Term | Meaning |
|---|---|
| `activity=1` | pesticide / active class |
| `activity=0` | non-pesticide / inactive class |
| `DMPNN` | directed message passing neural network, модель по molecular graph |
| `RDKit2D` | набор расчётных 2D molecular descriptors |
| `RMT` | Random Matrix Theory, метод отделения сигнальных компонент от шумовых в матрице residue energies |
| `RMT95` | 95 RMT-selected residue energy features: 45 GLC1 + 50 ACHE |
| `CatBoost(RMT95)` | CatBoost classifier на 95 RMT-selected features |
| `SVM(RMT95)` | SVM classifier на тех же 95 features, diagnostic model |
| `veto` | понижение вероятности DMPNN при слабой docking-поддержке |
| `FPR` | false positive rate = FP / (FP + TN) |
| `FP` | non-pesticide, ошибочно предсказанный как pesticide |
| `TP lost` | true pesticide, который DMPNN считал positive, но veto опустил ниже threshold |

## Структура папки

```text
./
  README.md
  reproduction_process.md
  select_rmt_residues.py
  run_split_veto.py
  run_30splits_veto.py
  sha256_manifest_30splits.csv

  data/
    data.csv
    splits/
      split_registry.csv
      split_summary.csv
      metadata.json
    dmpnn_predictions/
      split_XX_val_predictions.csv
      split_XX_test_predictions.csv
    rmt95_features/
      split_XX/
        rmt95_train.csv
        rmt95_test.csv
        feature_names.csv

  source_inputs/
    residue_matrices/
      *_residue_matrix.csv
    dmpnn_rdkit2d_cache/
      fp_rdkit2d.npy
      rdkit2d_names.npy

  source_code/
    rmt_filter.py
    train_models.py

  rmt_selection/
    GLC1_ranked_residues.csv
    ACHE_ranked_residues.csv
    GLC1_m_scan.csv
    ACHE_m_scan.csv
    selected_rmt95_features.csv
    diag.json

  models/
    dmpnn_rdkit2d/
      split_XX/
        dmpnn_model.pt
        feature_scaler.pkl
        global_feature_names.npy
        global_feature_names_after_filter.npy
        metrics.json
        val_predictions.csv
        test_predictions.csv
    catboost_rmt95/
      split_XX/
        model.cbm
        feature_names.csv
    svm_rmt95/
      split_XX/
        model.joblib
        feature_names.csv

  metrics/
    split_XX_metrics.csv
    split_XX_summary.json
    per_split_methods_rmt95.csv
    method_summary_30splits.csv
    split_summary_30splits.csv
    summary_30splits.json
    compare_to_entry52_product_rmt95.csv

  predictions/
    split_XX_test_predictions_with_veto.csv
    split_XX_catboost_rmt95_val.csv
    split_XX_svm_rmt95_val.csv
    pooled_30splits_predictions_with_veto.csv

  reports/
    split_XX_report.md
    30splits_report.md
```

## Самодостаточность папки

Эта папка самодостаточна для воспроизведения logged RMT95-veto (RMT_RTE_SEL-veto) эксперимента.

Она содержит:

| Bundle path | Для чего нужно |
|---|---|
| `data/data.csv` | основной датасет с SMILES, activity и docking scores |
| `data/splits/split_registry.csv` | фиксированные 30 split-ов |
| `source_inputs/residue_matrices/` | raw residue-energy matrices для RMT-отбора |
| `source_code/rmt_filter.py` | функции RMT и загрузки residue matrices |
| `source_code/train_models.py` | исходный training script для DMPNN/RDKit2D route |
| `source_inputs/dmpnn_rdkit2d_cache/fp_rdkit2d.npy` | RDKit2D descriptor cache |
| `models/dmpnn_rdkit2d/split_*/` | сохранённые DMPNN baseline models и predictions |
| `models/catboost_rmt95/split_*/` | сохранённые CatBoost(RMT95) models |
| `models/svm_rmt95/split_*/` | сохранённые SVM(RMT95) diagnostic models |

Практическая граница воспроизводимости:

| Уровень воспроизведения | Состояние |
|---|---|
| Пересчитать RMT selection из raw residue matrices | да |
| Заново собрать 95 RMT features | да |
| Заново обучить CatBoost(RMT95) на 30 split-ах | да |
| Заново обучить SVM(RMT95) на 30 split-ах | да |
| Заново посчитать DMPNN x CatBoost veto | да |
| Заново получить все logged FPR/ROC/PR metrics | да |
| Заново переобучить 30 DMPNN weights одной командой | нет, DMPNN weights сохранены, route/source сохранены, но отдельный wrapper не добавлен |

То есть для воспроизведения именно logged FPR-veto эксперимента папки достаточно. Если требовать переобучить DMPNN weights с нуля из SMILES, нужен дополнительный wrapper вокруг `source_code/train_models.py`.

## Требования к окружению

Скрипты рассчитаны на окружение проекта `canpest`.

Проверенное окружение в проекте:

| Package / tool | Purpose |
|---|---|
| Python `~/.pyenv/shims/python3` | основной интерпретатор |
| `pandas` | таблицы |
| `numpy` | матрицы |
| `scipy` | point-biserial correlation |
| `scikit-learn` | LogisticRegression, SVM, metrics, StandardScaler |
| `catboost` | CatBoostClassifier |
| `joblib` | сохранение SVM model |
| `rdkit` | DMPNN/RDKit route, если переобучать DMPNN |
| `torch` | DMPNN model weights, если переобучать DMPNN |

Все команды ниже запускать из корня этого репозитория:

```bash
cd <repo-root>
```

## Входные данные

### Основной датасет

Файл:

```text
data/data.csv
```

Содержит `5920` строк.

Ключевые колонки:

| Column | Meaning |
|---|---|
| `ligand_id` | ID молекулы |
| `SMILES` | химическая структура |
| `activity` | `1` pesticide, `0` non-pesticide, `NaN` unlabelled |
| `2imi` | GSTE2 docking score |
| `d8v7j0` | ACHE docking score |
| `3rif` | GLC1 docking score |
| `8sfy` | UGT202A2 docking score |
| `8udb` | GSTM12 docking score |
| `8v3d` | AGAMOR28 docking score |

Размеченная часть:

| Class | Count |
|---|---:|
| `activity=1` | 1680 |
| `activity=0` | 1491 |
| labelled total | 3171 |

### Split registry

Файл:

```text
data/splits/split_registry.csv
```

Содержит фиксированные split-ы:

| Property | Value |
|---|---:|
| Number of splits | 30 |
| Train per split | 2536 |
| Test per split | 635 |
| Train active | 1344 |
| Train inactive | 1192 |
| Test active | 336 |
| Test inactive | 299 |

Использовать надо именно этот файл. Нельзя генерировать новые split-ы, иначе результаты не будут совпадать.

### Residue matrices

Папка:

```text
source_inputs/residue_matrices/
```

Для RMT95 используются две матрицы:

| Protein | File key | Biological target |
|---|---|---|
| `GLC1` | `glc1` / `3RIF` | glutamate-gated chloride channel |
| `ACHE` | `D8V7J0` | acetylcholinesterase |

В папке сохранены и другие residue matrices, но для текущего RMT95-veto используются только `GLC1` и `ACHE`.

## Шаг 1. RMT-отбор 95 остатков

Скрипт:

```text
select_rmt_residues.py
```

Команда:

```bash
python3 select_rmt_residues.py
```

Что делает скрипт:

1. Загружает labelled molecules из `data/data.csv`.
2. Загружает residue matrices из `source_inputs/residue_matrices/`.
3. Для каждого белка строит матрицу `X`: строки = молекулы, колонки = residue energy features.
4. Нормирует матрицу и считает correlation matrix.
5. Делает eigen decomposition.
6. Считает Marcenko-Pastur threshold:

```text
q = p / n
lambda_plus = (1 + sqrt(q))^2
```

7. Компоненты с eigenvalue выше `lambda_plus` считаются signal components.
8. Для каждого остатка считает RMT prior `s_i`:

```text
s_i = sum(lambda_j * v_ij^2) over signal components j
```

9. Считает point-biserial correlation `r_pb` между residue energy и `activity`.
10. Нормирует `s_i` в `s_norm`.
11. Ранжирует остатки:

```text
rank_score_i = abs(r_pb_i) * s_norm_i
```

12. Подбирает `m_opt` через inner CV по PR-AUC.

Inner CV:

| Parameter | Value |
|---|---:|
| repeats | 3 |
| inner train fraction | 0.7 |
| max_m | 50 |
| metric | PR-AUC |
| probe model | LogisticRegression on cumulative residue-energy score |

Ожидаемый результат:

| Protein | n_residues | n_signal | m_opt | best PR-AUC | top1 |
|---|---:|---:|---:|---:|---|
| GLC1 | 126 | 14 | 45 | 0.7366 | `3rif_ser224` |
| ACHE | 132 | 15 | 50 | 0.7005 | `32264_val396` |

Итог:

```text
GLC1: 45 selected residues
ACHE: 50 selected residues
Total: 95 RMT95 features
```

Основные выходы:

| File | Meaning |
|---|---|
| `rmt_selection/GLC1_ranked_residues.csv` | полный ranking GLC1 residues |
| `rmt_selection/ACHE_ranked_residues.csv` | полный ranking ACHE residues |
| `rmt_selection/GLC1_m_scan.csv` | scan по m для GLC1 |
| `rmt_selection/ACHE_m_scan.csv` | scan по m для ACHE |
| `rmt_selection/selected_rmt95_features.csv` | выбранные 95 features |
| `rmt_selection/diag.json` | RMT diagnostics |

## Шаг 2. Один split для smoke-test

Скрипт:

```text
run_split_veto.py
```

Команда:

```bash
python3 run_split_veto.py --split-id split_00
```

Что делает скрипт:

1. Загружает RMT-selected features из `rmt_selection/selected_rmt95_features.csv`.
2. Загружает raw residue matrices.
3. Собирает `rmt95_train.csv` и `rmt95_test.csv` для указанного split-а.
4. Загружает сохранённые DMPNN validation/test predictions.
5. Исключает DMPNN validation ligand IDs из CatBoost/SVM train.
6. Обучает CatBoost(RMT95).
7. Обучает SVM(RMT95).
8. Сохраняет обе модели.
9. Считает product-veto:

```text
p_final = p_DMPNN * p_CatBoost_RMT95
```

10. Считает metrics при threshold `0.5`.

Ожидаемый результат для `split_00`:

| Method | FPR | FP | ROC AUC | PR AUC | Recall |
|---|---:|---:|---:|---:|---:|
| DMPNN | 11.04% | 33 | 0.9366 | 0.9488 | 0.8631 |
| CatBoost RMT95 standalone | 35.79% | 107 | 0.7815 | 0.7844 | 0.7798 |
| SVM RMT95 standalone | 31.77% | 95 | 0.7876 | 0.7797 | 0.7619 |
| DMPNN x CatBoost RMT95 | 4.68% | 14 | 0.9265 | 0.9265 | 0.6964 |

Для `split_00` veto удаляет:

```text
FP removed: 19 / 33 = 57.6%
TP lost: 56 / 290 = 19.3%
FPR: 11.04% -> 4.68%
```

## Шаг 3. Полный прогон на 30 split-ах

Скрипт:

```text
run_30splits_veto.py
```

Команда:

```bash
python3 run_30splits_veto.py
```

Что делает скрипт для каждого `split_00` ... `split_29`:

1. Берёт train/test ligand IDs из `data/splits/split_registry.csv`.
2. Строит per-split RMT95 feature matrices.
3. Загружает DMPNN validation/test predictions.
4. Делит train на CatBoost/SVM train и validation protocol:

```text
CatBoost/SVM train = split train - DMPNN val ligand IDs
CatBoost eval_set = DMPNN val ligand IDs
SVM train = split train - DMPNN val ligand IDs
test = split test
```

5. Обучает CatBoost(RMT95).
6. Обучает SVM(RMT95).
7. Сохраняет модели.
8. Считает DMPNN baseline metrics.
9. Считает CatBoost standalone metrics.
10. Считает SVM standalone metrics.
11. Считает DMPNN x CatBoost(RMT95) product-veto metrics.
12. Пишет per-split predictions и summary.
13. Собирает aggregate по 30 split-ам.

## Модели и гиперпараметры

### DMPNN(rdkit2d)

DMPNN baseline был обучен раньше и сохранён в bundle.

Путь:

```text
models/dmpnn_rdkit2d/split_*/
```

Важные файлы на каждый split:

| File | Meaning |
|---|---|
| `dmpnn_model.pt` | Torch weights |
| `feature_scaler.pkl` | scaler/filter для RDKit2D features |
| `global_feature_names.npy` | все global feature names |
| `global_feature_names_after_filter.npy` | features после train-only constant filter |
| `metrics.json` | DMPNN split metrics |
| `val_predictions.csv` | validation predictions |
| `test_predictions.csv` | test predictions |

Используемая probability column:

```text
blend_proba
```

Для этого rdkit2d baseline blend effectively equals DMPNN route in the saved predictions.

### CatBoost(RMT95)

Гиперпараметры:

| Parameter | Value |
|---|---:|
| iterations | 500 |
| learning_rate | 0.05 |
| depth | 6 |
| loss_function | `Logloss` |
| eval_metric | `AUC` |
| early_stopping_rounds | 30 |
| random_seed | `42 + split_index` |
| thread_count | 1 |

Сохраняется как:

```text
models/catboost_rmt95/split_XX/model.cbm
```

CatBoost используется в veto formula.

### SVM(RMT95)

Pipeline:

```text
StandardScaler() -> SVC(probability=True)
```

Гиперпараметры:

| Parameter | Value |
|---|---|
| kernel | `rbf` |
| C | `1.0` |
| gamma | `scale` |
| class_weight | `balanced` |
| probability | `True` |
| random_state | `42 + split_index` |

Сохраняется как:

```text
models/svm_rmt95/split_XX/model.joblib
```

SVM не используется в veto. Он нужен, чтобы показать, что residue-level RMT95 features дают сходный самостоятельный signal в CatBoost/SVM diagnostic benchmark.

## Как считается FPR и veto

Threshold для classification:

```text
threshold = 0.5
```

Baseline DMPNN prediction:

```text
dmpnn_pred = p_DMPNN >= 0.5
```

Product-veto prediction:

```text
p_final = p_DMPNN * p_CatBoost_RMT95
veto_pred = p_final >= 0.5
```

False positive:

```text
activity = 0 and predicted class = 1
```

FPR:

```text
FPR = FP / (FP + TN)
```

FP removed:

```text
DMPNN predicted positive, activity=0, but product-veto predicted negative
```

TP lost:

```text
DMPNN predicted positive, activity=1, but product-veto predicted negative
```

## Главные output files

| File | Meaning |
|---|---|
| `metrics/per_split_methods_rmt95.csv` | все metrics по split/method |
| `metrics/method_summary_30splits.csv` | aggregate metrics по 30 split-ам |
| `metrics/split_summary_30splits.csv` | краткая per-split summary |
| `metrics/summary_30splits.json` | machine-readable full summary |
| `metrics/compare_to_entry52_product_rmt95.csv` | проверка совпадения с Entry 52 |
| `predictions/pooled_30splits_predictions_with_veto.csv` | 19050 test predictions по всем split-ам |
| `reports/30splits_report.md` | человекочитаемый итоговый отчёт |
| `sha256_manifest_30splits.csv` | checksum manifest |

## Содержимое pooled predictions

Файл:

```text
predictions/pooled_30splits_predictions_with_veto.csv
```

Строк: `19050` = `30 splits * 635 test molecules`.

Ключевые колонки:

| Column | Meaning |
|---|---|
| `split_id` | split ID |
| `ligand_id` | molecule ID |
| `activity` | true label |
| `p_dmpnn` | DMPNN probability |
| `p_catboost_rmt95` | CatBoost support probability |
| `p_svm_rmt95` | SVM diagnostic probability |
| `p_final_product_rmt95` | product-veto probability |
| `dmpnn_pred` | DMPNN class at threshold 0.5 |
| `veto_pred` | product-veto class at threshold 0.5 |
| `fp_before` | false positive before veto |
| `fp_after` | false positive after veto |
| `tp_lost_by_veto` | true positive removed by veto |
| `vetoed` | DMPNN positive changed to negative by veto |

## Проверка результата

После полного запуска выполнить:

```bash
python3 - <<'PY'
import pandas as pd
from pathlib import Path

base = Path('.')
summary = pd.read_csv(base / 'metrics/method_summary_30splits.csv')
print(summary[['method','n_splits','fpr_mean','roc_auc_mean','fp_total','fp_removed_total','tp_lost_total']].to_string(index=False))

comp = pd.read_csv(base / 'metrics/compare_to_entry52_product_rmt95.csv')
print('max_abs_fpr_diff', comp.fpr_diff.abs().max())
print('max_abs_fp_diff', comp.fp_diff.abs().max())
print('catboost_models', len(list((base / 'models/catboost_rmt95').glob('split_*/model.cbm'))))
print('svm_models', len(list((base / 'models/svm_rmt95').glob('split_*/model.joblib'))))
print('dmpnn_models', len(list((base / 'models/dmpnn_rdkit2d').glob('split_*/dmpnn_model.pt'))))
print('pooled_rows', len(pd.read_csv(base / 'predictions/pooled_30splits_predictions_with_veto.csv')))
PY
```

Ожидаемый вывод по смыслу:

```text
dmpnn                    fpr_mean ~= 0.121962, fp_total = 1094
catboost_rmt95_standalone fpr_mean ~= 0.356856
svm_rmt95_standalone      fpr_mean ~= 0.334783
product_rmt95_veto        fpr_mean ~= 0.051394, fp_total = 461, fp_removed_total = 633

max_abs_fpr_diff ~= 1e-16
max_abs_fp_diff = 0
catboost_models = 30
svm_models = 30
dmpnn_models = 30
pooled_rows = 19050
```

`max_abs_fpr_diff ~= 1e-16` означает, что новое воспроизведение совпадает с историческим Entry 52 на уровне floating-point roundoff.

## Как проверить checksum manifest

Manifest:

```text
sha256_manifest_30splits.csv
```

Проверка:

```bash
python3 - <<'PY'
from pathlib import Path
import hashlib
import pandas as pd

base = Path('.')
manifest = pd.read_csv(base / 'sha256_manifest_30splits.csv')
bad = []
for _, row in manifest.iterrows():
    path = base / row['path']
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    if h.hexdigest() != row['sha256']:
        bad.append(row['path'])
print('bad files:', bad)
print('checked files:', len(manifest))
PY
```

Ожидаемо:

```text
bad files: []
```

## Полная команда воспроизведения

Минимальный полный rerun:

```bash
cd <repo-root>
python3 select_rmt_residues.py
python3 run_30splits_veto.py
```

Если нужен только быстрый smoke-test:

```bash
python3 run_split_veto.py --split-id split_00
```

Если нужно прогнать несколько split-ов:

```bash
python3 run_30splits_veto.py --split-ids split_00,split_01,split_02
```

## Связь с Entry 44, Entry 52 и Entry 53

Entry 44:

```text
RMT + point-biserial selection -> GLC1 top-45 + ACHE top-50 = RMT95
```

Эта папка воспроизводит Entry 44 через:

```text
select_rmt_residues.py
```

Entry 52:

```text
product_rmt95_prev = p_DMPNN * p_CatBoost_RMT95
```

Эта папка воспроизводит Entry 52 product RMT95 branch. Сравнение сохранено в:

```text
metrics/compare_to_entry52_product_rmt95.csv
```

Entry 53:

```text
manuscript wording for FPR reduction
```

Entry 53 в тексте рукописи использует как главный результат более поздний balanced mix с RTE390. Но RMT95-only branch остаётся ключевой компактной mechanistic-veto демонстрацией: 95 признаков дают почти тот же порядок снижения FPR.

## Claims Supported

Поддерживаемые утверждения:

1. DMPNN(rdkit2d) имеет сильное baseline качество: ROC AUC `0.9283`.
2. RMT95 residue features сами по себе несут сигнал: CatBoost ROC AUC `0.7886`, SVM ROC AUC `0.7851`.
3. CatBoost/SVM standalone слабее DMPNN и не должны заменять DMPNN.
4. CatBoost(RMT95) полезен как veto/support module после DMPNN.
5. Product-veto снижает FPR с `12.20%` до `5.14%`.
6. Относительное снижение FPR составляет `57.9%`.
7. Снижение FPR имеет цену: recall падает с `0.8514` до `0.6801`, а TP lost составляет `20.1%`.

## Claims To Avoid

Нельзя утверждать:

1. Нельзя говорить, что CatBoost(RMT95) или SVM(RMT95) лучше DMPNN как standalone classifier.
2. Нельзя говорить, что veto повышает ROC AUC. ROC AUC падает с `0.9283` до `0.9167`.
3. Нельзя скрывать потерю TP: `1727` true positives опускаются ниже threshold.
4. Нельзя называть cannabis/unlabelled veto rate `FPR`, потому что FPR требует true negative labels.
5. Нельзя переносить result `5.14% FPR` на unlabelled cannabis set без true labels.
6. Нельзя смешивать этот RMT95-only rerun с RTE390 geometric mix; RTE390 здесь не используется.

## Типичные ошибки при воспроизведении

### Ошибка 1. Сгенерировать новые split-ы

Нельзя. Надо использовать:

```text
data/splits/split_registry.csv
```

Иначе FPR, FP, TP и все метрики изменятся.

### Ошибка 2. Обучить CatBoost на полном train без DMPNN validation exclusion

Нельзя. Правильный protocol:

```text
CatBoost train = split train - DMPNN val IDs
CatBoost eval_set = DMPNN val IDs
```

Это нужно, чтобы support model validation соответствовала DMPNN validation protocol.

### Ошибка 3. Разрешить docking model повышать DMPNN probability

Нельзя для этого эксперимента. Product-veto должен быть:

```text
p_final <= p_DMPNN
```

Поэтому используется умножение, а не weighted average.

### Ошибка 4. Использовать SVM в veto formula

В этом experiment SVM не используется в veto. Он diagnostic model.

Правильная formula:

```text
p_final = p_DMPNN * p_CatBoost_RMT95
```

### Ошибка 5. Добавить RTE390

В этой папке RTE390 намеренно исключён. Если добавить RTE390, это будет другой experiment.

### Ошибка 6. Запускать не из корня репозитория

Рекомендуемый запуск:

```bash
cd <repo-root>
python3 run_30splits_veto.py
```

## Troubleshooting

### `ModuleNotFoundError: catboost`

Окружение не содержит CatBoost. Нужен `catboost`, потому что `model.cbm` обучается заново.

### `ModuleNotFoundError: rdkit`

Для текущего veto rerun RDKit обычно не нужен, потому что DMPNN predictions уже сохранены. Но RDKit нужен, если переобучать DMPNN через `source_code/train_models.py`.

### `FileNotFoundError` для residue matrices

Проверить, что существует:

```text
source_inputs/residue_matrices/
```

### `FileNotFoundError` для DMPNN predictions

Проверить, что существуют:

```text
models/dmpnn_rdkit2d/split_XX/val_predictions.csv
models/dmpnn_rdkit2d/split_XX/test_predictions.csv
```

### Результаты не совпадают с README

Проверить три вещи:

1. Используется `threshold = 0.5`.
2. Используется fixed `split_registry.csv` из `data/splits/`.
3. CatBoost seed равен `42 + split_index`.

## Проверенная финальная консистентность

Последняя проверка bundle:

```text
CatBoost models: 30
SVM models: 30
DMPNN snapshots: 30
RMT train matrices: 30
Residue matrices: 6
Pooled test prediction rows: 19050
Entry 52 FP/TP diff: 0
```

Финальная метрика:

```text
DMPNN FPR = 0.12196209587513936
Product RMT95 veto FPR = 0.05139353400222965
Relative FPR reduction = 0.5786106032906765
```
