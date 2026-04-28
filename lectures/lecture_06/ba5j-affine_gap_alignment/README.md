# Affine Gap Alignment (BA5J)

## Задача

Найти наилучшее **глобальное выравнивание** двух аминокислотных последовательностей, где стоимость гэпа зависит не от длины гэпа пропорционально, а имеет фиксированную стоимость «открытия» и меньшую стоимость «продления».

**Вход:** два пептида
**Выход:** максимальный счёт; две выровненные строки

Формат входного файла Rosalind:
```
PRTEINS
PRTWPSEIN
```

**Пример:**
```
Вход:       Вывод:
PRTEINS     8
PRTWPSEIN   PRT---EINS
            PRTWPSEIN-
```

Проверяем счёт вручную (BLOSUM62, sigma=11, epsilon=1):
- P vs P: 7
- R vs R: 5
- T vs T: 5
- `---` vs `WPS`: гэп длиной 3 открыт в v, стоит -(11 + 2*1) = -13
- E vs E: 5
- I vs I: 4
- N vs N: 6
- S vs `-`: гэп длиной 1, стоит -11
- Итого: 7+5+5-13+5+4+6-11 = **8** ✓

## Теория

### Проблема линейного штрафа за гэп

В простом (линейном) выравнивании гэп длиной k стоит `k * sigma`. Это приводит к биологически нереалистичному результату: алгоритм предпочитает много коротких гэпов одному длинному. На практике одна делеция (потеря сразу 10 нуклеотидов) намного вероятнее, чем 10 отдельных делеций.

**Аффинный штраф** решает эту проблему: гэп длиной k стоит `sigma + (k-1) * epsilon`, где:
- `sigma` — штраф за **открытие** гэпа (дорого: новый гэп — редкое событие)
- `epsilon` — штраф за **продление** (дёшево: раз уж открылся, один символ больше — не так страшно)

С sigma=11, epsilon=1: гэп длиной 1 стоит 11, длиной 3 стоит 13 (а не 33 как при линейном). Длинные гэпы становятся относительно «дешевле».

### Три матрицы динамического программирования

Стандартный DP с одной матрицей не умеет различать «мы продолжаем гэп» и «мы открываем новый гэп». Нужно три матрицы, каждая из которых отвечает за своё «состояние» выравнивания:

**`lower[i][j]`** — лучший счёт, если последний шаг — **гэп в v** (ребро влево, j увеличился, i не изменился). Иными словами, символ `w[j-1]` выровнен против `-` в v.

**`upper[i][j]`** — лучший счёт, если последний шаг — **гэп в w** (ребро вверх, i увеличился, j не изменился). Символ `v[i-1]` выровнен против `-` в w.

**`middle[i][j]`** — лучший счёт, если последний шаг — **совпадение или замена** (диагональный шаг). Оба символа `v[i-1]` и `w[j-1]` выровнены друг против друга.

Переходы между матрицами:
```
lower[i][j] = max(
    lower[i][j-1]  - epsilon,   // продление гэпа: уже в lower, идём дальше влево
    middle[i][j-1] - sigma      // открытие гэпа: из middle открываем новый гэп в v
)

upper[i][j] = max(
    upper[i-1][j]  - epsilon,   // продление гэпа: уже в upper, идём дальше вверх
    middle[i-1][j] - sigma      // открытие гэпа: из middle открываем новый гэп в w
)

middle[i][j] = max(
    lower[i][j],                              // закрываем гэп в v, переходим к совпадению
    upper[i][j],                              // закрываем гэп в w, переходим к совпадению
    middle[i-1][j-1] + BLOSUM62[v[i]][w[j]]  // диагональ: совпадение/замена
)
```

### Матрица BLOSUM62

BLOSUM62 построена на основе реальных выравниваний блоков белковых последовательностей с идентичностью около 62%. Клетка `BLOSUM62[a][b]` — log-odds score: положительный, если замена a→b встречается чаще случайной, отрицательный — реже. Используется в BLAST по умолчанию для белков умеренного сходства.

### Инициализация

Первая строка (i=0): мы выровняли 0 символов v с j символами w — это один гэп в v длиной j:
```
lower[0][j]  = -sigma - (j-1) * epsilon  (гэп открыт один раз, продлён j-1 раз)
middle[0][j] = lower[0][j]               (best достигается через lower)
```

Первая колонка (j=0): аналогично для гэпов в w:
```
upper[i][0]  = -sigma - (i-1) * epsilon
middle[i][0] = upper[i][0]
```

## Применение

**BLAST** использует именно BLOSUM62 с sigma=11, epsilon=1 — те же параметры, что в этой задаче. Аффинный гэп позволяет лучше находить белки с вставками/делециями доменов.

**EMBOSS Needle** — инструмент глобального выравнивания с аффинными гэпами. Применяется для детального сравнения ортологичных белков из разных организмов.

## Алгоритм

```
Инициализация:
  middle[0][0] = 0
  lower[0][j]  = -sigma - (j-1)*epsilon,  middle[0][j] = lower[0][j]
  upper[i][0]  = -sigma - (i-1)*epsilon,  middle[i][0] = upper[i][0]

for i = 1..n:
    for j = 1..m:
        lower[i][j]  = max(lower[i][j-1]-epsilon, middle[i][j-1]-sigma)
        upper[i][j]  = max(upper[i-1][j]-epsilon, middle[i-1][j]-sigma)
        middle[i][j] = max(lower[i][j], upper[i][j],
                           middle[i-1][j-1] + BLOSUM62[v[i]][w[j]])

Ответ: max(lower[n][m], upper[n][m], middle[n][m])
Backtrack от той матрицы, где достигнут максимум.
```

## Реализация на C++

```cpp
// 1. Пять констант для backtrack-матриц.
//    Их значения (0..4) кодируют, ОТКУДА пришли в текущую ячейку:
static const int8_t FROM_LOWER = 0;  // lower←lower:  продлили гэп в v
static const int8_t FROM_UPPER = 1;  // upper←upper:  продлили гэп в w
static const int8_t FROM_DIAG  = 2;  // middle←diag:  совпадение/замена
static const int8_t OPEN_LOWER = 3;  // lower←middle: открыли новый гэп в v
static const int8_t OPEN_UPPER = 4;  // upper←middle: открыли новый гэп в w

// 2. Три матрицы DP. Инициализируем значением NEG_INF = INT_MIN/2,
//    а не нулём, чтобы "недостижимые" ячейки не влияли на max.
//    (INT_MIN/2 — безопасный "минус бесконечность": вычитание sigma
//    из него не приведёт к переполнению.)
std::vector<std::vector<int>> lower(n+1, std::vector<int>(m+1, NEG_INF));
std::vector<std::vector<int>> middle(n+1, std::vector<int>(m+1, NEG_INF));
std::vector<std::vector<int>> upper(n+1, std::vector<int>(m+1, NEG_INF));

// 3. Три backtrack-матрицы типа int8_t (1 байт на ячейку вместо 4).
//    backL[i][j] — откуда пришли в lower[i][j]
//    backU[i][j] — откуда пришли в upper[i][j]
//    backM[i][j] — откуда пришли в middle[i][j]
std::vector<std::vector<int8_t>> backL(n+1, std::vector<int8_t>(m+1, -1));
std::vector<std::vector<int8_t>> backU(n+1, std::vector<int8_t>(m+1, -1));
std::vector<std::vector<int8_t>> backM(n+1, std::vector<int8_t>(m+1, -1));

// 4. Заполнение lower[i][j]: пришли из lower[i][j-1] (продление)
//    или из middle[i][j-1] (открытие нового гэпа).
//    NEG_INF-защита: если предшественник недостижим, не берём его.
int lFromL = (lower[i][j-1]  != NEG_INF) ? lower[i][j-1]  - epsilon : NEG_INF;
int lFromM = (middle[i][j-1] != NEG_INF) ? middle[i][j-1] - sigma   : NEG_INF;
if (lFromL >= lFromM) { lower[i][j] = lFromL; backL[i][j] = FROM_LOWER; }
else                  { lower[i][j] = lFromM; backL[i][j] = OPEN_LOWER; }

// 5. Аналогично для upper[i][j]: пришли из upper[i-1][j] или middle[i-1][j].
int uFromU = (upper[i-1][j]  != NEG_INF) ? upper[i-1][j]  - epsilon : NEG_INF;
int uFromM = (middle[i-1][j] != NEG_INF) ? middle[i-1][j] - sigma   : NEG_INF;
if (uFromU >= uFromM) { upper[i][j] = uFromU; backU[i][j] = FROM_UPPER; }
else                  { upper[i][j] = uFromM; backU[i][j] = OPEN_UPPER; }

// 6. middle[i][j]: лучшее из трёх вариантов.
//    Диагональный предшественник middle[i-1][j-1] + blosum(v[i-1], w[j-1]).
int mFromD = (middle[i-1][j-1] != NEG_INF)
               ? middle[i-1][j-1] + blosum(v[i-1], w[j-1])
               : NEG_INF;
middle[i][j] = std::max({lower[i][j], upper[i][j], mFromD});
// backM записывает, из какого предшественника пришли.
if      (mFromD >= lower[i][j] && mFromD >= upper[i][j]) backM[i][j] = FROM_DIAG;
else if (lower[i][j] >= upper[i][j])                     backM[i][j] = FROM_LOWER;
else                                                      backM[i][j] = FROM_UPPER;

// 7. Определяем стартовую матрицу backtrack: та, в которой dp[n][m] максимально.
int8_t startState;
if (middle[n][m] >= lower[n][m] && middle[n][m] >= upper[n][m])
    { bestScore = middle[n][m]; startState = FROM_DIAG; }  // стартуем в middle
else if (lower[n][m] >= upper[n][m])
    { bestScore = lower[n][m];  startState = FROM_LOWER; } // стартуем в lower
else
    { bestScore = upper[n][m];  startState = FROM_UPPER; } // стартуем в upper

// 8. curState — текущая матрица при backtrack: 0=lower, 1=upper, 2=middle.
//    Начальное значение: FROM_LOWER→0, FROM_UPPER→1, FROM_DIAG→2 (middle).
int curState = (startState == FROM_LOWER) ? 0 : (startState == FROM_UPPER) ? 1 : 2;

while (i > 0 || j > 0) {
    if (curState == 0) {
        // В lower: последний символ — гэп в v.
        alignV += '-';
        alignW += w[j-1];
        --j;
        // Читаем backL[i][j+1] (ячейку ДО шага): узнаём, продление это или открытие.
        // Если FROM_LOWER — продолжаем в lower. Если OPEN_LOWER — переходим в middle.
        int8_t bt = backL[i][j+1];
        curState = (bt == FROM_LOWER) ? 0 : 2;
    } else if (curState == 1) {
        // В upper: последний символ — гэп в w.
        alignV += v[i-1];
        alignW += '-';
        --i;
        int8_t bt = backU[i+1][j];
        curState = (bt == FROM_UPPER) ? 1 : 2;
    } else {
        // В middle: читаем backM[i][j].
        if (backM[i][j] == FROM_DIAG) {
            // Диагональ: оба символа реальные.
            alignV += v[i-1]; alignW += w[j-1]; --i; --j;
            // curState остаётся 2 (middle)
        } else if (backM[i][j] == FROM_LOWER) {
            // Переходим в lower — на следующей итерации обработаем гэп в v.
            curState = 0;
        } else {
            // Переходим в upper — на следующей итерации обработаем гэп в w.
            curState = 1;
        }
    }
}
```

## Разбор кода

**Почему три матрицы, а не одна?**
Для аффинного штрафа нужно знать, продолжается ли гэп или открывается новый. Одна матрица не хранит эту информацию. Lower и upper — «промежуточные» матрицы: они отслеживают, что мы сейчас «внутри» гэпа. Middle — «нейтральное» состояние между гэпами.

**`int8_t` вместо `int` для backtrack-матриц.** У нас 5 возможных значений (0..4). Каждое вмещается в 1 байт. Три матрицы по 1000×1000 ячеек = 3 МБ вместо 12 МБ. При n=m=1000 это заметная экономия.

**`NEG_INF = INT_MIN / 2`.** Если бы мы взяли `INT_MIN` и вычли `sigma`, получили бы переполнение знакового int (undefined behavior в C++). Деление на 2 даёт достаточно большое отрицательное число, при этом `NEG_INF - sigma` всё ещё помещается в int.

**Почему backtrack читает `backL[i][j+1]`, а не `backL[i][j]`?**
Когда мы в lower на шаге `(i, j+1)`, мы делаем `--j`, переходя в `(i, j)`. Нас интересует, как мы попали в `(i, j+1)` — это записано в `backL[i][j+1]`. После декремента j этот индекс становится `backL[i][j+1]` в новых переменных.

**Переход middle → lower без движения.** Когда `backM[i][j] == FROM_LOWER`, мы просто меняем `curState = 0`, не двигая i и j. Это корректно: следующая итерация будет в lower на той же позиции `(i, j)` и сразу сделает шаг влево.

## Запуск

```bash
./ba5j-affine_gap_alignment                              # sample PRTEINS/PRTWPSEIN
./ba5j-affine_gap_alignment --input rosalind_ba5j.txt   # Rosalind вход
./ba5j-affine_gap_alignment -i in.txt -o out.txt
```

## Анализ сложности

**Время: O(n * m)** — три раздельных двойных цикла (фактически один с тремя блоками), каждая ячейка вычисляется за O(1). Примерно в 3 раза медленнее, чем линейный гэп, но асимптотика та же. Узкое место: двойной цикл `for i... for j` в `affineGapAlignment`.

**Память: O(n * m)** — шесть матриц: three DP (`int`) + три backtrack (`int8_t`). Для n=m=1000: 3×4 МБ + 3×1 МБ = ~15 МБ.
