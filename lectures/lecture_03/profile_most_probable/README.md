# Profile-Most Probable k-mer (BA2C)

## Задача

Найти k-мер в строке с максимальной вероятностью по заданной профильной матрице.

**Вход:** строка `Text`, число `k`, профильная матрица `4×k`
**Выход:** k-мер из `Text` с наибольшей вероятностью

Где:
- `Text` — строка, в которой ищется лучший k-мер
- `k` — длина искомого k-мера
- `Profile` — матрица вероятностей размера `4×k`
- строки `Profile` идут в порядке `A`, `C`, `G`, `T`

Формат входного файла Rosalind:
```text
Text
k
row_A
row_C
row_G
row_T
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `Text k "row_A" "row_C" "row_G" "row_T"`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `profile_most_probable_output.txt`

**Пример:**
```
Text: ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
k: 5
Profile:
  A: 0.2 0.2 0.3 0.2 0.3
  C: 0.4 0.3 0.1 0.5 0.1
  G: 0.3 0.3 0.5 0.2 0.4
  T: 0.1 0.2 0.1 0.1 0.2
Вывод: CCGAG
```

## Теория

### Профильная матрица

**Profile** — матрица 4×k, где строки соответствуют нуклеотидам (A, C, G, T), а столбцы — позициям. Элемент Profile[c][j] — вероятность нуклеотида c на позиции j.

**Вероятность k-мера по профилю:**
```
Pr(Pattern | Profile) = ∏ Profile[Pattern[i]][i]
                         i=0..k-1
```

**Пример для CCGAG:**
```
Pr(CCGAG) = Profile[C][0] × Profile[C][1] × Profile[G][2] × Profile[A][3] × Profile[G][4]
          = 0.4 × 0.3 × 0.5 × 0.2 × 0.4
          = 0.0048
```

**Применение:**
- Основа для жадного и рандомизированного поиска мотивов
- Профиль кодирует «консенсусные предпочтения» набора мотивов
- Profile-most probable k-мер — лучший кандидат для продолжения мотива

## Алгоритм

### Псевдокод
```
ProfileMostProbable(Text, k, Profile):
    bestKmer ← Text[0..k-1]
    bestProb ← -1
    for i ← 0 to |Text| - k:
        kmer ← Text[i..i+k-1]
        prob ← Pr(kmer | Profile)
        if prob > bestProb:
            bestProb ← prob
            bestKmer ← kmer
    return bestKmer
```

## Реализация на C++

```cpp
double probability(const std::string& kmer,
                   const std::vector<std::vector<double>>& profile) {
    double prob = 1.0;
    for (size_t i = 0; i < kmer.size(); ++i) {
        prob *= profile[nucleotideIndex(kmer[i])][i];  // (1)
    }
    return prob;
}

std::string profileMostProbable(const std::string& text, size_t k,
                                const std::vector<std::vector<double>>& profile) {
    double maxProb = -1.0;  // (2)
    std::string bestKmer = text.substr(0, k);
    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string kmer = text.substr(i, k);
        double prob = probability(kmer, profile);
        if (prob > maxProb) {
            maxProb = prob;
            bestKmer = kmer;
        }
    }
    return bestKmer;
}
```

### Особенности реализации

**1. `std::vector<std::vector<double>>` — двумерная матрица**
```cpp
const std::vector<std::vector<double>>& profile
```
- «Вектор векторов» — способ представления 2D-матрицы в C++
- `profile[row][col]` — доступ: сначала к внешнему вектору (строка), потом к внутреннему (столбец)
- `profile` имеет 4 строки (A, C, G, T) и k столбцов (позиции)
- Каждая строка — отдельный `vector<double>` в куче (не непрерывная память!)
- **Альтернативы:** `std::array<std::array<double, K>, 4>` (на стеке, но K должен быть constexpr), или плоский `vector<double>` с индексом `row * k + col`

**2. Тип `double` — вещественные числа**
```cpp
double prob = 1.0;
prob *= profile[nucleotideIndex(kmer[i])][i];
```
- `double` — 64-битное число с плавающей точкой (IEEE 754)
- Точность: ~15-17 значащих десятичных цифр
- `1.0` — литерал `double` (не `1.0f`, что было бы `float`)
- `*=` — составной оператор умножения: `prob = prob * value`
- **Потеря точности при умножении:** произведение k вероятностей < 1 быстро стремится к 0. Для k > 50 может потребоваться логарифмическая шкала: `logProb += log(profile[...][i])`

**3. Двойная индексация**
```cpp
profile[nucleotideIndex(kmer[i])][i]
```
- `kmer[i]` — символ на позиции i (тип `char`)
- `nucleotideIndex(kmer[i])` — маппинг A→0, C→1, G→2, T→3 (номер строки)
- `profile[row][col]` — вероятность нуклеотида `row` на позиции `col`
- Весь вызов: «вероятность символа `kmer[i]` на позиции `i` по профилю»

**4. Начальное значение -1.0 для поиска максимума**
```cpp
double maxProb = -1.0;
```
- Все вероятности ≥ 0.0, поэтому любой k-мер лучше начального значения
- Гарантирует, что `bestKmer` будет обновлён хотя бы один раз
- **Альтернатива:** вычислить вероятность первого k-мера и использовать её как начальное

**Важно:** при нулевой вероятности хотя бы на одной позиции весь k-мер получает вероятность 0. Это мотивирует использование псевдосчётчиков (см. BA2E).

## Запуск

```bash
# Rosalind-style аргументы
./profile_most_probable --output answer.txt \
  "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT" 5 \
  "0.2 0.2 0.3 0.2 0.3" \
  "0.4 0.3 0.1 0.5 0.1" \
  "0.3 0.3 0.5 0.2 0.4" \
  "0.1 0.2 0.1 0.1 0.2"

# Входной файл в формате Rosalind
./profile_most_probable --input rosalind_ba2c.txt --output rosalind_ba2c.out

# Короткая форма: один путь трактуется как входной файл
./profile_most_probable rosalind_ba2c.txt
```

## Сложность

### Время: O((n - k) × k)
- Один проход по строке, для каждого k-мера — вычисление вероятности O(k)

### Память: O(k)
- Хранение текущего лучшего k-мера
