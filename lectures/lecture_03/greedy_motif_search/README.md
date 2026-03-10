# Greedy Motif Search with Pseudocounts (BA2E)

## Задача

Найти набор мотивов жадным алгоритмом с использованием псевдосчётчиков Лапласа.

**Вход:** целые числа `k`, `t` и набор строк `Dna`
**Выход:** набор k-меров (по одному из каждой строки), минимизирующий счёт

Где:
- `k` — длина искомого мотива
- `t` — количество строк в коллекции `Dna`
- `Dna` — набор ДНК-строк

Формат входного файла Rosalind:
```text
k t
dna1
dna2
...
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `k t dna1 dna2 ... dna_t`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `greedy_motif_search_output.txt`

**Пример:**
```
k: 3, t: 5
Dna: GGCGTTCAGGCA  AAGAATCAGTCA  CAAGGAGTTCGC  CACGTCAATCAC  CAATAATATTCG
Вывод:
TTC
ATC
TTC
ATC
TTC
```

## Теория

### Счёт мотивов (Score)

Для набора мотивов Score = суммарное число несовпадений с консенсусом:
```
Score(Motifs) = Σ (t - max_count[j])  для j = 0..k-1
```

### Псевдосчётчики Лапласа

Без псевдосчётчиков: если нуклеотид не встречается ни разу на позиции, его вероятность = 0, и любой k-мер с таким нуклеотидом будет иметь вероятность 0.

**С псевдосчётчиками:** добавляем 1 к каждому счётчику перед нормализацией:
```
Profile[c][j] = (count[c][j] + 1) / (t + 4)
```

Это предотвращает нулевые вероятности и делает алгоритм более устойчивым.

### Жадный подход

1. Фиксируем k-мер из первой строки
2. Для каждой следующей строки строим профиль из уже выбранных мотивов
3. Находим Profile-most probable k-мер в следующей строке
4. Перебираем все стартовые k-меры из первой строки, выбираем лучший набор

## Алгоритм

### Псевдокод
```
GreedyMotifSearch(Dna, k, t):
    BestMotifs ← первые k-меры каждой строки
    for each k-mer Motif_0 in Dna[0]:
        Motifs ← {Motif_0}
        for j ← 1 to t-1:
            Profile ← ProfileWithPseudocounts(Motifs)
            Motifs ← Motifs ∪ {ProfileMostProbable(Dna[j], k, Profile)}
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs ← Motifs
    return BestMotifs
```

## Реализация на C++

```cpp
std::vector<std::vector<double>> buildProfileWithPseudocounts(
        const std::vector<std::string>& motifs) {
    size_t k = motifs[0].size();
    size_t t = motifs.size();
    std::vector<std::vector<double>> profile(4, std::vector<double>(k, 1.0));  // (1)

    for (const auto& motif : motifs)
        for (size_t j = 0; j < k; ++j)
            ++profile[nucleotideIndex(motif[j])][j];

    double total = static_cast<double>(t) + 4.0;  // (2)
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < k; ++j)
            profile[i][j] /= total;

    return profile;
}
```

### Особенности реализации

**1. Инициализация матрицы единицами — псевдосчётчики**
```cpp
std::vector<std::vector<double>> profile(4, std::vector<double>(k, 1.0));
```
- Конструктор `vector(n, val)`: создаёт 4 строки, каждая — вектор из k элементов 1.0
- `1.0` — начальный псевдосчётчик Лапласа для каждого нуклеотида на каждой позиции
- Без псевдосчётчиков: `...(k, 0.0)` — нулевая вероятность убивает весь k-мер
- Математически: добавляем по 1 виртуальному наблюдению каждого нуклеотида

**2. Нормализация: `static_cast<double>(t) + 4.0`**
```cpp
double total = static_cast<double>(t) + 4.0;
profile[i][j] /= total;
```
- `t` имеет тип `size_t` (unsigned integer)
- `static_cast<double>(t)` — явное преобразование в `double` перед арифметикой
- `+ 4.0` — 4 псевдонаблюдения (по одному на A, C, G, T)
- `/=` — составной оператор деления: `profile[i][j] = profile[i][j] / total`
- Результат: `profile[i][j] = (count + 1) / (t + 4)` — всегда > 0

**3. C-style массив для подсчёта в `score()`**
```cpp
int count[4] = {0, 0, 0, 0};
```
- Фиксированный массив из 4 элементов на стеке (не в куче)
- `= {0, 0, 0, 0}` — aggregate initialization: все элементы = 0
- Можно короче: `int count[4] = {};` — zero initialization (все нули)
- Быстрее, чем `vector<int>(4, 0)`: нет аллокации, нет деструктора
- Подходит, когда размер известен на этапе компиляции и мал

**4. Декомпозиция на функции**
```cpp
auto profile = buildProfileWithPseudocounts(motifs);
motifs.push_back(profileMostProbable(dna[j], k, profile));
```
- Каждый шаг алгоритма — отдельная функция:
  - `buildProfileWithPseudocounts` — матрица → профиль
  - `profileMostProbable` — текст + профиль → лучший k-мер
  - `score` — мотивы → числовая оценка
- Функции переиспользуются в BA2F (RandomizedMotifSearch) и BA2G (GibbsSampler)
- Принцип **single responsibility**: каждая функция решает одну задачу

## Запуск

```bash
# Rosalind-style аргументы
./greedy_motif_search --output answer.txt 3 5 \
  "GGCGTTCAGGCA" "AAGAATCAGTCA" "CAAGGAGTTCGC" "CACGTCAATCAC" "CAATAATATTCG"

# Входной файл в формате Rosalind
./greedy_motif_search --input rosalind_ba2e.txt --output rosalind_ba2e.out

# Короткая форма: один путь трактуется как входной файл
./greedy_motif_search rosalind_ba2e.txt
```

## Сложность

### Время: O(n × t × n × k)
- n - k + 1 стартовых k-меров из первой строки
- Для каждого: t-1 итераций построения профиля O(t×k) и поиска O(n×k)

### Память: O(t × k + 4 × k)
- Набор мотивов + профильная матрица
