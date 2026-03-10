# Randomized Motif Search (BA2F)

## Задача

Найти набор мотивов рандомизированным алгоритмом.

**Вход:** целые числа `k`, `t` и набор строк `Dna`
**Выход:** набор мотивов, найденный за 1000 запусков `RandomizedMotifSearch`

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
- при запуске от аргументов или sample dataset без `--output` — в `randomized_motif_search_output.txt`

**Пример:**
```
k: 8, t: 5
Dna: CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
     GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
     ...
Вывод:
TCTCGGGG
CCAAGGTG
TACAGGCG
TTCAGGTG
TCCACGTG
```

## Теория

### Рандомизированные алгоритмы

В отличие от жадного алгоритма, RandomizedMotifSearch:
- Начинает со **случайного** набора мотивов
- Итеративно улучшает: строит профиль → находит лучшие мотивы → повторяет
- Останавливается, когда улучшение невозможно (локальный минимум)
- Запускается многократно (1000 раз) для поиска глобального минимума

**Преимущество:** не зависит от порядка строк (в отличие от жадного алгоритма)

### Алгоритм Монте-Карло

RandomizedMotifSearch — алгоритм Монте-Карло: не гарантирует оптимальный ответ, но с высокой вероятностью находит близкий к оптимальному при достаточном числе запусков.

## Алгоритм

### Псевдокод
```
RandomizedMotifSearch(Dna, k, t):
    случайно выбрать Motifs (по одному k-меру из каждой строки)
    BestMotifs ← Motifs
    while true:
        Profile ← ProfileWithPseudocounts(Motifs)
        for i ← 0 to t-1:
            Motifs[i] ← ProfileMostProbable(Dna[i], k, Profile)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs ← Motifs
        else:
            return BestMotifs

// Основной цикл
BestResult ← запуск с худшим возможным Score
repeat 1000 times:
    result ← RandomizedMotifSearch(Dna, k, t)
    if Score(result) < Score(BestResult):
        BestResult ← result
return BestResult
```

## Реализация на C++

```cpp
std::vector<std::string> randomizedMotifSearchOnce(
        const std::vector<std::string>& dna,
        size_t k, size_t t, std::mt19937& rng) {
    // Случайные начальные k-меры
    std::vector<std::string> motifs;
    for (size_t i = 0; i < t; ++i) {
        std::uniform_int_distribution<size_t> dist(0, dna[i].size() - k);
        motifs.push_back(dna[i].substr(dist(rng), k));
    }

    std::vector<std::string> bestMotifs = motifs;

    while (true) {
        auto profile = buildProfileWithPseudocounts(motifs);  // (1)
        for (size_t i = 0; i < t; ++i)
            motifs[i] = profileMostProbable(dna[i], k, profile);  // (2)

        if (score(motifs) < score(bestMotifs))  // (3)
            bestMotifs = motifs;
        else
            return bestMotifs;  // (4)
    }
}
```

### Особенности реализации

**1. `#include <random>` — генерация случайных чисел в C++11**
```cpp
#include <random>
std::mt19937 rng(42);  // генератор с seed = 42
```
- `std::mt19937` — генератор псевдослучайных чисел **Mersenne Twister**
- Период 2^19937 - 1 (≈ 10^6001) — достаточен для любых практических задач
- `rng(42)` — конструктор с seed (начальное состояние): одинаковый seed → одинаковая последовательность
- **Фиксированный seed = 42:** результаты воспроизводимы при каждом запуске
- Для непредсказуемого seed: `std::mt19937 rng(std::random_device{}());`
- **Не используйте `rand()`** из C: плохое качество, глобальное состояние, нет потокобезопасности

**2. `std::uniform_int_distribution` — равномерное распределение**
```cpp
std::uniform_int_distribution<size_t> dist(0, dna[i].size() - k);
size_t pos = dist(rng);
```
- `dist(a, b)` — создаёт распределение на [a, b] (включая оба конца!)
- `dist(rng)` — генерирует случайное число: принимает генератор как аргумент
- Разделение генератора и распределения: один `mt19937` + много `distribution`
- `dist` можно вызывать многократно — каждый вызов даёт новое число

**3. `while (true)` с `return` — цикл-с-выходом**
```cpp
while (true) {
    auto profile = buildProfileWithPseudocounts(motifs);
    for (size_t i = 0; i < t; ++i)
        motifs[i] = profileMostProbable(dna[i], k, profile);

    if (score(motifs) < score(bestMotifs))
        bestMotifs = motifs;
    else
        return bestMotifs;  // выход из функции
}
```
- `while (true)` — бесконечный цикл, выход только через `return` или `break`
- Паттерн «улучшение до сходимости»: повторяем, пока результат улучшается
- `return` сразу выходит из функции (не только из цикла, в отличие от `break`)
- Гарантия завершения: Score ∈ [0, t×k] и монотонно убывает → конечное число итераций

**4. Передача `std::mt19937&` по ссылке**
```cpp
std::vector<std::string> randomizedMotifSearchOnce(..., std::mt19937& rng)
```
- `rng` передаётся по **неконстантной ссылке** (`&`, не `const&`)
- Каждый вызов `dist(rng)` изменяет внутреннее состояние генератора
- Если передать по значению — каждый вызов функции получит копию с тем же состоянием → одинаковые «случайные» числа
- Один `rng` на все 1000 запусков обеспечивает различные начальные мотивы

## Запуск

```bash
# Rosalind-style аргументы
./randomized_motif_search --output answer.txt 8 5 \
  "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA" \
  "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG" \
  "TAGTACCGAGACCGAAAGAAGTATACAGGCGT" \
  "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC" \
  "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"

# Входной файл в формате Rosalind
./randomized_motif_search --input rosalind_ba2f.txt --output rosalind_ba2f.out

# Короткая форма: один путь трактуется как входной файл
./randomized_motif_search rosalind_ba2f.txt
```

## Сложность

### Время: O(1000 × I × t × n × k)
- 1000 запусков
- I — среднее число итераций до сходимости (обычно небольшое)
- На каждой итерации: построение профиля O(t×k) + поиск O(t×n×k)

### Память: O(t × k)
- Наборы мотивов + профильная матрица
