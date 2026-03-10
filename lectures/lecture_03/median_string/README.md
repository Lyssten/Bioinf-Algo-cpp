# Median String (BA2B)

## Задача

Найти k-мер, минимизирующий суммарное расстояние до коллекции ДНК-строк.

**Вход:** целое число `k` и набор строк `Dna`
**Выход:** k-мер `Pattern`, минимизирующий `d(Pattern, Dna)`

Где:
- `k` — длина искомого мотива
- `Dna` — коллекция ДНК-строк
- `Pattern` — кандидат на медианную строку длины `k`

Формат входного файла Rosalind:
```text
k
dna1
dna2
...
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `k dna1 dna2 ...`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `median_string_output.txt`

**Пример:**
```
k: 3
Dna: AAATTGACGCAT  GACGACCACGTT  CGTCAGCGCCTG  GCTGAGCACCGG  AGTACGGGACAG
Вывод: GAC
```

## Теория

### Медианная строка

**Median String** — k-мер, который лучше всего «объясняет» все строки коллекции. Он минимизирует суммарное расстояние:
```
Median = argmin d(Pattern, Dna)
         Pattern
```

**Связь с задачей мотивов:**
- Поиск медианной строки эквивалентен поиску наилучшего мотива
- Вместо перебора по строкам (как в Motif Enumeration), перебираем по паттернам
- Пространство поиска: 4^k паттернов (фиксировано) vs C(n,1)^t комбинаций мотивов

**Применение:** более эффективный способ найти консенсусный мотив, когда k невелико

## Алгоритм

### Идея
Полный перебор всех 4^k возможных k-меров. Для каждого вычисляем d(Pattern, Dna) и запоминаем минимум.

### Псевдокод
```
MedianString(Dna, k):
    bestDist ← ∞
    median ← ""
    for i ← 0 to 4^k - 1:
        pattern ← NumberToPattern(i, k)
        dist ← DistanceBetweenPatternAndStrings(pattern, Dna)
        if dist < bestDist:
            bestDist ← dist
            median ← pattern
    return median
```

## Реализация на C++

```cpp
std::string medianString(const std::vector<std::string>& dna, size_t k) {
    int minDist = INT_MAX;
    std::string median;

    int64_t total = 1;
    for (size_t i = 0; i < k; ++i) total *= 4;  // (1)

    for (int64_t i = 0; i < total; ++i) {
        std::string pattern = numberToPattern(i, k);  // (2)
        int dist = distanceBetweenPatternAndStrings(pattern, dna);  // (3)
        if (dist < minDist) {
            minDist = dist;
            median = pattern;
        }
    }
    return median;
}
```

### Особенности реализации

**1. Итеративное вычисление 4^k вместо `std::pow`**
```cpp
int64_t total = 1;
for (size_t i = 0; i < k; ++i) total *= 4;
```
- `std::pow(4, k)` возвращает `double` — возможна потеря точности при больших k
- Итеративное умножение: `1 → 4 → 16 → 64 → ...` — целочисленная арифметика, точна
- `int64_t` позволяет k до 31 (4^31 < 2^63)
- **Альтернатива:** `int64_t total = 1LL << (2 * k);` — битовый сдвиг, ещё быстрее

**2. Переиспользование функций (композиция)**
```cpp
std::string pattern = numberToPattern(i, k);
int dist = distanceBetweenPatternAndStrings(pattern, dna);
```
- `numberToPattern` (из BA1M) и `distanceBetweenPatternAndStrings` (из BA2H) — отдельные функции
- Код не дублируется, а переиспользуется — принцип DRY (Don't Repeat Yourself)
- В данном проекте каждое решение самодостаточно (копия функций в `main.cpp`), но в реальном проекте они были бы в общем `.h` файле

**3. Перебор всех k-меров через числовое кодирование**
```cpp
for (int64_t i = 0; i < total; ++i) {
    std::string pattern = numberToPattern(i, k);
```
- Вместо рекурсивной генерации всех комбинаций ACGT — простой цикл от 0 до 4^k-1
- `numberToPattern` декодирует число обратно в строку
- Лексикографический порядок: i=0 → "AAA...A", i=1 → "AAA...C", ..., i=4^k-1 → "TTT...T"
- Подход **enumeration by encoding** — типичен для задач с ограниченным алфавитом

**Примечание:** при наличии нескольких оптимальных ответов возвращается лексикографически наименьший (первый с минимальным расстоянием).

## Запуск

```bash
# Rosalind-style аргументы
./median_string --output answer.txt 3 \
  "AAATTGACGCAT" "GACGACCACGTT" "CGTCAGCGCCTG" "GCTGAGCACCGG" "AGTACGGGACAG"

# Входной файл в формате Rosalind
./median_string --input rosalind_ba2b.txt --output rosalind_ba2b.out

# Короткая форма: один путь трактуется как входной файл
./median_string rosalind_ba2b.txt
```

## Сложность

### Время: O(4^k × t × n × k)
- 4^k паттернов
- Для каждого: d(Pattern, Dna) = O(t × n × k)
- Практично при k ≤ 12-14

### Память: O(k)
- Хранится только текущий лучший паттерн
