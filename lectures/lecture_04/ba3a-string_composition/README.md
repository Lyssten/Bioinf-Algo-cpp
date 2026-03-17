# String Composition (BA3A)

## Задача

Найти все k-меры строки Text, отсортированные лексикографически.

**Вход:** целое число `k` и строка `Text`
**Выход:** все k-меры строки `Text` в лексикографическом порядке

Формат входного файла Rosalind:
```text
k
Text
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `k Text`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `string_composition_output.txt`

**Пример:**
```
k: 5
Text: CAATCCAAC
Вывод:
AATCC
ATCCA
CAATC
CCAAC
TCCAA
```

## Теория

### k-мерная композиция строки

**k-мерная композиция** (Composition_k) строки Text — это коллекция всех подстрок длины k, извлечённых скользящим окном и отсортированных лексикографически:

```
Composition_k(Text) = sort({Text[i..i+k-1] : 0 ≤ i ≤ |Text| - k})
```

**Биологический смысл:**
- Секвенаторы нового поколения (NGS) разбивают геном на короткие фрагменты — **риды** (reads)
- k-мерная композиция моделирует идеальный случай: каждый возможный k-мер прочитан ровно один раз
- Задача сборки генома — восстановить исходную строку по её k-мерной композиции

### Скользящее окно

Из строки длины n извлекается ровно (n - k + 1) подстрок длины k. Каждая подстрока начинается в позиции i, где 0 <= i <= n - k.

## Алгоритм

### Псевдокод
```
StringComposition(k, Text):
    kmers ← пустой список
    for i ← 0 to |Text| - k:
        добавить Text[i..i+k-1] в kmers
    sort(kmers)
    return kmers
```

## Реализация на C++

```cpp
std::vector<std::string> solve(const InputData& input) {
    std::vector<std::string> kmers;

    for (size_t i = 0; i + input.k <= input.text.size(); ++i) {  // (1)
        kmers.push_back(input.text.substr(i, input.k));           // (2)
    }

    std::sort(kmers.begin(), kmers.end());                        // (3)
    return kmers;
}
```

### Особенности реализации

**1. `i + input.k <= input.text.size()` — безопасное условие цикла**
```cpp
for (size_t i = 0; i + input.k <= input.text.size(); ++i) {
```
- Все значения имеют тип `size_t` (беззнаковый)
- Условие `i + k <= size` эквивалентно `i <= size - k`, но безопаснее: если `size < k`, то `size - k` переполнит `size_t` (станет огромным положительным числом), а `i + k <= size` корректно вернёт `false`
- Количество итераций: `n - k + 1`, где n = `text.size()`

**2. `substr(i, k)` — извлечение подстроки**
```cpp
kmers.push_back(input.text.substr(i, input.k));
```
- `std::string::substr(pos, count)` — возвращает новую строку длины `count`, начиная с позиции `pos`
- Каждый вызов выделяет память в куче для новой строки
- `push_back` перемещает строку в вектор (move semantics в C++17: возврат из `substr` — prvalue, вызывается move-конструктор)

**3. `std::sort` — лексикографическая сортировка**
```cpp
std::sort(kmers.begin(), kmers.end());
```
- `std::sort` из `<algorithm>` — интросорт (O(n log n) в среднем и худшем)
- Для `std::string` по умолчанию используется `operator<`, который сравнивает строки лексикографически (посимвольно)
- Сравнение двух строк длины k — O(k), итого сортировка O(m * k * log m), где m = n - k + 1
- `begin()` / `end()` возвращают итераторы произвольного доступа — необходимое требование `std::sort`

**4. `std::vector<std::string>` — динамический массив строк**
```cpp
std::vector<std::string> kmers;
```
- Хранит строки последовательно в памяти (указатели на данные в куче)
- `push_back` — амортизированное O(1): при нехватке ёмкости вектор удваивает capacity
- Возврат вектора из функции — гарантированный move (copy elision / NRVO в C++17)

## Запуск

```bash
# Rosalind-style аргументы
./ba3a-string_composition --output answer.txt 5 "CAATCCAAC"

# Входной файл в формате Rosalind
./ba3a-string_composition --input rosalind_ba3a.txt --output rosalind_ba3a.out

# Короткая форма: один путь трактуется как входной файл
./ba3a-string_composition rosalind_ba3a.txt
```

## Сложность

### Время: O(m * k * log m)
- m = n - k + 1 — количество k-меров
- Извлечение каждого k-мера — O(k)
- Сортировка m строк длины k — O(m * k * log m)
- Итого: O(m * k * log m)

### Память: O(m * k)
- Хранение m строк длины k в векторе
