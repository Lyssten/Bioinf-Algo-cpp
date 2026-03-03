# Distance Between Pattern and Strings (BA2H)

## Задача

Вычислить суммарное расстояние между паттерном и коллекцией строк.

**Вход:** строка Pattern и набор строк Dna
**Выход:** d(Pattern, Dna) — сумма минимальных расстояний Хэмминга

**Пример:**
```
Pattern: AAA
Dna: TTACCTTAAC  GATATCTGTC  ACGGCGTTCG  CCCTAAAGAG  CGTCAGAGGT
Вывод: 5
```

## Теория

### d(Pattern, Text)

Минимальное расстояние Хэмминга между Pattern и всеми k-мерами в Text:
```
d(Pattern, Text) = min { HammingDistance(Pattern, Text[i..i+k-1]) }
                   i=0..|Text|-k
```

### d(Pattern, Dna)

Сумма d(Pattern, Text) по всем строкам коллекции:
```
d(Pattern, Dna) = Σ d(Pattern, Text_i)
```

**Применение:**
- Целевая функция для задачи Median String
- Мера качества мотива: чем меньше d(Pattern, Dna), тем лучше Pattern описывает консервативный сигнал

## Алгоритм

### Псевдокод
```
DistanceBetweenPatternAndStrings(Pattern, Dna):
    k ← |Pattern|
    distance ← 0
    for each string Text in Dna:
        minDist ← ∞
        for i ← 0 to |Text| - k:
            dist ← HammingDistance(Pattern, Text[i..i+k-1])
            minDist ← min(minDist, dist)
        distance ← distance + minDist
    return distance
```

## Реализация на C++

```cpp
int distanceBetweenPatternAndStrings(const std::string& pattern,
                                      const std::vector<std::string>& dna) {
    size_t k = pattern.size();
    int totalDist = 0;

    for (const auto& text : dna) {
        int minDist = INT_MAX;  // (1)
        for (size_t i = 0; i <= text.size() - k; ++i) {
            int dist = hammingDistance(pattern, text.substr(i, k));
            minDist = std::min(minDist, dist);  // (2)
        }
        totalDist += minDist;  // (3)
    }
    return totalDist;
}
```

### Особенности реализации

**1. `INT_MAX` из `<climits>` — сентинельное значение**
```cpp
#include <climits>
int minDist = INT_MAX;  // 2147483647 на 32-битном int
```
- `INT_MAX` — максимальное значение типа `int`, определено в `<climits>`
- Используется как «бесконечность»: любое реальное расстояние Хэмминга будет меньше
- Гарантирует, что первое найденное значение заменит начальное
- **Альтернативы:** `std::numeric_limits<int>::max()` из `<limits>` (более C++-стиль)

**2. `std::min()` из `<algorithm>`**
```cpp
minDist = std::min(minDist, dist);
```
- Возвращает меньшее из двух значений
- Шаблонная функция: `template<class T> const T& min(const T& a, const T& b)`
- Оба аргумента должны быть одного типа (иначе — ошибка компиляции)
- **Альтернатива:** `if (dist < minDist) minDist = dist;` — эквивалентно, но менее компактно

**3. Range-based for по вектору строк**
```cpp
for (const auto& text : dna) {
    ...
}
```
- `dna` — `const std::vector<std::string>&`
- `text` на каждой итерации — `const std::string&` (ссылка на элемент вектора)
- Нет копирования строк, нет индексной переменной
- Эквивалент: `for (size_t idx = 0; idx < dna.size(); ++idx) { const auto& text = dna[idx]; ... }`

**4. `std::istringstream` — парсинг строки как потока**
```cpp
#include <sstream>
std::istringstream iss(argv[2]);
std::string s;
while (iss >> s)
    dna.push_back(s);
```
- `istringstream` оборачивает строку в поток ввода (как `std::cin`)
- `iss >> s` — извлекает слова, разделённые пробелами
- `while (iss >> s)` — цикл до конца строки (оператор `>>` возвращает `false` при EOF)
- Позволяет передать `"AAA BBB CCC"` как один аргумент CLI и разобрать на отдельные строки

## Запуск

```bash
./distance_pattern_strings "AAA" "TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT"
# Вывод: 5
```

## Сложность

### Время: O(t × n × k)
- t строк, каждая длины n, для каждой — скользящее окно с Хэммингом O(k)

### Память: O(1)
- Только счётчики, входные данные по ссылке
