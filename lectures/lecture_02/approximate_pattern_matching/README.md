# Approximate Pattern Matching (BA1H)

## Задача

Найти все позиции приближённого вхождения паттерна в строку.

**Вход:** строки Pattern и Text, целое число d
**Выход:** все позиции в Text, где Pattern встречается с ≤ d несовпадениями

**Пример:**
```
Pattern: ATTCTGGA
Text: CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
d: 3
Вывод: 6 7 26 27 78
```

## Теория

### Приближённое сопоставление

k-мер Pattern **приближённо** совпадает с подстрокой Text, если расстояние Хэмминга между ними не превышает d:
```
HammingDistance(Pattern, Text[i..i+k-1]) ≤ d
```

**Применение в биоинформатике:**
- Поиск мотивов с мутациями в ДНК
- Выявление регуляторных сайтов (транскрипционные факторы связываются с вариациями одного мотива)
- Поиск точки начала репликации (oriC)

### Отличие от точного совпадения
- **Точное:** HammingDistance = 0
- **Приближённое:** HammingDistance ≤ d (допускаем d мутаций)

## Алгоритм

### Идея
Скользящее окно размером |Pattern| по тексту. Для каждой позиции проверяем, не превышает ли расстояние Хэмминга порога d.

### Псевдокод
```
ApproximatePatternMatching(Text, Pattern, d):
    positions ← пустой список
    for i ← 0 to |Text| - |Pattern|:
        if HammingDistance(Text[i..i+|Pattern|-1], Pattern) ≤ d:
            добавить i в positions
    return positions
```

## Реализация на C++

```cpp
int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

std::vector<int> approximatePatternMatching(
        const std::string& text, const std::string& pattern, size_t d) {
    std::vector<int> positions;

    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {  // (1)
        std::string window = text.substr(i, pattern.size());  // (2)
        if (hammingDistance(window, pattern) <= static_cast<int>(d))  // (3)
            positions.push_back(static_cast<int>(i));
    }
    return positions;
}
```

### Особенности реализации

**1. Границы скользящего окна**
```cpp
for (size_t i = 0; i <= text.size() - pattern.size(); ++i)
```
- Последняя допустимая позиция: `|Text| - |Pattern|`
- **Осторожно:** если `pattern.size() > text.size()`, вычитание `size_t` даёт underflow → огромное число
- Поэтому в функции есть проверка `if (pattern.size() > text.size()) return positions;`

**2. `std::string::substr()` — создание копии подстроки**
```cpp
std::string window = text.substr(i, pattern.size());
```
- `substr(pos, len)` возвращает **новую строку** — аллокация памяти O(k)
- На каждой итерации создаётся временная строка → k байт аллокации
- Для коротких строк (< 15-22 символа) работает SSO (Small String Optimization) — хранение на стеке без `new`
- Альтернатива без копирования: `std::string_view` (C++17) — см. ниже

**3. `std::vector<int>` — динамический массив**
```cpp
std::vector<int> positions;
...
positions.push_back(static_cast<int>(i));
```
- `vector` — массив переменного размера, хранит элементы в непрерывной памяти
- `push_back()` добавляет элемент в конец — **амортизированно O(1)**
- Стратегия роста: при нехватке места выделяет в 2 раза больше (capacity doubles)
- Пример: добавление 100 элементов → ~7 реаллокаций (1→2→4→8→16→32→64→128)

**4. `static_cast<int>()` — безопасное приведение типов**
```cpp
if (hammingDistance(window, pattern) <= static_cast<int>(d))
    positions.push_back(static_cast<int>(i));
```
- `static_cast<T>(expr)` — явное приведение типа в C++, проверяемое компилятором
- Зачем: `hammingDistance` возвращает `int`, а `d` и `i` имеют тип `size_t` (unsigned)
- Сравнение `signed <= unsigned` вызывает предупреждение компилятора (implicit conversion)
- C-style cast `(int)d` тоже работает, но `static_cast` нагляднее и безопаснее
- **Правило:** `static_cast` для числовых конверсий, `reinterpret_cast` для указателей, `dynamic_cast` для полиморфизма

## Альтернативные реализации

**С ранним выходом из hammingDistance:**
```cpp
int hammingDistanceWithThreshold(const std::string& s1, const std::string& s2, int maxDist) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i] && ++dist > maxDist)
            return dist;  // ранний выход
    }
    return dist;
}
```
- `++dist > maxDist` — сначала инкремент (prefix), потом сравнение (одно выражение вместо двух)
- Позволяет не считать полное расстояние, если порог уже превышен

**С `std::string_view` (C++17) — без копирования подстроки:**
```cpp
#include <string_view>

std::vector<int> approximatePatternMatching(
        std::string_view text, std::string_view pattern, size_t d) {
    std::vector<int> positions;
    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {
        auto window = text.substr(i, pattern.size());  // O(1)! Без аллокации
        if (hammingDistance(window, pattern) <= d)
            positions.push_back(i);
    }
    return positions;
}
```
- `string_view` — невладеющая «обёртка» над строкой: хранит только указатель + длину
- `string_view::substr()` возвращает новый `string_view` за O(1) — нет копирования
- Существенное ускорение при большом количестве итераций

## Запуск

```bash
# Без аргументов — тестовые данные
./approximate_pattern_matching

# Аргументы: <Pattern> <Text> <d>
./approximate_pattern_matching "ATTCTGGA" "CGCCCGAAT..." 3
# Вывод: 6 7 26 27 78
```

## Сложность

### Время: O((|Text| - |Pattern|) × |Pattern|)
- Внешний цикл: |Text| - |Pattern| + 1 позиций
- На каждой позиции: hammingDistance — O(|Pattern|)
- Итого: O(n × k), где n = |Text|, k = |Pattern|

### Память: O(m + k)
- Вектор позиций: O(m), где m — количество совпадений
- Подстрока окна: O(k)
