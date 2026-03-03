# Frequency Array (BA1K)

## Задача

Построить массив частот всех k-меров в строке.

**Вход:** строка Text и целое число k
**Выход:** массив частот длины 4^k, где элемент с индексом i — количество вхождений i-го k-мера (в лексикографическом порядке) в Text

**Пример:**
```
Text: ACGCGGCTCTGAAA
k: 2
Вывод: 2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0
```

## Теория

### Массив частот

Массив частот — вектор длины 4^k, где индекс = PatternToNumber(k-мер), значение = число вхождений.

**Для k=2 (16 элементов):**
```
Индекс: 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
k-мер:  AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
Частота: 2  1  0  0  0  0  2  2  1  2  1  0  0  1  1  0
```

**Проверка для "ACGCGGCTCTGAAA":**
- AA: позиции 12,13 → 2 вхождения
- AC: позиция 0 → 1
- CG: позиции 1,4 → 2
- CT: позиции 7,9 → 2
- GA: позиция 11 → 1
- GC: позиции 3,6 → 2
- GG: позиция 5 → 1
- TC: позиция 8 → 1
- TG: позиция 10 → 1

**Применение в биоинформатике:**
- Эффективный подсчёт частот всех k-меров за один проход
- Основа для поиска наиболее частых k-меров
- Сравнение геномных сигнатур между организмами

## Алгоритм

### Идея
Проходим по тексту скользящим окном размера k. Для каждого k-мера вычисляем PatternToNumber и инкрементируем соответствующий элемент массива.

### Псевдокод
```
ComputeFrequencyArray(Text, k):
    FrequencyArray ← массив из 4^k нулей
    for i ← 0 to |Text| - k:
        pattern ← Text[i..i+k-1]
        j ← PatternToNumber(pattern)
        FrequencyArray[j] ← FrequencyArray[j] + 1
    return FrequencyArray
```

## Реализация на C++

```cpp
int64_t patternToNumber(const std::string& pattern) {
    int64_t number = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        number = number * 4 + nucleotideToNumber(pattern[i]);
    }
    return number;
}

std::vector<int> computeFrequencyArray(const std::string& text, size_t k) {
    size_t arraySize = static_cast<size_t>(std::pow(4, k));  // (1)
    std::vector<int> freqArray(arraySize, 0);  // (2)

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        int64_t index = patternToNumber(pattern);  // (3)
        ++freqArray[index];  // (4)
    }
    return freqArray;
}
```

### Особенности реализации

**1. `std::pow` и приведение к целому типу**
```cpp
#include <cmath>  // для std::pow
size_t arraySize = static_cast<size_t>(std::pow(4, k));
```
- `std::pow(4, k)` возвращает `double` (даже для целых аргументов!)
- `double` имеет 53 бита мантиссы → точен до 2^53 ≈ 9×10¹⁵
- Для k ≤ 26: 4^26 ≈ 4.5×10¹⁵ < 2^53, результат точен
- `static_cast<size_t>()` усекает дробную часть (double → unsigned integer)
- **Альтернатива без float:** `size_t arraySize = 1; for (...) arraySize *= 4;`

**2. Конструктор `std::vector` с инициализацией**
```cpp
std::vector<int> freqArray(arraySize, 0);
```
- `vector<T>(n, val)` — создаёт вектор из `n` элементов, каждый = `val`
- Все 4^k ячеек заполнены нулями
- Память выделяется одним блоком `new int[arraySize]` — непрерывный массив
- Доступ по индексу `freqArray[i]` — O(1), как к обычному массиву
- **Отличие от `std::array`:** `vector` размещается в куче (heap), `array<T, N>` — на стеке, но `N` должен быть constexpr

**3. Прямая адресация через PatternToNumber**
```cpp
int64_t index = patternToNumber(pattern);
++freqArray[index];
```
- Каждый k-мер имеет уникальный числовой индекс в [0, 4^k)
- Вместо хеш-таблицы (unordered_map) используем **массив прямой адресации**
- Преимущество: O(1) доступ без хеширования и коллизий
- Недостаток: фиксированный размер 4^k (для k=10 → ~1M, для k=15 → ~1G элементов)
- Подход называется **frequency array** или **counting sort** в общем случае

**4. Оператор `++` (prefix increment)**
```cpp
++freqArray[index];
```
- `++x` (prefix) инкрементирует и возвращает новое значение
- `x++` (postfix) возвращает старое значение, потом инкрементирует
- Для примитивов разница минимальна, но для итераторов prefix быстрее (нет копии)
- Здесь результат не используется, поэтому выбор — дело стиля

## Альтернативные реализации

**Без `std::pow` — вычисление 4^k сдвигом:**
```cpp
size_t arraySize = static_cast<size_t>(1) << (2 * k);  // 1 << 2k = 4^k
```
- `1 << n` = 2^n, поэтому `1 << 2k` = 2^(2k) = 4^k
- Работает только для оснований-степеней двойки
- Целочисленная операция — без потерь точности

## Запуск

```bash
# Без аргументов — тестовые данные
./frequency_array

# Аргументы: <Text> <k>
./frequency_array "ACGCGGCTCTGAAA" 2
# Вывод: 2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0
```

## Сложность

### Время: O((|Text| - k) × k)
- Внешний цикл: |Text| - k + 1 итераций
- patternToNumber: O(k) на каждую итерацию
- Доступ к массиву: O(1)

### Память: O(4^k)
- Массив частот: 4^k элементов
- Растёт экспоненциально с k (k=10 → ~1M, k=15 → ~1G)
