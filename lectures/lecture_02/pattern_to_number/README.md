# Pattern To Number (BA1L)

## Задача

Преобразовать ДНК-строку в число.

**Вход:** строка Pattern из алфавита {A, C, G, T}
**Выход:** PatternToNumber(Pattern)

**Пример:**
```
Pattern: AGT
Вывод: 11
```

## Теория

### Кодирование k-меров числами

Каждый нуклеотид кодируется числом от 0 до 3:
```
A = 0, C = 1, G = 2, T = 3
```

ДНК-строка рассматривается как число в **четверичной системе счисления** (основание 4):
```
PatternToNumber("AGT") = 0 × 4² + 2 × 4¹ + 3 × 4⁰ = 0 + 8 + 3 = 11
```

**Применение в биоинформатике:**
- Быстрое хеширование k-меров для массивов частот
- Компактное хранение k-меров (число вместо строки)
- Основа для алгоритмов подсчёта частот (Frequency Array)

### Формула

```
PatternToNumber(Pattern) = Σ NucleotideToNumber(Pattern[i]) × 4^(k-1-i)
                           i=0..k-1
```

Или рекурсивно:
```
PatternToNumber(Pattern) = 4 × PatternToNumber(Prefix(Pattern)) + NucleotideToNumber(LastSymbol(Pattern))
```

## Алгоритм

### Идея
Итеративно обрабатываем строку слева направо, накапливая результат по схеме Горнера: `number = number * 4 + digit`.

### Псевдокод
```
PatternToNumber(Pattern):
    number ← 0
    for i ← 0 to |Pattern|-1:
        number ← number × 4 + NucleotideToNumber(Pattern[i])
    return number
```

## Реализация на C++

```cpp
int nucleotideToNumber(char c) {
    switch (c) {
        case 'A': return 0;  // (1)
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

int64_t patternToNumber(const std::string& pattern) {
    int64_t number = 0;  // (2)
    for (size_t i = 0; i < pattern.size(); ++i) {  // (3)
        number = number * 4 + nucleotideToNumber(pattern[i]);  // (4)
    }
    return number;
}
```

### Особенности реализации

**1. Оператор `switch`**
```cpp
switch (c) {
    case 'A': return 0;
    ...
    default: return -1;
}
```
- `switch` по `char` компилируется в **jump table** — O(1) переход по значению символа
- `'A'`, `'C'` и т.д. — символьные литералы (тип `char`, одинарные кавычки, в отличие от `"A"` — строковый литерал `const char*`)
- `default` обрабатывает некорректный ввод — защита от UB (undefined behavior)
- `return` внутри `case` заменяет `break` — управление сразу покидает функцию

**2. Тип `int64_t` из `<cstdint>`**
```cpp
#include <cstdint>  // гарантированные размеры целых типов
int64_t number = 0;
```
- `int` обычно 32 бита → максимум 2^31 - 1 ≈ 2.1×10⁹
- Для k > 15: 4^16 = 4.3×10⁹ > INT_MAX — переполнение!
- `int64_t` — ровно 64 бита, поддерживает k до 31 (4^31 ≈ 4.6×10¹⁸ < 2^63)
- Альтернатива: `long long` (≥ 64 бит, но не гарантирует ровно 64)

**3. Передача `const std::string&`**
```cpp
int64_t patternToNumber(const std::string& pattern)
```
- `&` — передача по ссылке: функция работает с оригинальной строкой, без копирования
- `const` — гарантия, что функция не изменит строку
- Без `&` при каждом вызове создавалась бы копия строки — O(k) лишней работы и памяти
- Идиома C++: входные параметры всегда передаём как `const T&`

**4. Схема Горнера**
```cpp
number = number * 4 + nucleotideToNumber(pattern[i]);
```
- Вычисляет полином a₀×4^(k-1) + a₁×4^(k-2) + ... + aₖ₋₁ без возведения в степень
- Эквивалентно: `((a₀ × 4 + a₁) × 4 + a₂) × 4 + ...`
- Один проход слева направо: k умножений вместо k умножений + k возведений в степень

**5. Доступ по индексу `pattern[i]`**
```cpp
pattern[i]  // vs pattern.at(i)
```
- `operator[]` — без проверки границ, быстрее
- `.at(i)` — с проверкой, бросает `std::out_of_range` при выходе за границы
- В алгоритмическом коде границы контролируются циклом, поэтому `[]` безопасен и предпочтителен

## Альтернативные реализации

**С использованием `std::unordered_map`:**
```cpp
int64_t patternToNumber(const std::string& pattern) {
    static const std::unordered_map<char, int> encoding = {
        {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
    };
    int64_t number = 0;
    for (char c : pattern)
        number = number * 4 + encoding.at(c);
    return number;
}
```
- `static const` — инициализируется один раз за время жизни программы
- `encoding.at(c)` — с проверкой наличия ключа
- Менее эффективно, чем `switch` (хеширование vs jump table)

**С использованием битовых операций:**
```cpp
int64_t patternToNumber(const std::string& pattern) {
    int64_t number = 0;
    for (char c : pattern) {
        number <<= 2;  // сдвиг на 2 бита = умножение на 4
        number |= (c >> 1) & 3;  // A=0, C=1, T=3, G=2 (не лексикографический!)
    }
    return number;
}
```
- `<<= 2` эквивалентно `*= 4`, но на уровне битов
- Работает только для ASCII-кодировки нуклеотидов

## Запуск

```bash
# Без аргументов — тестовые данные
./pattern_to_number

# Аргументы: <Pattern>
./pattern_to_number "AGT"
# Вывод: 11
```

## Сложность

### Время: O(k)
- Один проход по строке длины k
- На каждой итерации: умножение, сложение, маппинг — O(1)

### Память: O(1)
- Только переменная `number` и индекс цикла
- Входная строка по константной ссылке
