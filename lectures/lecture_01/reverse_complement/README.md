# Reverse Complement

## Задача

Построить обратный комплемент ДНК-последовательности.

**Вход:** строка `Text` (символы A, T, C, G)  
**Выход:** обратный комплемент

**Пример:**
```
Text: AAAACCCGGT
Вывод: ACCGGGTTTT
```

## Теория

### Комплементарность нуклеотидов
В ДНК нуклеотиды образуют пары:
- **A (аденин)** ↔ **T (тимин)**
- **C (цитозин)** ↔ **G (гуанин)**

### Обратный комплемент
1. Разворачиваем последовательность
2. Заменяем каждый нуклеотид на комплементарный

**Биологический смысл:** ДНК состоит из двух комплементарных цепей. Если одна цепь `5'-AAAACCCGGT-3'`, то вторая (читаемая в обратном направлении) будет `5'-ACCGGGTTTT-3'`.

## Алгоритм

### Псевдокод
```
ReverseComplement(Text):
    result ← пустая строка
    for i ← |Text|-1 downto 0:
        nucleotide ← Text[i]
        complement ← Complement(nucleotide)
        добавить complement к result
    return result

Complement(n):
    if n = 'A': return 'T'
    if n = 'T': return 'A'
    if n = 'C': return 'G'
    if n = 'G': return 'C'
```

## Реализация на C++

```cpp
std::string reverseComplement(const std::string& text) {
    std::string result;
    result.reserve(text.size());  // (1)

    for (auto it = text.rbegin(); it != text.rend(); ++it) {  // (2)
        switch (*it) {  // (3)
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default: result += *it;  // (4)
        }
    }
    return result;
}
```

### Особенности реализации

**1. `reserve()` — предварительное выделение памяти**
```cpp
result.reserve(text.size());
```
- Резервирует память под `text.size()` символов
- Избегает многократных реаллокаций при росте строки
- **Важно:** не изменяет `size()`, только `capacity()`

**2. Reverse iterators**
```cpp
for (auto it = text.rbegin(); it != text.rend(); ++it)
```
- `rbegin()` — итератор на последний элемент
- `rend()` — итератор на позицию перед первым
- Движение от конца к началу без явного разворота строки

**3. `switch` вместо `if-else`**
```cpp
switch (*it) { ... }
```
- Компилятор может оптимизировать в таблицу переходов (jump table)
- Быстрее цепочки `if-else` при большом количестве вариантов

**4. Обработка неизвестных символов**
```cpp
default: result += *it;
```
- Символы кроме A, T, C, G копируются как есть
- Полезно для работы с ambiguous-кодами (например, N = любой нуклеотид)

## Альтернативы

**Таблица поиска (lookup table):**
```cpp
const char complement[256] = {
    ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C'
};
result += complement[static_cast<unsigned char>(*it)];
```

**`std::transform` с лямбдой:**
```cpp
std::string result = text;
std::reverse(result.begin(), result.end());
std::transform(result.begin(), result.end(), result.begin(),
    [](char c) {
        switch(c) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default: return c;
        }
    });
```

## Запуск

```bash
# Без аргументов — тестовые данные
./reverse_complement

# Аргументы: <Text>
./reverse_complement "AAAACCCGGT"
# Вывод: ACCGGGTTTT
```

## Сложность

### Время: O(n)

**Обоснование:**
- Один проход по строке: **n** итераций (где n = |Text|)
- На каждой итерации:
  - Разыменование итератора — **O(1)**
  - `switch` — **O(1)** (таблица переходов)
  - `result += char` — **O(1)** амортизированное (с геометрическим ростом capacity)
- Итого: **O(n)**

**Пример:** Text = 1000 нуклеотидов → ~1000 операций

### Память: O(n)

**Обоснование:**
- Результат — строка длиной n
- `reserve(n)` выделяет память сразу, без реаллокаций
- Итераторы — O(1)
- Итого: **O(n)** для результата

**Пример:** Text = 1000 символов → ~1000 байт результата, где n = |Text|
