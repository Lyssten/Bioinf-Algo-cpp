# Frequent Words with Mismatches and Reverse Complements (BA1J)

## Задача

Найти наиболее частые k-меры с учётом mismatches и обратных комплементов.

**Вход:** строка Text, целые числа k и d
**Выход:** все k-меры Pattern, максимизирующие Count_d(Text, Pattern) + Count_d(Text, ReverseComplement(Pattern))

**Пример:**
```
Text: ACGTTGCATGTCGCATGATGCATGAGAGCT
k: 4, d: 1
Вывод: ATGT ACAT
```

## Теория

### Расширение задачи Frequent Words with Mismatches

В задаче BA1I мы искали наиболее частые k-меры с ≤ d мутациями. Здесь мы дополнительно учитываем **обратный комплемент**: вхождения Pattern и ReverseComplement(Pattern) складываются.

**Почему это важно:**
- ДНК двуцепочечная — мотив на одной цепи читается как обратный комплемент на другой
- DnaA-боксы (сайты связывания белка DnaA) встречаются на обеих цепях
- Суммарная частота прямого + обратного мотива — более точный индикатор функциональности

**Пример:**
```
Pattern = "ATGT", ReverseComplement = "ACAT"
Если ATGT встречается 3 раза с ≤1 mismatch,
а ACAT встречается 2 раза с ≤1 mismatch,
то суммарный счёт ATGT = 3 + 2 = 5
```

**Важно:** Pattern не обязан присутствовать в Text точно — достаточно приближённых вхождений.

## Алгоритм

### Идея
Для каждого k-мера в тексте генерируем соседей (≤ d mismatches) как для прямой, так и для обратной комплементарной строки. Подсчитываем суммарные частоты и находим максимум.

### Псевдокод
```
FrequentWordsWithMismatchesAndComplements(Text, k, d):
    freq ← пустой словарь
    maxCount ← 0

    for i ← 0 to |Text| - k:
        pattern ← Text[i..i+k-1]
        rc ← ReverseComplement(pattern)
        neighbors ← Neighbors(pattern, d) ∪ Neighbors(rc, d)
        for each neighbor in neighbors:
            freq[neighbor] ← freq[neighbor] + 1
            maxCount ← max(maxCount, freq[neighbor])

    result ← {p : freq[p] = maxCount}
    return result
```

## Реализация на C++

```cpp
std::string reverseComplement(const std::string& pattern) {
    std::string result;
    result.reserve(pattern.size());
    for (auto it = pattern.rbegin(); it != pattern.rend(); ++it) {
        result += complement(*it);  // A↔T, C↔G
    }
    return result;
}

std::set<std::string> frequentWordsWithMismatchesAndComplements(
        const std::string& text, size_t k, size_t d) {
    std::unordered_map<std::string, int> freq;
    int maxCount = 0;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        std::string rc = reverseComplement(pattern);  // (1)

        std::set<std::string> neighbors;
        generateNeighbors(pattern, d, neighbors);  // (2)
        generateNeighbors(rc, d, neighbors);  // (3)

        for (const auto& neighbor : neighbors) {
            int count = ++freq[neighbor];  // (4)
            maxCount = std::max(maxCount, count);
        }
    }

    std::set<std::string> result;
    for (const auto& [p, count] : freq) {
        if (count == maxCount)
            result.insert(p);
    }
    return result;
}
```

### Особенности реализации

**1. `std::string::reserve()` — предвыделение памяти**
```cpp
std::string result;
result.reserve(pattern.size());
```
- `reserve(n)` выделяет память под `n` символов без изменения длины строки
- Без `reserve`: каждый `+=` может вызвать реаллокацию (копирование всей строки)
- С `reserve`: `+=` просто записывает в заранее выделенный буфер — O(1) каждый
- `size()` остаётся 0, `capacity()` становится ≥ n

**2. Обратные итераторы `rbegin()` / `rend()`**
```cpp
for (auto it = pattern.rbegin(); it != pattern.rend(); ++it) {
    result += complement(*it);
}
```
- `rbegin()` указывает на последний элемент, `rend()` — перед первым
- `++it` на обратном итераторе двигает **к началу** строки
- `*it` разыменовывает итератор → `char`
- `auto` выводит тип `std::string::const_reverse_iterator` (длинный, но auto скрывает)
- Результат: строка читается задом наперёд, каждый символ заменяется комплементом

**3. `std::unordered_map` — хеш-таблица**
```cpp
std::unordered_map<std::string, int> freq;
int count = ++freq[neighbor];
```
- Хранит пары ключ-значение с хеш-доступом O(1) в среднем
- `freq[key]` — если ключа нет, **создаёт его** со значением по умолчанию (`int` → 0)
- `++freq[neighbor]` — создаёт запись 0 (если нет) → инкрементирует → возвращает новое значение
- Один оператор заменяет: `if (!freq.count(key)) freq[key] = 0; freq[key]++; int count = freq[key];`

**4. `unordered_map::count()` для проверки наличия ключа**
```cpp
if (combined.count(p)) continue;
```
- `count(key)` возвращает 0 или 1 (ключ есть/нет)
- Эффективнее, чем `find(key) != end()` — не создаёт итератор
- С C++20: `contains(key)` — ещё нагляднее

**5. Structured bindings (C++17)**
```cpp
for (const auto& [p, count] : freq) {
    ...
}
```
- Деструктурирующее связывание: `[p, count]` привязывается к `std::pair<const std::string, int>`
- `p` → `pair.first` (ключ), `count` → `pair.second` (значение)
- `const auto&` — без копирования пар
- До C++17: `for (const auto& pair : freq) { pair.first; pair.second; }`

### Отличие от BA1I (без комплементов)
В BA1I генерируются соседи только для `pattern`. Здесь алгоритм разделён на два шага:
1. Подсчёт `freq[P]` = Count_d(Text, P) для всех P
2. Суммирование `freq[P] + freq[RC(P)]` для нахождения максимума

## Запуск

```bash
# Без аргументов — тестовые данные
./frequent_words_mismatches_complements

# Аргументы: <Text> <k> <d>
./frequent_words_mismatches_complements "ACGTTGCATGTCGCATGATGCATGAGAGCT" 4 1
# Вывод: ATGT ACAT
```

## Сложность

### Время: O((|Text| - k) × k × N(k, d))
- Внешний цикл: |Text| - k + 1 позиций
- reverseComplement: O(k)
- generateNeighbors (дважды): O(k × N(k, d)) каждый
- N(k, d) = Σ C(k,i) × 3^i для i=0..d

### Память: O(|Text| × N(k, d))
- unordered_map: до O((|Text| - k) × N(k, d)) записей
- set соседей на итерацию: O(N(k, d))
