# Frequent Words with Mismatches

## Задача

Найти наиболее частые k-меры с не более чем d несовпадениями.

**Вход:** `Text`, `k`, `d`  
**Выход:** k-меры с максимальным Countd(Text, Pattern)

**Пример:**
```
Text: ACGTTGCATGTCGCATGATGCATGAGAGCT
k=4, d=1
Вывод: ATGC ATGT GATG
```

## Теория

### Мismatch (несовпадение)
Позиция i в k-мерах p и q — **несовпадение**, если p[i] ≠ q[i].

### Countd(Text, Pattern)
Общее количество вхождений Pattern в Text с **не более чем d несовпадениями**.

**Пример:**
```
Text: AACAAGCTGATAAACATTTAAAGAG
Pattern: AAAAA, d=1

Вхождения с ≤1 несовпадением:
AACAA — 1 несовпадение (позиция 2: C≠A)
ATAAA — 1 несовпадение (позиция 1: T≠A)
AAACA — 1 несовпадение (позиция 4: C≠A)
AAAGA — 1 несовпадение (позиция 3: G≠A)

Count₁(Text, AAAAA) = 4
```

### Важное отличие
Pattern **не обязан точно встречаться** в Text! AAAAA может быть самым частым 5-мером с d=1, даже если точно AAAAA нигде нет.

## Алгоритм

### Идея
1. Для каждого k-мера в тексте генерируем все k-меры на расстоянии ≤ d
2. Считаем частоты этих "соседей"
3. Находим максимум

### Псевдокод
```
FrequentWordsWithMismatches(Text, k, d):
    freq ← пустая хэш-таблица
    maxCount ← 0
    
    for i ← 0 to |Text| - k:
        pattern ← Text[i..i+k-1]
        neighbors ← все k-меры с ≤ d mismatches от pattern
        
        for each neighbor in neighbors:
            freq[neighbor] ← freq[neighbor] + 1
            maxCount ← max(maxCount, freq[neighbor])
    
    frequentPatterns ← пустое множество
    for each pattern in freq:
        if freq[pattern] = maxCount:
            добавить pattern в frequentPatterns
    
    return frequentPatterns
```

## Реализация на C++

### Функция расстояния Хэмминга
```cpp
int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}
```

### Генерация соседей
```cpp
void generateNeighbors(const std::string& pattern, size_t d, 
                       std::set<std::string>& result) {
    // Рекурсивная генерация всех k-меров с ≤ d mismatches
    auto generate = [&](auto&& self, size_t pos, std::string& current, 
                        size_t mismatches) -> void {
        if (mismatches > d) return;
        if (pos == current.size()) {
            result.insert(current);
            return;
        }
        
        char orig = pattern[pos];
        for (char c : {'A', 'C', 'G', 'T'}) {
            current[pos] = c;
            size_t newMismatches = mismatches + (c != orig ? 1 : 0);
            if (newMismatches <= d)
                self(self, pos + 1, current, newMismatches);
        }
        current[pos] = orig;
    };
    
    std::string mutablePattern = pattern;
    generate(generate, 0, mutablePattern, 0);
}
```

### Основная функция
```cpp
std::set<std::string> frequentWordsWithMismatches(const std::string& text, 
                                                   size_t k, size_t d) {
    std::unordered_map<std::string, int> freq;
    int maxCount = 0;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        std::set<std::string> neighbors;
        generateNeighbors(pattern, d, neighbors);  // (1)
        
        for (const auto& neighbor : neighbors) {  // (2)
            int count = ++freq[neighbor];
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

**1. Генерация соседей**
```cpp
generateNeighbors(pattern, d, neighbors);
```
- Рекурсивный перебор всех комбинаций
- Для k=4, d=1: 4 (сам паттерн) + 4×3 = 16 соседей
- Для k=4, d=2: значительно больше

**2. Рекурсивная лямбда (C++14)**
```cpp
auto generate = [&](auto&& self, ...) -> void {
    ...
    self(self, pos + 1, current, newMismatches);
};
```
- Лямбда вызывает саму себя через параметр `self`
- Нужно для рекурсии внутри лямбды

**3. Подсчёт частот соседей**
```cpp
for (const auto& neighbor : neighbors) {
    int count = ++freq[neighbor];
    maxCount = std::max(maxCount, count);
}
```
- Каждый сосед учитывается один раз на позицию
- `maxCount` обновляется "на лету"

## Количество соседей

Для k-мера длины k с ≤ d mismatches:
```
N(k, d) = Σ C(k,i) × 3^i  для i от 0 до d
```

| k | d | Количество соседей |
|---|---|-------------------|
| 4 | 0 | 1 |
| 4 | 1 | 1 + 12 = 13 |
| 4 | 2 | 1 + 12 + 54 = 67 |
| 9 | 2 | 1 + 27 + 324 = 352 |

## Запуск

```bash
# Без аргументов — тестовые данные
./frequent_words_mismatches

# Аргументы: <Text> <k> <d>
./frequent_words_mismatches "ACGTTGCATGTCGCATGATGCATGAGAGCT" 4 1
# Вывод: ATGC ATGT GATG
```

## Сложность

### Время: O(n² × k × N(k,d))

**Обоснование:**
- Внешний цикл: **O(n - k + 1)** = **O(n)** итераций (где n = |Text|)
- На каждой итерации:
  - `substr()` — **O(k)**
  - `generateNeighbors()` — генерация всех k-меров с ≤ d mismatches
    - Количество соседей: **N(k,d) = Σ C(k,i) × 3^i** для i от 0 до d
    - Для k=4, d=1: N(4,1) = 1 + 12 = 13 соседей
    - Для k=9, d=2: N(9,2) = 1 + 27 + 324 = 352 соседа
  - Вставка соседей в `freq`: **O(N(k,d) × k)** (хэширование каждого)
- Итого: **O(n × k × N(k,d))**

**Пример:** n = 1000, k = 4, d = 1  
→ 1000 × 4 × 13 = 52,000 операций

**Пример:** n = 1000, k = 9, d = 2  
→ 1000 × 9 × 352 = 3,168,000 операций

### Память: O(u × k)

**Обоснование:**
- `freq` — хэш-таблица всех возможных k-меров с mismatches
- u = число уникальных k-меров (может быть до 4^k в теории)
- На практике: **O(n × N(k,d) × k)**
- `neighbors` — множество соседей для одного паттерна: **O(N(k,d) × k)**
- Итого: **O(n × N(k,d) × k)**

**Пример:** k = 4, d = 1, все соседи уникальны  
→ 1000 × 13 × 4 = 52,000 байт

### Таблица N(k,d) — количество соседей

| k | d=0 | d=1 | d=2 | d=3 |
|---|-----|-----|-----|-----|
| 4 | 1 | 13 | 67 | 253 |
| 6 | 1 | 19 | 163 | 877 |
| 9 | 1 | 28 | 352 | 2800 |
