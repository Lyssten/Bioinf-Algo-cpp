# Pattern Matching

## Задача

Найти все позиции вхождения паттерна в геном.

**Вход:** `Pattern` и `Genome`  
**Выход:** позиции (0-based) где Pattern встречается в Genome

**Пример:**
```
Pattern: ATAT
Genome: GATATATGCATATACTT
Вывод: 1 3 9
```

## Алгоритм

### Идея
Классический **алгоритм поиска подстроки** скользящим окном. Проверяем каждую позицию генома на совпадение с паттерном.

### Псевдокод
```
PatternMatching(Pattern, Genome):
    positions ← пустой список
    
    for i ← 0 to |Genome| - |Pattern|:
        if Genome[i..i+|Pattern|-1] = Pattern:
            добавить i к positions
    
    return positions
```

### Визуализация
```
Genome:  G A T A T A T G C A T A T A C T T
Index:   0 1 2 3 4 5 6 7 8 9 ...

Pattern: ATAT

i=1:  G[ATAT]ATGCATATACTT ✓ → positions = [1]
i=3:  GAT[ATAT]GCATATACTT ✓ → positions = [1, 3]
i=9:  GATATATGC[ATAT]ACTT ✓ → positions = [1, 3, 9]
```

## Реализация на C++

```cpp
std::vector<int> findPatternPositions(const std::string& genome, 
                                       const std::string& pattern) {
    std::vector<int> positions;  // (1)
    
    if (pattern.empty() || pattern.size() > genome.size())
        return positions;

    for (size_t i = 0; i <= genome.size() - pattern.size(); ++i) {  // (2)
        if (genome.substr(i, pattern.size()) == pattern)  // (3)
            positions.push_back(i);  // (4)
    }
    return positions;
}
```

### Особенности реализации

**1. `std::vector<int>` для позиций**
```cpp
std::vector<int> positions;
```
- Динамический массив целых чисел
- Автоматическое управление памятью
- **O(1)** амортизированное добавление в конец

**2. Границы цикла**
```cpp
i <= genome.size() - pattern.size()
```
- Последняя возможная позиция начала паттерна
- Если `genome.size() = 17`, `pattern.size() = 4`, то последняя позиция = 13
- Проверка `pattern.size() > genome.size()` предотвращает underflow для `size_t`

**3. `substr()` для сравнения**
```cpp
genome.substr(i, pattern.size()) == pattern
```
- Создаёт временную строку длиной `pattern.size()`
- Сравнивает с паттерном посимвольно
- **O(k)** где k = длина паттерна

**4. `push_back()` для добавления**
```cpp
positions.push_back(i);
```
- Добавляет элемент в конец вектора
- При превышении `capacity()` происходит реаллокация (обычно ×2)

## Вывод результата

```cpp
for (size_t i = 0; i < positions.size(); ++i) {
    if (i > 0) std::cout << " ";  // Пробел между числами
    std::cout << positions[i];
}
```

**Альтернатива с range-based for:**
```cpp
bool first = true;
for (int pos : positions) {
    if (!first) std::cout << " ";
    std::cout << pos;
    first = false;
}
```

## Продвинутые алгоритмы

Для больших геномов (миллионы оснований) наивный алгоритм медленный. Быстрые алгоритмы:

| Алгоритм | Время | Память |
|----------|-------|--------|
| Наивный | O(n×m) | O(1) |
| Кнута-Морриса-Пратта | O(n+m) | O(m) |
| Бойера-Мура | O(n/m) в среднем | O(σ) |
| Суффиксное дерево | O(n+m) | O(n) |

где n = |Genome|, m = |Pattern|, σ = размер алфавита

## Запуск

```bash
# Без аргументов — тестовые данные
./pattern_matching

# Аргументы: <Pattern> <Genome>
./pattern_matching "ATAT" "GATATATGCATATACTT"
# Вывод: 1 3 9
```

## Сложность

### Время: O(|Genome| × |Pattern|)

**Обоснование:**
- Внешний цикл: **O(|Genome| - |Pattern| + 1)** = **O(|Genome|)** итераций
- На каждой итерации:
  - `substr()` — **O(|Pattern|)** копирование
  - Сравнение строк — **O(|Pattern|)** посимвольное
  - `push_back()` — **O(1)** амортизированное
- В худшем случае (много совпадений): **O(|Genome| × |Pattern|)**

**Пример:** Genome = 1,000,000 (бактерия), Pattern = 10  
→ ~1,000,000 × 10 = 10,000,000 операций

**Лучший случай:** Pattern не встречается, сравнение обрывается рано → **O(|Genome|)**

### Память: O(m)

**Обоснование:**
- `vector<int>` для позиций: **O(m)** где m = число вхождений
- Временная строка от `substr()`: **O(|Pattern|)**
- Итого: **O(m + |Pattern|)**

**Пример:** Pattern встречается 100 раз → ~400 байт (int = 4 байта)

**Худший случай:** Pattern = "A", Genome = "AAAA..." → m = |Genome|, память **O(|Genome|)**
