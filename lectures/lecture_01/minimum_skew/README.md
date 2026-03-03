# Minimum Skew

## Задача

Найти позиции, где skew генома минимален.

**Skew** = (#G) − (#C) в префиксе генома.

**Вход:** `Genome`  
**Выход:** позиции i, где Skew(Prefixᵢ) минимален

**Пример:**
```
Genome: CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG
Вывод: 53 97
```

## Теория

### Skew диаграмма
**Skew(Genome)** — разность между количеством G и C в префиксе:
```
Skew(Prefixᵢ) = count(G, Prefixᵢ) - count(C, Prefixᵢ)
```

**Свойства:**
- Skew(Prefix₀) = 0 (пустой префикс)
- При чтении G: skew увеличивается на 1
- При чтении C: skew уменьшается на 1
- A и T не влияют на skew

### Биологический смысл
Минимум skew часто соответствует **origin of replication** (oriC) — точке начала репликации бактериальной хромосомы.

### Пример расчёта
```
Text:   C A T G G G C A T
i:      0 1 2 3 4 5 6 7 8 9
Skew:   0 -1 -1 -1 0 1 2 1 1 1

Минимум: -1 на позициях 1, 2, 3
```

## Алгоритм

### Идея
Один проход по геному с поддержкой текущего skew и отслеживанием минимума.

### Псевдокод
```
MinimumSkew(Genome):
    skew ← 0
    minSkew ← 0
    positions ← [0]  # Skew(0) = 0
    
    for i ← 1 to |Genome|:
        if Genome[i] = 'G':
            skew ← skew + 1
        else if Genome[i] = 'C':
            skew ← skew - 1
        
        if skew < minSkew:
            minSkew ← skew
            positions ← [i]  # новый минимум
        else if skew = minSkew:
            добавить i к positions
    
    return positions
```

## Реализация на C++

```cpp
std::vector<int> minimumSkew(const std::string& genome) {
    std::vector<int> positions;
    int skew = 0, minSkew = 0;

    positions.push_back(0);  // (1) Skew(0) = 0

    for (size_t i = 0; i < genome.size(); ++i) {  // (2)
        if (genome[i] == 'G') ++skew;  // (3)
        else if (genome[i] == 'C') --skew;

        if (skew < minSkew) {  // (4)
            minSkew = skew;
            positions.clear();  // (5)
            positions.push_back(i + 1);  // (6)
        } else if (skew == minSkew) {
            positions.push_back(i + 1);
        }
    }
    return positions;
}
```

### Особенности реализации

**1. Начальная позиция 0**
```cpp
positions.push_back(0);
```
- Skew пустого префикса = 0
- Это может быть минимумом (если дальше skew только растёт)

**2. Индексация в цикле**
```cpp
for (size_t i = 0; i < genome.size(); ++i)
```
- i = индекс символа в строке (0-based)
- Обрабатываем genome[0], genome[1], ..., genome[n-1]

**3. Условные инкременты**
```cpp
if (genome[i] == 'G') ++skew;
else if (genome[i] == 'C') --skew;
```
- Префиксный `++skew` быстрее чем `skew++` (хотя компилятор оптимизирует)
- `else if` — взаимно исключающие условия

**4. Новый минимум**
```cpp
if (skew < minSkew) {
    minSkew = skew;
    positions.clear();
    positions.push_back(i + 1);
}
```
- Нашли новый минимум — очищаем предыдущие позиции
- Сохраняем только текущую

**5. 1-based индексация результата**
```cpp
positions.push_back(i + 1);
```
- После обработки genome[i] получаем Skew(Prefixᵢ₊₁)
- В биоинформатике позиции обычно 1-based

**6. Равенство минимуму**
```cpp
else if (skew == minSkew) {
    positions.push_back(i + 1);
}
```
- Может быть несколько позиций с одинаковым минимумом

## Запуск

```bash
# Без аргументов — тестовые данные
./minimum_skew

# Аргументы: <Genome>
./minimum_skew "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG"
# Вывод: 53 97
```

## Сложность

### Время: O(n)

**Обоснование:**
- Один проход по геному: **n** итераций (где n = |Genome|)
- На каждой итерации:
  - Доступ к символу `genome[i]` — **O(1)**
  - Сравнение символов — **O(1)**
  - Инкремент/декремент skew — **O(1)**
  - Сравнение с минимумом — **O(1)**
  - `push_back()` / `clear()` — **O(1)** амортизированное
- Итого: **O(n)**

**Пример:** Genome = 4,600,000 (E. coli) → ~4,600,000 операций

### Память: O(m)

**Обоснование:**
- `vector<int>` для позиций минимума: **O(m)**
- Переменные `skew`, `minSkew`: **O(1)**
- Итого: **O(m)** где m = число позиций с минимальным skew

**Пример:** 2 позиции минимума → 2 × 4 = 8 байт

**Худший случай:** skew постоянно равен минимуму → m = n, память **O(n)**

**Типичный случай:** m << n (несколько позиций оригина репликации)
