# d-Neighborhood (BA1N)

## Задача

Сгенерировать d-окрестность строки — множество всех k-меров с расстоянием Хэмминга ≤ d от данного паттерна.

**Вход:** строка Pattern и целое число d
**Выход:** все строки из d-окрестности Pattern

**Пример:**
```
Pattern: ACG
d: 1
Вывод:
CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG
```

## Теория

### d-окрестность

**Neighbors(Pattern, d)** — множество всех k-меров, отличающихся от Pattern не более чем в d позициях:
```
Neighbors(Pattern, d) = {Pattern' : HammingDistance(Pattern, Pattern') ≤ d}
```

**Пример для Pattern="ACG", d=1:**
```
Исходный: ACG (0 mismatches)

Замена позиции 0:  CCG, GCG, TCG  (1 mismatch)
Замена позиции 1:  AAG, AGG, ATG  (1 mismatch)
Замена позиции 2:  ACA, ACC, ACT  (1 mismatch)

Итого: 10 строк (1 + 3×3)
```

**Размер d-окрестности:**
```
|Neighbors(k, d)| = Σ C(k, i) × 3^i  для i = 0..d
```
- d=0: 1 (только сам паттерн)
- d=1, k=3: 1 + 9 = 10
- d=2, k=3: 1 + 9 + 27 = 37

**Применение в биоинформатике:**
- Поиск мотивов с мутациями
- Перечисление всех возможных вариантов после d точечных мутаций
- Основа для алгоритмов Frequent Words with Mismatches

## Алгоритм

### Идея
Рекурсивно проходим по позициям паттерна. На каждой позиции пробуем все 4 нуклеотида. Если замена увеличивает число несовпадений — учитываем это. Если число несовпадений превысило d — обрезаем ветвь рекурсии.

### Псевдокод
```
Neighbors(Pattern, d):
    result ← ∅
    Generate(Pattern, 0, Pattern, 0, d, result)
    return result

Generate(original, pos, current, mismatches, d, result):
    if mismatches > d: return
    if pos = |current|:
        добавить current в result
        return
    for each nucleotide c in {A, C, G, T}:
        current[pos] ← c
        newMismatches ← mismatches + (c ≠ original[pos] ? 1 : 0)
        Generate(original, pos+1, current, newMismatches, d, result)
```

## Реализация на C++

```cpp
std::set<std::string> neighbors(const std::string& pattern, size_t d) {
    std::set<std::string> result;

    auto generate = [&](auto&& self, size_t pos,
                        std::string& current, size_t mismatches) -> void {
        if (mismatches > d) return;  // (1)
        if (pos == current.size()) {  // (2)
            result.insert(current);
            return;
        }

        char orig = pattern[pos];
        for (char c : {'A', 'C', 'G', 'T'}) {  // (3)
            current[pos] = c;
            size_t newMismatches = mismatches + (c != orig ? 1 : 0);  // (4)
            if (newMismatches <= d)
                self(self, pos + 1, current, newMismatches);
        }
        current[pos] = orig;  // (5)
    };

    std::string mutablePattern = pattern;
    generate(generate, 0, mutablePattern, 0);
    return result;
}
```

### Особенности реализации

**1. Рекурсивная лямбда (C++14) — разбор синтаксиса**
```cpp
auto generate = [&](auto&& self, size_t pos,
                    std::string& current, size_t mismatches) -> void {
    ...
    self(self, pos + 1, current, newMismatches);  // рекурсивный вызов
};
generate(generate, 0, mutablePattern, 0);  // первый вызов
```
Разберём каждый элемент:
- **`auto generate =`** — компилятор выводит тип лямбды (каждая лямбда имеет уникальный, невыразимый тип)
- **`[&]`** — захват всех внешних переменных по ссылке (`d`, `pattern`, `result`). Лямбда «видит» контекст функции
- **`auto&& self`** — forwarding reference на саму лямбду. Нужен для рекурсии: лямбда не может вызвать себя по имени (тип ещё не определён), поэтому передаёт себя как параметр
- **`-> void`** — явное указание возвращаемого типа. Без него компилятор не может вывести тип из рекурсивного вызова
- **`generate(generate, 0, ...)`** — при вызове передаём лямбду саму в себя первым аргументом

**Альтернатива — обычная функция:** но тогда нет доступа к `d`, `pattern`, `result` без передачи их как параметров.

**2. `std::set<std::string>` — упорядоченное множество**
```cpp
std::set<std::string> result;
result.insert(current);
```
- Реализован как **красно-чёрное дерево** (самобалансирующееся BST)
- Гарантии: элементы уникальны и автоматически отсортированы
- `insert()` — O(k × log n), где k = длина строки (сравнение), n = размер множества
- Для строк: сравнение лексикографическое (`"AAG" < "ACA" < "ACG"`)
- **vs `std::unordered_set`**: хеш-таблица, O(k) амортизированно, но без сортировки

**3. Range-based for по initializer list**
```cpp
for (char c : {'A', 'C', 'G', 'T'})
```
- `{'A', 'C', 'G', 'T'}` — `std::initializer_list<char>`, создаётся на стеке
- Range-based for (C++11) эквивалентен: `for (auto it = list.begin(); it != list.end(); ++it)`
- Более читаемо, чем `const char nucs[] = "ACGT"; for (int i = 0; i < 4; ++i)`

**4. Тернарный оператор**
```cpp
size_t newMismatches = mismatches + (c != orig ? 1 : 0);
```
- `условие ? значение_если_true : значение_если_false`
- Компактная замена `if-else`: `(c != orig)` → если заменили символ, добавляем 1 к счётчику
- Скобки вокруг `(c != orig ? 1 : 0)` — для ясности, технически не обязательны

**5. Backtracking — восстановление состояния**
```cpp
current[pos] = c;          // пробуем символ c
self(self, pos + 1, ...);  // рекурсия
current[pos] = orig;       // восстанавливаем исходный символ
```
- Паттерн backtracking: изменяем → рекурсия → откат
- `current` передаётся по ссылке (`std::string&`), изменения видны на всех уровнях рекурсии
- Без восстановления: после возврата из рекурсии `current` содержал бы «мусор» от дочерних вызовов
- Экономит память: одна строка `current` вместо копии на каждом уровне

**6. Отсечение (pruning)**
```cpp
if (mismatches > d) return;
if (newMismatches <= d)
    self(self, pos + 1, current, newMismatches);
```
- Двойная проверка: на входе в функцию и перед рекурсивным вызовом
- Без pruning: генерируем все 4^k строк; с pruning: только Σ C(k,i)×3^i
- Для k=10, d=2: без pruning 1M строк → с pruning 991 строка (ускорение в 1000 раз)

## Запуск

```bash
# Без аргументов — тестовые данные
./d_neighborhood

# Аргументы: <Pattern> <d>
./d_neighborhood "ACG" 1
# Вывод: каждая строка на отдельной строке
```

## Сложность

### Время: O(k × |Neighbors(k, d)|)
- Генерация каждого соседа: O(k) (глубина рекурсии)
- Количество соседей: Σ C(k,i) × 3^i для i=0..d
- Вставка в set: O(k × log(|result|))

### Память: O(k × |Neighbors(k, d)|)
- Хранение всех соседей в множестве
- Стек рекурсии: O(k)
