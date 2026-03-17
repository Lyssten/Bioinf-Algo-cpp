# De Bruijn Graph of k-mers (BA3E)

## Задача

Построить граф де Брёйна по заданной коллекции k-меров.

**Вход:** коллекция k-меров (по одному на строку)
**Выход:** список смежности графа де Брёйна (отсортированный по узлам)

Каждый k-мер становится ребром из своего префикса (первые k-1 символов) в суффикс (последние k-1 символов). Дубликаты k-меров порождают дубликаты рёбер.

Формат входного файла Rosalind:
```text
kmer1
kmer2
...
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `kmer1 kmer2 ...`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `de_bruijn_kmers_output.txt`

**Пример:**
```
Вход:
GAGG
CAGG
GGGG
GGGA
CAGG
AGGG
GGAG

Выход:
AGG -> GGG
CAG -> AGG,AGG
GAG -> AGG
GGA -> GAG
GGG -> GGA,GGG
```

## Теория

### Граф де Брёйна по коллекции k-меров

Граф де Брёйна `DeBruijn(Patterns)` строится по коллекции k-меров следующим образом:
- **Узлы** — все уникальные (k-1)-меры, являющиеся префиксами или суффиксами k-меров из коллекции
- **Рёбра** — каждый k-мер `p` порождает направленное ребро из `Prefix(p)` в `Suffix(p)`

Где:
- `Prefix(p)` = первые k-1 символов k-мера
- `Suffix(p)` = последние k-1 символов k-мера

**Биологический смысл:**
- Секвенирование ДНК разбивает геном на короткие фрагменты (риды)
- Граф де Брёйна — ключевая структура данных для сборки генома из ридов
- В отличие от overlap-графа (где узлы — риды), здесь узлы — (k-1)-меры, а рёбра — k-меры
- Задача сборки генома сводится к поиску Эйлерова пути в графе де Брёйна

**Ключевое свойство:** дубликаты k-меров порождают **мультирёбра** (параллельные рёбра между одной парой узлов). Это важно, так как в реальном секвенировании одинаковые k-меры встречаются многократно.

## Алгоритм

### Псевдокод
```
DeBruijnFromKmers(Patterns):
    graph ← пустой словарь (узел → список целей)
    for each kmer in Patterns:
        prefix ← kmer[0..k-2]
        suffix ← kmer[1..k-1]
        graph[prefix].append(suffix)
    отсортировать узлы и списки целей
    return graph
```

## Реализация на C++

```cpp
AdjList solve(const InputData& input) {
    AdjList graph;                                          // (1)

    for (const auto& kmer : input.kmers) {                  // (2)
        std::string prefix = kmer.substr(0, kmer.size() - 1);  // (3)
        std::string suffix = kmer.substr(1);                    // (4)
        graph[prefix].push_back(suffix);                    // (5)
    }

    for (auto& [node, targets] : graph) {                   // (6)
        std::sort(targets.begin(), targets.end());
    }

    return graph;
}
```

### Особенности реализации

**1. `std::map<std::string, std::vector<std::string>>` — упорядоченный словарь**
```cpp
using AdjList = std::map<std::string, std::vector<std::string>>;
```
- `std::map` хранит пары ключ-значение, отсортированные по ключу
- Ключ — имя узла (`std::string`), значение — список целей (`std::vector<std::string>`)
- Автоматическая сортировка узлов при итерации — не нужно сортировать вручную
- `std::unordered_map` был бы быстрее (O(1) vs O(log n)), но не гарантирует порядок

**2. `const auto&` в range-based for**
```cpp
for (const auto& kmer : input.kmers) {
```
- `auto` выводит тип элемента: здесь `const std::string&`
- `const &` — не копируем строку, не модифицируем
- Без `&`: на каждой итерации создавалась бы копия строки

**3. `substr(0, kmer.size() - 1)` — префикс k-мера**
```cpp
std::string prefix = kmer.substr(0, kmer.size() - 1);
```
- `substr(pos, len)` возвращает подстроку длины `len`, начиная с позиции `pos`
- Для k-мера длины k: префикс — первые k-1 символов
- `kmer.size() - 1` безопасно, так как `validateInput` гарантирует длину >= 2

**4. `substr(1)` — суффикс k-мера**
```cpp
std::string suffix = kmer.substr(1);
```
- `substr(pos)` без второго аргумента — от позиции `pos` до конца строки
- Суффикс — последние k-1 символов (с позиции 1 до конца)

**5. `graph[prefix].push_back(suffix)` — добавление ребра**
```cpp
graph[prefix].push_back(suffix);
```
- `operator[]` в `std::map`: если ключ не существует — создаёт пустой `vector<string>`
- `.push_back(suffix)` — добавляет суффикс в конец вектора
- Дубликаты k-меров автоматически создают дубликаты рёбер (мультирёбра)

**6. Structured bindings (C++17)**
```cpp
for (auto& [node, targets] : graph) {
    std::sort(targets.begin(), targets.end());
}
```
- `auto& [node, targets]` распаковывает `std::pair<const std::string, std::vector<std::string>>&`
- `node` — ссылка на ключ (const), `targets` — ссылка на значение (mutable)
- `&` без `const` — нужен для модификации `targets` через `std::sort`
- Альтернатива без structured bindings: `it->first`, `it->second`

## Запуск

```bash
# Rosalind-style аргументы
./ba3e-de_bruijn_kmers --output answer.txt \
  GAGG CAGG GGGG GGGA CAGG AGGG GGAG

# Входной файл в формате Rosalind
./ba3e-de_bruijn_kmers --input rosalind_ba3e.txt --output rosalind_ba3e.out

# Короткая форма: один путь трактуется как входной файл
./ba3e-de_bruijn_kmers rosalind_ba3e.txt
```

## Сложность

### Время: O(n * k + n * log n)
- n = количество k-меров, k = длина k-мера
- `substr` для каждого k-мера — O(k), всего O(n * k)
- Вставка в `std::map` — O(log m) на каждую операцию, где m — число уникальных узлов
- Сортировка списков смежности — суммарно O(n * log n) в худшем случае

### Память: O(n * k)
- Хранение графа: до n рёбер, каждое ребро — два (k-1)-мера
- `std::map` — дополнительные O(m) на узлы дерева
