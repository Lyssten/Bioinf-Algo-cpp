# Overlap Graph (BA3C)

## Задача

Построить граф перекрытий для коллекции k-меров.

**Вход:** коллекция k-меров `Patterns` (по одному на строку)
**Выход:** список смежности графа перекрытий (отсортированный по исходной вершине)

Ребро ведёт от k-мера A к k-меру B, если **суффикс** A (последние k-1 символов) совпадает с **префиксом** B (первые k-1 символов). Допускаются петли (ребро из вершины в себя).

Формат входного файла Rosalind:
```text
ATGCG
GCATG
CATGC
AGGCA
GGCAT
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `kmer1 kmer2 ...`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `overlap_graph_output.txt`

**Пример:**
```
Вход: ATGCG  GCATG  CATGC  AGGCA  GGCAT
Вывод:
AGGCA -> GGCAT
CATGC -> ATGCG
GCATG -> CATGC
GGCAT -> GCATG
```

## Теория

### Граф перекрытий

Граф перекрытий (overlap graph) — ориентированный граф, в котором:
- **вершины** — k-меры из входной коллекции
- **рёбра** — направленные связи от k-мера A к k-меру B, если `Suffix(A) = Prefix(B)`

Где:
- `Suffix(A)` — строка A без первого символа: `A[1..k-1]`
- `Prefix(B)` — строка B без последнего символа: `B[0..k-2]`

**Биологический смысл:**
- Граф перекрытий — базовая структура данных для сборки генома из коротких чтений (ридов)
- Ребро A → B означает, что рид B может идти сразу после рида A в геноме
- Гамильтонов путь в графе перекрытий соответствует реконструкции исходной последовательности
- На практике используют граф де Брёйна (более эффективный), но граф перекрытий — важная концептуальная основа

### Связь с задачей сборки генома

При секвенировании геном разбивается на короткие фрагменты (k-меры). Задача сборки — восстановить исходную последовательность. Граф перекрытий показывает, какие фрагменты могут следовать друг за другом, образуя перекрывающиеся участки длиной k-1.

## Алгоритм

### Псевдокод
```
OverlapGraph(Patterns):
    adjacency ← пустой словарь (отсортированный)
    k ← длина первого k-мера
    for each kmer_a in Patterns:
        suffix ← kmer_a[1..k-1]
        for each kmer_b in Patterns:
            prefix ← kmer_b[0..k-2]
            if suffix == prefix:
                adjacency[kmer_a].append(kmer_b)
    return adjacency
```

## Реализация на C++

```cpp
std::map<std::string, std::vector<std::string>> solve(const std::vector<std::string>& kmers) {
    std::map<std::string, std::vector<std::string>> adjacency;

    if (kmers.empty())
        return adjacency;

    size_t k = kmers[0].size();

    for (const auto& a : kmers) {
        std::string suffix = a.substr(1);              // (1)
        for (const auto& b : kmers) {
            std::string prefix = b.substr(0, k - 1);   // (2)
            if (suffix == prefix) {
                adjacency[a].push_back(b);              // (3)
            }
        }
    }

    return adjacency;
}
```

### Особенности реализации

**1. `a.substr(1)` — суффикс k-мера**
```cpp
std::string suffix = a.substr(1);
```
- `substr(pos)` без второго аргумента возвращает подстроку от позиции `pos` до конца
- Для k-мера длины k получаем строку длины k-1, начиная со второго символа
- Результат — новый `std::string`, аллокация в куче

**2. `b.substr(0, k - 1)` — префикс k-мера**
```cpp
std::string prefix = b.substr(0, k - 1);
```
- `substr(pos, len)` возвращает подстроку длины `len`, начиная с позиции `pos`
- Первые k-1 символов k-мера B

**3. `std::map` — отсортированный ассоциативный контейнер**
```cpp
std::map<std::string, std::vector<std::string>> adjacency;
adjacency[a].push_back(b);
```
- `std::map` хранит пары (ключ, значение), отсортированные по ключу (лексикографически для строк)
- `adjacency[a]` — при первом обращении автоматически создаёт пустой `vector<string>` (value-initialization)
- `push_back(b)` — добавляет целевую вершину в конец вектора смежности
- Итерация по `std::map` гарантирует порядок по ключу — вывод автоматически отсортирован
- В отличие от `std::unordered_map`, использует красно-чёрное дерево: O(log n) вставка и поиск

**4. Structured bindings в форматировании**
```cpp
for (const auto& [node, targets] : adjacency) {
```
- C++17 structured bindings: `[node, targets]` связывает первый и второй элемент пары
- `node` — `const std::string&` (ключ карты)
- `targets` — `const std::vector<std::string>&` (список смежности)
- Эквивалент C++11: `it->first`, `it->second`

**5. Обработка пустого входа**
```cpp
if (kmers.empty())
    return adjacency;
```
- Граничный случай: пустая коллекция → пустой граф
- Защита от `kmers[0].size()` при пустом векторе (undefined behavior)

## Запуск

```bash
# Rosalind-style аргументы
./ba3c-overlap_graph --output answer.txt \
  "ATGCG" "GCATG" "CATGC" "AGGCA" "GGCAT"

# Входной файл в формате Rosalind
./ba3c-overlap_graph --input rosalind_ba3c.txt --output rosalind_ba3c.out

# Короткая форма: один путь трактуется как входной файл
./ba3c-overlap_graph rosalind_ba3c.txt
```

## Сложность

### Время: O(n^2 * k)
- n = количество k-меров
- Для каждой пары (a, b) сравниваем суффикс и префикс длины k-1 → O(k)
- Всего n^2 пар → O(n^2 * k)

### Память: O(n^2 + n * k)
- O(n * k) — хранение входных k-меров
- O(n^2) — в худшем случае каждая вершина связана с каждой (полный граф)
