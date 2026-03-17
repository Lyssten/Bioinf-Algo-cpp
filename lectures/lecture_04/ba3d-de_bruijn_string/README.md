# De Bruijn Graph of a String (BA3D)

## Задача

Построить граф де Брёйна строки Text для заданного k.

**Вход:** целое число `k` и строка `Text`
**Выход:** граф де Брёйна DeBruijn_k(Text) в виде списка смежности

Формат входного файла Rosalind:
```text
k
Text
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `k Text`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `de_bruijn_string_output.txt`

**Пример:**
```
k: 4
Text: AAGATTCTCTAC
Вывод:
AAG -> AGA
AGA -> GAT
ATT -> TTC
CTA -> TAC
CTC -> TCT
GAT -> ATT
TCT -> CTA,CTC
TTC -> TCT
```

## Теория

### Граф де Брёйна строки

Граф де Брёйна DeBruijn_k(Text) строится следующим образом:
1. Извлекаем все k-меры из строки Text
2. Для каждого k-мера создаём ребро из его **префикса** (первые k-1 символов) в его **суффикс** (последние k-1 символов)
3. Одинаковые узлы объединяются (merge)

В отличие от overlap-графа, где узлы — это k-меры, в графе де Брёйна узлы — это (k-1)-меры, а рёбра соответствуют k-мерам.

**Биологический смысл:**
- Граф де Брёйна — основа современных геномных ассемблеров (Velvet, SPAdes, MEGAHIT)
- Позволяет эффективно работать с большим количеством коротких ридов (reads)
- Задача сборки генома сводится к поиску эйлерова пути в графе де Брёйна

### Связь с эйлеровыми путями

Ключевое свойство: если строка Text может быть реконструирована из своих k-меров, то в графе де Брёйна DeBruijn_k(Text) существует эйлеров путь, проходящий через каждое ребро ровно один раз.

## Алгоритм

### Псевдокод
```
DeBruijnGraph(Text, k):
    graph ← пустой словарь (node → list of targets)
    for i from 0 to |Text| - k:
        kmer ← Text[i..i+k-1]
        prefix ← kmer[0..k-2]    // первые k-1 символов
        suffix ← kmer[1..k-1]    // последние k-1 символов
        graph[prefix].append(suffix)
    return graph
```

## Реализация на C++

```cpp
using AdjList = std::map<std::string, std::vector<std::string>>;

AdjList solve(const InputData& input) {
    AdjList graph;

    for (size_t i = 0; i + input.k <= input.text.size(); ++i) {  // (1)
        std::string kmer = input.text.substr(i, input.k);
        std::string prefix = kmer.substr(0, input.k - 1);         // (2)
        std::string suffix = kmer.substr(1, input.k - 1);         // (3)
        graph[prefix].push_back(suffix);                           // (4)
    }

    return graph;
}
```

### Особенности реализации

**1. `i + input.k <= input.text.size()` — безопасная граница цикла**
```cpp
for (size_t i = 0; i + input.k <= input.text.size(); ++i) {
```
- Условие `i + k <= size` эквивалентно `i <= size - k`, но безопасно для беззнаковых типов
- При `size < k` вычитание `size - k` дало бы переполнение `size_t` (обход через максимум)
- Сложение `i + k` всегда корректно, т.к. оба значения неотрицательны и не превышают `size`

**2. `substr(0, k-1)` — префикс k-мера**
```cpp
std::string prefix = kmer.substr(0, input.k - 1);
```
- `substr(pos, len)` возвращает подстроку длины `len`, начиная с позиции `pos`
- Префикс — первые k-1 символов k-мера
- Возвращает новый `std::string` (аллокация в куче для длинных строк)

**3. `substr(1, k-1)` — суффикс k-мера**
```cpp
std::string suffix = kmer.substr(1, input.k - 1);
```
- Суффикс — последние k-1 символов, т.е. начиная с позиции 1

**4. `std::map<std::string, std::vector<std::string>>` — упорядоченный граф**
```cpp
graph[prefix].push_back(suffix);
```
- `std::map` хранит ключи в отсортированном порядке (лексикографически)
- `operator[]` создаёт пустой `vector<string>` при первом обращении к новому ключу
- `push_back` добавляет суффикс в список смежности
- Порядок суффиксов в каждом списке соответствует порядку их появления в строке
- В отличие от `unordered_map`, `map` гарантирует сортировку при итерации — это нужно для формата вывода Rosalind

**5. Structured bindings при форматировании**
```cpp
for (const auto& [node, targets] : graph) {
```
- C++17 structured bindings: `[node, targets]` распаковывает `std::pair<const std::string, std::vector<std::string>>`
- `node` — ключ (имя узла), `targets` — значение (список смежных узлов)
- `const auto&` — без копирования пары

## Запуск

```bash
# Rosalind-style аргументы
./ba3d-de_bruijn_string --output answer.txt 4 "AAGATTCTCTAC"

# Входной файл в формате Rosalind
./ba3d-de_bruijn_string --input rosalind_ba3d.txt --output rosalind_ba3d.out

# Короткая форма: один путь трактуется как входной файл
./ba3d-de_bruijn_string rosalind_ba3d.txt
```

## Сложность

### Время: O(n * k * log n)
- n = |Text| - k + 1 — количество k-меров
- Для каждого k-мера: два `substr` за O(k) и вставка в `map` за O(k * log n)
- Итого: O(n * k * log n)

### Память: O(n * k)
- Граф хранит до n рёбер, каждое ребро — две строки длины k-1
- Количество уникальных узлов <= 2n
