# Motif Enumeration (BA2A)

## Задача

Найти все (k,d)-мотивы в коллекции ДНК-строк.

**Вход:** целые числа `k`, `d` и набор строк `Dna`
**Выход:** все k-меры, встречающиеся в каждой строке из `Dna` с ≤ `d` несовпадениями

Где:
- `k` — длина искомого мотива
- `d` — максимально допустимое число несовпадений
- `Dna` — коллекция ДНК-строк

Формат входного файла Rosalind:
```text
k d
dna1
dna2
...
```

Программа поддерживает два режима:
- Rosalind-style аргументы: `k d dna1 dna2 ...`
- путь к входному файлу в формате Rosalind

Результат печатается в `stdout` и сохраняется в файл:
- при `--output <path>` — в указанный путь
- при запуске от входного файла без `--output` — в `<input>.out`
- при запуске от аргументов или sample dataset без `--output` — в `motif_enumeration_output.txt`

**Пример:**
```
k: 3, d: 1
Dna: ATTTGGC  TGCCTTA  CGGTATC  GAAAATT
Вывод: ATA ATT GTT TTT
```

## Теория

### (k,d)-мотив

k-мер Pattern называется **(k,d)-мотивом** коллекции Dna, если он встречается **в каждой строке** с не более чем d мутациями:
```
∀ s ∈ Dna : ∃ позиция i : HammingDistance(Pattern, s[i..i+k-1]) ≤ d
```

**Биологический смысл:**
- Регуляторные мотивы (DnaA-боксы, сайты связывания транскрипционных факторов) сохраняются во всех последовательностях, но с мутациями
- (k,d)-мотив моделирует консервативный сигнал с допустимыми вариациями

### Подход: перебор с порождением соседей

Brute-force стратегия:
1. Для каждого k-мера в первой строке генерируем все соседей (≤ d mismatches)
2. Для каждого кандидата проверяем, есть ли он в каждой из оставшихся строк с ≤ d mismatches

## Алгоритм

### Псевдокод
```
MotifEnumeration(Dna, k, d):
    motifs ← ∅
    for each k-mer pattern in Dna[0]:
        for each neighbor in Neighbors(pattern, d):
            if neighbor appears in every string of Dna with ≤ d mismatches:
                motifs ← motifs ∪ {neighbor}
    return motifs
```

## Реализация на C++

```cpp
std::set<std::string> motifEnumeration(const std::vector<std::string>& dna,
                                        size_t k, size_t d) {
    std::set<std::string> motifs;

    for (size_t i = 0; i <= dna[0].size() - k; ++i) {
        std::string kmer = dna[0].substr(i, k);
        std::set<std::string> neighbors;
        generateNeighbors(kmer, d, neighbors);  // (1)

        for (const auto& candidate : neighbors) {
            bool inAll = true;
            for (size_t j = 1; j < dna.size(); ++j) {  // (2)
                if (!appearsWithMismatches(dna[j], candidate, d)) {
                    inAll = false;
                    break;  // (3)
                }
            }
            if (inAll)
                motifs.insert(candidate);
        }
    }
    return motifs;
}
```

### Особенности реализации

**1. `const std::vector<std::string>&` — вектор строк по ссылке**
```cpp
std::set<std::string> motifEnumeration(const std::vector<std::string>& dna, ...)
```
- `std::vector<std::string>` — динамический массив строк
- Каждый элемент `dna[i]` — отдельный `std::string` (своя память в куче)
- `const &` — передаём ссылку, не копируем весь массив строк
- Доступ по индексу `dna[0]`, `dna[j]` — O(1)
- `.size()` — количество строк в коллекции

**2. `const auto&` в range-based for**
```cpp
for (const auto& candidate : neighbors) {
```
- `auto` выводит тип элемента контейнера: здесь `const std::string&`
- `const` — не модифицируем элемент
- `&` — берём по ссылке, без копирования строки
- Без `&`: на каждой итерации создавалась бы копия строки — лишняя аллокация

**3. Флаг `bool` и ранний `break`**
```cpp
bool inAll = true;
for (size_t j = 1; j < dna.size(); ++j) {
    if (!appearsWithMismatches(dna[j], candidate, d)) {
        inAll = false;
        break;
    }
}
if (inAll) motifs.insert(candidate);
```
- `break` немедленно выходит из **ближайшего** цикла `for`
- Как только кандидат не найден в одной строке — проверять остальные бессмысленно
- Паттерн «все удовлетворяют условию»: инициализируем `true`, ставим `false` при первом нарушении
- **Альтернатива C++11:** `std::all_of(dna.begin()+1, dna.end(), [&](const auto& s) { return appearsWithMismatches(s, candidate, d); })`

**4. Парсинг переменного числа аргументов**
```cpp
if (argc >= 3) {
    k = std::stoul(argv[1]);
    d = std::stoul(argv[2]);
    dna.clear();
    for (int i = 3; i < argc; ++i)
        dna.push_back(argv[i]);
}
```
- `argc` — количество аргументов (включая имя программы)
- `argv` — массив C-строк (`char*`), `argv[0]` = имя программы
- `dna.clear()` — очищает вектор (size → 0, но capacity сохраняется)
- `argv[i]` неявно конвертируется в `std::string` через конструктор `string(const char*)`

## Запуск

```bash
# Rosalind-style аргументы
./motif_enumeration --output answer.txt 3 1 \
  "ATTTGGC" "TGCCTTA" "CGGTATC" "GAAAATT"

# Входной файл в формате Rosalind
./motif_enumeration --input rosalind_ba2a.txt --output rosalind_ba2a.out

# Короткая форма: один путь трактуется как входной файл
./motif_enumeration rosalind_ba2a.txt
```

## Сложность

### Время: O(n × N(k,d) × t × n × k)
- n = длина строки, t = количество строк
- N(k,d) = количество соседей каждого k-мера
- Для каждого соседа проверяем все строки, в каждой — скользящее окно

### Память: O(N(k,d) + |motifs|)
