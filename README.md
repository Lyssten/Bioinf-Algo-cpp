# Bioinformatics Algorithms (C++)

Коллекция алгоритмов биоинформатики на C++, сгруппированная по занятиям электива.

## Структура

```
cpp/
├── CMakeLists.txt
├── README.md
├── CLAUDE.md
└── lectures/
    ├── lecture_01/           # Строковые алгоритмы
    │   ├── pattern_count/
    │   ├── reverse_complement/
    │   ├── frequent_words/
    │   ├── pattern_matching/
    │   ├── clump_finding/
    │   ├── minimum_skew/
    │   ├── hamming_distance/
    │   └── frequent_words_mismatches/
    ├── lecture_02/           # Кодирование и приближённое сопоставление
    │   ├── pattern_to_number/
    │   ├── number_to_pattern/
    │   ├── approximate_pattern_matching/
    │   ├── frequency_array/
    │   ├── d_neighborhood/
    │   └── frequent_words_mismatches_complements/
    └── lecture_03/           # Поиск мотивов
        ├── motif_enumeration/
        ├── median_string/
        ├── profile_most_probable/
        ├── distance_pattern_strings/
        ├── greedy_motif_search/
        ├── randomized_motif_search/
        └── gibbs_sampler/
```

## Сборка

Все задачи собираются одной командой из корня проекта:

```bash
# 1. Перейти в корень проекта
cd cpp

# 2. Создать директорию для сборки и собрать
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

После сборки исполняемые файлы появляются внутри `build/` с зеркальной структурой:
```
build/lectures/lecture_XX/<task>/<task>
```

Директория `build/` добавлена в `.gitignore` — после клонирования репозитория её нужно пересобрать.

## Запуск

Все команды выполняются из директории `build/`:

```bash
cd build

# Lecture 01
./lectures/lecture_01/pattern_count/pattern_count "GCGCG" "GCG"
./lectures/lecture_01/reverse_complement/reverse_complement "AAAACCCGGT"
./lectures/lecture_01/frequent_words/frequent_words "ACGTTGCATGTCGCATGATGCATGAGAGCT" 4
./lectures/lecture_01/pattern_matching/pattern_matching "ATAT" "GATATATGCATATACTT"
./lectures/lecture_01/clump_finding/clump_finding "<genome>" 5 75 4
./lectures/lecture_01/minimum_skew/minimum_skew "<genome>"
./lectures/lecture_01/hamming_distance/hamming_distance "GGGCCGTTGGT" "GGACCGTTGAC"
./lectures/lecture_01/frequent_words_mismatches/frequent_words_mismatches "<text>" 4 1

# Lecture 02
./lectures/lecture_02/pattern_to_number/pattern_to_number "AGT"
./lectures/lecture_02/number_to_pattern/number_to_pattern 45 4
./lectures/lecture_02/approximate_pattern_matching/approximate_pattern_matching "ATTCTGGA" "<text>" 3
./lectures/lecture_02/frequency_array/frequency_array "ACGCGGCTCTGAAA" 2
./lectures/lecture_02/d_neighborhood/d_neighborhood "ACG" 1
./lectures/lecture_02/frequent_words_mismatches_complements/frequent_words_mismatches_complements "<text>" 4 1

# Lecture 03
./lectures/lecture_03/motif_enumeration/motif_enumeration 3 1 "ATTTGGC" "TGCCTTA" "CGGTATC" "GAAAATT"
./lectures/lecture_03/median_string/median_string 3 "AAATTGACGCAT" "GACGACCACGTT" "CGTCAGCGCCTG"
./lectures/lecture_03/profile_most_probable/profile_most_probable "<text>" 5 "0.2 0.2 ..." "0.4 0.3 ..." "0.3 0.3 ..." "0.1 0.2 ..."
./lectures/lecture_03/distance_pattern_strings/distance_pattern_strings "AAA" "TTACCTTAAC GATATCTGTC ..."
./lectures/lecture_03/greedy_motif_search/greedy_motif_search 3 5 "dna1" "dna2" ...
./lectures/lecture_03/randomized_motif_search/randomized_motif_search 8 5 "dna1" "dna2" ...
./lectures/lecture_03/gibbs_sampler/gibbs_sampler 8 5 100 "dna1" "dna2" ...
```

Без аргументов каждая программа запускается с тестовыми данными из условия задачи.

## Добавление нового занятия

1. Создайте папку `lectures/lecture_XX/`
2. Добавьте задачи в подпапки
3. Для каждой задачи создайте `CMakeLists.txt`:
   ```cmake
   add_executable(<name> main.cpp)
   ```
4. Добавьте `add_subdirectory` в корневой `CMakeLists.txt`

## Занятия

### Lecture 01: Строковые алгоритмы в биоинформатике

| Задача | Описание |
|--------|----------|
| [Pattern Count](lectures/lecture_01/pattern_count/README.md) | Подсчёт всех вхождений паттерна |
| [Reverse Complement](lectures/lecture_01/reverse_complement/README.md) | Обратный комплемент ДНК |
| [Frequent Words](lectures/lecture_01/frequent_words/README.md) | Частые k-меры |
| [Pattern Matching](lectures/lecture_01/pattern_matching/README.md) | Позиции паттерна в геноме |
| [Clump Finding](lectures/lecture_01/clump_finding/README.md) | (L,t)-клумбы k-меров |
| [Minimum Skew](lectures/lecture_01/minimum_skew/README.md) | Минимум skew (#G − #C) |
| [Hamming Distance](lectures/lecture_01/hamming_distance/README.md) | Расстояние Хэмминга |
| [Frequent Words with Mismatches](lectures/lecture_01/frequent_words_mismatches/README.md) | Частые k-меры с mismatches |

### Lecture 02: Кодирование и приближённое сопоставление

| Задача | Описание |
|--------|----------|
| [Pattern To Number](lectures/lecture_02/pattern_to_number/README.md) | ДНК-строка → число (base-4) |
| [Number To Pattern](lectures/lecture_02/number_to_pattern/README.md) | Число → ДНК-строка |
| [Approximate Pattern Matching](lectures/lecture_02/approximate_pattern_matching/README.md) | Поиск с ≤d несовпадениями |
| [Frequency Array](lectures/lecture_02/frequency_array/README.md) | Массив частот всех k-меров |
| [d-Neighborhood](lectures/lecture_02/d_neighborhood/README.md) | d-окрестность строки |
| [Frequent Words with Mismatches and Complements](lectures/lecture_02/frequent_words_mismatches_complements/README.md) | Частые k-меры + reverse complement |

### Lecture 03: Поиск мотивов

| Задача | Описание |
|--------|----------|
| [Motif Enumeration](lectures/lecture_03/motif_enumeration/README.md) | Перебор (k,d)-мотивов |
| [Median String](lectures/lecture_03/median_string/README.md) | k-мер с мин. расстоянием до Dna |
| [Profile-Most Probable](lectures/lecture_03/profile_most_probable/README.md) | k-мер с макс. вероятностью по профилю |
| [Distance Between Pattern and Strings](lectures/lecture_03/distance_pattern_strings/README.md) | d(Pattern, Dna) |
| [Greedy Motif Search](lectures/lecture_03/greedy_motif_search/README.md) | Жадный поиск с псевдосчётчиками |
| [Randomized Motif Search](lectures/lecture_03/randomized_motif_search/README.md) | Рандомизированный поиск мотивов |
| [Gibbs Sampler](lectures/lecture_03/gibbs_sampler/README.md) | Семплирование Гиббса |
