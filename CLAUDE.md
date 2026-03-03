# AI Agent Context — Bioinformatics Algorithms (C++)

## Обзор проекта

Учебный репозиторий с решениями задач по биоинформатике (электив). Задачи из [Rosalind](https://rosalind.info/), сгруппированы по лекциям.

## Структура

```
cpp/
├── CMakeLists.txt              # Корневой CMake (C++17, добавляет все subdirectory)
├── README.md                    # Документация проекта (русский)
├── CLAUDE.md                    # ← этот файл
├── .gitignore                   # build/, IDE, *.o, *.so, *.exe
└── lectures/
    ├── LECTURE_TEMPLATE.md       # Шаблон для новых лекций
    ├── lecture_01/               # Строковые алгоритмы (8 задач)
    ├── lecture_02/               # Кодирование и приближённое сопоставление (6 задач)
    └── lecture_03/               # Поиск мотивов (7 задач)
```

Каждая задача — отдельная директория:
```
task_name/
├── main.cpp          # Полное самодостаточное решение
├── CMakeLists.txt    # Одна строка: add_executable(task_name main.cpp)
└── README.md         # Описание на русском: теория, псевдокод, реализация, C++ синтаксис, сложность
```

## Соглашения по коду

### Стандарт: C++17
Используемые фичи: structured bindings, `auto` в лямбдах, `std::string_view`, initializer lists.

### Именование
- Директории: `snake_case` (например, `frequent_words_mismatches`)
- Функции: `camelCase` (например, `patternCount`, `hammingDistance`)
- Исполняемые файлы: совпадают с именем директории

### Типы данных
| Назначение | Тип |
|-----------|-----|
| Размеры, индексы циклов | `size_t` |
| Счётчики, расстояния | `int` |
| Большие числа (PatternToNumber) | `int64_t` |
| ДНК-строки | `const std::string&` (вход) / `std::string` (результат) |
| Частоты | `std::unordered_map<std::string, int>` |
| Уникальные результаты (отсортированные) | `std::set<std::string>` |
| Позиции | `std::vector<int>` |
| Профильная матрица | `std::vector<std::vector<double>>` |

### Паттерн main()
```cpp
int main(int argc, char* argv[]) {
    // Значения по умолчанию (sample dataset)
    std::string text = "SAMPLE";
    size_t k = 4;

    // Парсинг CLI
    if (argc == N) {
        text = argv[1];
        k = std::stoul(argv[2]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    // Вычисление
    auto result = solve(text, k);

    // Вывод (пробелы между элементами, endl в конце)
    for (size_t i = 0; i < result.size(); ++i) {
        if (i > 0) std::cout << " ";
        std::cout << result[i];
    }
    std::cout << std::endl;
    return 0;
}
```

### Общие утилиты (копируются в каждый main.cpp)
- `hammingDistance(s1, s2)` → `int`
- `generateNeighbors(pattern, d, result)` → заполняет `set<string>`
- `nucleotideToNumber(char)` / `numberToNucleotide(int)` — маппинг A↔0, C↔1, G↔2, T↔3
- `patternToNumber(pattern)` / `numberToPattern(index, k)` — кодирование k-меров числами
- `reverseComplement(pattern)` → `string`
- `buildProfileWithPseudocounts(motifs)` → `vector<vector<double>>`
- `profileMostProbable(text, k, profile)` → `string`
- `score(motifs)` → `int`

Нет общих заголовочных файлов. Каждый main.cpp самодостаточен.

## Задачи по лекциям

### Lecture 01: Строковые алгоритмы
| Директория | Rosalind | Суть |
|-----------|----------|------|
| pattern_count | BA1A | Подсчёт вхождений паттерна |
| reverse_complement | BA1C | Обратный комплемент ДНК |
| frequent_words | BA1B | Наиболее частые k-меры |
| pattern_matching | BA1D | Позиции паттерна в геноме |
| clump_finding | BA1E | (L,t)-клумбы k-меров |
| minimum_skew | BA1F | Минимум skew (G-C) |
| hamming_distance | BA1G | Расстояние Хэмминга |
| frequent_words_mismatches | BA1I | Частые k-меры с ≤d мутациями |

### Lecture 02: Кодирование и приближённое сопоставление
| Директория | Rosalind | Суть |
|-----------|----------|------|
| pattern_to_number | BA1L | ДНК → число (base-4) |
| number_to_pattern | BA1M | Число → ДНК (обратно) |
| approximate_pattern_matching | BA1H | Позиции с ≤d mismatches |
| frequency_array | BA1K | Массив частот 4^k |
| d_neighborhood | BA1N | d-окрестность строки |
| frequent_words_mismatches_complements | BA1J | Частые k-меры + reverse complement |

### Lecture 03: Поиск мотивов
| Директория | Rosalind | Суть |
|-----------|----------|------|
| motif_enumeration | BA2A | Brute-force (k,d)-мотивы |
| median_string | BA2B | Медианная строка (min distance) |
| profile_most_probable | BA2C | k-мер с макс. вероятностью по профилю |
| distance_pattern_strings | BA2H | d(Pattern, Dna) |
| greedy_motif_search | BA2E | Жадный поиск с псевдосчётчиками |
| randomized_motif_search | BA2F | Рандомизированный поиск (1000 запусков) |
| gibbs_sampler | BA2G | Семплирование Гиббса (20 рестартов) |

## Сборка

```bash
cd cpp && mkdir -p build && cd build && cmake .. && make -j$(nproc)
```

Бинарники: `build/lectures/lecture_XX/task_name/task_name`

## Добавление новой лекции

1. Создать `lectures/lecture_XX/task_name/` с `main.cpp`, `CMakeLists.txt`, `README.md`
2. Добавить `add_subdirectory(...)` в корневой `CMakeLists.txt`
3. Обновить корневой `README.md`
4. README на русском, с подробными пояснениями C++ синтаксиса (см. существующие README как образец)

## Стиль README

Каждый README содержит разделы:
1. **Задача** — вход/выход/пример
2. **Теория** — математика и биологический смысл
3. **Алгоритм** — идея + псевдокод
4. **Реализация на C++** — код с пронумерованными комментариями
5. **Особенности реализации** — подробные пояснения C++ синтаксиса (типы, контейнеры, паттерны)
6. **Альтернативные реализации** (опционально) — другие подходы на C++
7. **Запуск** — примеры CLI
8. **Сложность** — время и память с обоснованием
