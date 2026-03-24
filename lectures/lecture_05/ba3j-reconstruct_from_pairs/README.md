# Reconstruct a String from its Paired Composition (BA3J)

## Задача

Реконструировать строку по её парной композиции (paired composition).

**Вход:** целые числа `k` и `d`, коллекция (k,d)-меров (парных чтений)
**Выход:** строка `Text`, парная композиция которой совпадает с данной коллекцией

Формат входного файла Rosalind:
```text
k d
read1_first|read1_second
read2_first|read2_second
...
```

(k,d)-мер в позиции i строки Text --- это пара (Text[i..i+k-1], Text[i+k+d..i+2k+d-1]). Расстояние d --- это разрыв (gap) между концом первого чтения и началом второго.

**Пример:**
```
Вход:
4 2
GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT

Вывод: GTGGTCGTGAGATGTTGA
```

## Теория

### От одиночных чтений к парным

В задачах BA3D--BA3F мы строили граф де Брёйна из обычных k-меров и искали эйлеров путь. Проблема: для реальных геномов с повторами длиннее k граф содержит множество эйлеровых путей, и мы не можем определить, какой из них соответствует настоящему геному.

Парные чтения (paired reads) решают эту проблему. Секвенатор читает два фрагмента ДНК на известном расстоянии друг от друга:
```
Геном:   ...GTGGTCGTGAGATGTTGA...
              ├─────┤  ├─────┤
              TGGT      TGAG
              ←──k──→  ←──k──→
                    ←d→
```

Пара (TGGT, TGAG) --- это (k,d)-мер при k=4, d=2. Параметр d = **расстояние между концом первого и началом второго** чтения (длина вставки минус 2k).

### Парный граф де Брёйна

Обычный граф де Брёйна: k-мер `AAGA` --- ребро из `AAG` в `AGA` (строковые узлы).

Парный граф де Брёйна: (k,d)-мер `(AAGA, TTGC)` --- ребро из парного узла `(AAG|TTG)` в `(AGA|TGC)`. Узлы --- это пары (k-1)-меров, записанные через разделитель `|`.

| | Обычный де Брёйн | Парный де Брёйн |
|---|---|---|
| **Узел** | (k-1)-мер | пара (k-1)-меров: `prefix1\|prefix2` |
| **Ребро** | k-мер | (k,d)-мер: пара k-меров |
| **Реконструкция** | эйлеров путь | эйлеров путь + верификация перекрытия |

Парные узлы несут **больше информации**: два фрагмента генома на расстоянии d+k дополнительно ограничивают возможные обходы графа. Если повтор короче d+k, парные чтения его «перешагивают» --- разрешая неоднозначность.

### Верификация перекрытия

После нахождения эйлерова пути мы «собираем» две строки: firstStr (из первых компонент узлов) и secondStr (из вторых). Строка secondStr начинается на позиции k+d результата. Хвост firstStr и начало secondStr должны совпасть:

```
Результат: GTGGTCGTGAGATGTTGA
firstStr:  GTGGTCGTGAGAT.....   (длина = |path| + k - 2)
secondStr: .......GAGATGTTGA    (начинается с позиции k+d = 6)
                   ↑↑↑↑↑
                   overlap — должны совпасть
```

Если хоть один символ не совпадает --- парная композиция не порождена одной строкой.

## Применение

### В биологии

- **Illumina paired-end секвенирование.** Большинство современных Illumina-запусков --- paired-end: каждый фрагмент ДНК читается с обоих концов. Типичные параметры: длина чтения 150 bp, длина вставки (insert size) 300--500 bp, следовательно d = insert_size - 2 * read_length = 0--200 bp. Для генома *E. coli* (4.6 Мб) при покрытии 50x секвенатор HiSeq генерирует ~15 миллионов пар чтений. Парные чтения превращают граф де Брёйна с тысячами эйлеровых путей в граф с единственным (или немногими) решениями.

- **Разрешение повторов.** Геном *E. coli* содержит 7 копий rRNA-оперона (~5.5 Кб каждая). Одиночные чтения по 150 bp не различают копии --- граф де Брёйна содержит «пузырь», допускающий 7! = 5040 маршрутов. Парные чтения с insert size > 5.5 Кб перешагивают повтор, фиксируя правильный порядок. На практике **SPAdes** использует парные чтения для построения скаффолдов: контиги (BA3K) соединяются в более длинные последовательности, используя пары, один конец которых попадает в один контиг, а другой --- в соседний.

- **Mate-pair библиотеки.** Для крупных вставок (2--20 Кб) используется mate-pair протокол. Программа **ALLPATHS-LG** (Broad Institute, 2011) объединяет paired-end (180 bp) и mate-pair (3 Кб) данные: первые дают точную локальную сборку, вторые --- дальнюю связность. Для генома мыши (2.7 Гб) ALLPATHS-LG достигал N50 > 20 Мб --- почти хромосомный уровень.

### В других областях

- **Стерео-триангуляция.** Две камеры (как два чтения) фиксируют одну и ту же сцену с известным смещением (baseline, аналог insert size). По парам точек на двух изображениях восстанавливается 3D-структура. В автономных автомобилях (Tesla, Waymo) стереокамеры с baseline 12--30 см дают точность глубины ~1% на дистанции до 50 м.

- **Мульти-сенсорная фузия данных.** Радар + лидар с известным смещением на крыше автомобиля: каждый сенсор даёт «чтение» окружения, а расстояние между ними фиксировано (аналог d). Совмещение двух независимых наблюдений устраняет неоднозначности --- точно как парные чтения устраняют повторы в геноме.

## Алгоритм

### Идея

Алгоритм повторяет схему сборки генома из BA3D--BA3F, но на парном графе:

1. **Строим парный граф де Брёйна.** Каждый (k,d)-мер `(GAGA, TTGA)` порождает ребро из `GAG|TTG` в `AGA|TGA`.
2. **Ищем эйлеров путь.** Находим стартовую вершину (out-degree - in-degree = 1), добавляем фиктивное ребро от финиша к старту, ищем эйлеров цикл алгоритмом Хирхольцера (BA3F), затем разрезаем цикл в нужном месте.
3. **Собираем две строки.** Идём по узлам пути: из первых компонент собираем firstStr, из вторых --- secondStr.
4. **Проверяем и объединяем.** Хвост firstStr (начиная с позиции k+d) должен совпадать с началом secondStr. Если совпадение есть --- конкатенируем: result = firstStr + secondStr[overlap..].

Почему работает: каждый узел пути хранит пару (k-1)-меров, связанных расстоянием d. При «склейке» firstStr и secondStr в одну строку эта связь гарантирует корректное заполнение разрыва между парами.

### Псевдокод
```
ReconstructFromPairs(k, d, reads):
    graph ← пустой парный граф
    for each (a, b) in reads:
        prefixNode ← prefix(a) + "|" + prefix(b)
        suffixNode ← suffix(a) + "|" + suffix(b)
        graph[prefixNode].append(suffixNode)

    path ← EulerianPath(graph)            // BA3F с доработкой для пути

    firstStr ← first(path[0])
    secondStr ← second(path[0])
    for i from 1 to |path| - 1:
        firstStr += last_char(first(path[i]))
        secondStr += last_char(second(path[i]))

    overlap ← |firstStr| - k - d
    for i from 0 to overlap - 1:
        if firstStr[k + d + i] ≠ secondStr[i]:
            error "Mismatch"

    return firstStr + secondStr[overlap..]
```

## Реализация на C++

```cpp
Graph buildPairedDeBruijn(const std::vector<PairedRead>& reads) {
    Graph graph;
    for (const auto& read : reads) {
        std::string prefixNode = makeNode(                              // (1)
            read.first.substr(0, read.first.size() - 1),
            read.second.substr(0, read.second.size() - 1));
        std::string suffixNode = makeNode(                              // (2)
            read.first.substr(1),
            read.second.substr(1));
        graph[prefixNode].push_back(suffixNode);                        // (3)
    }
    return graph;
}

std::vector<std::string> eulerianPath(Graph& graph) {
    std::map<std::string, int> inDeg, outDeg;                           // (4)
    std::set<std::string> allNodes;

    for (auto& [v, neighbors] : graph) {
        allNodes.insert(v);
        outDeg[v] = static_cast<int>(neighbors.size());
        for (auto& u : neighbors) {
            inDeg[u]++;
            allNodes.insert(u);
        }
    }

    std::string startNode, endNode;                                     // (5)
    for (auto& v : allNodes) {
        int diff = outDeg[v] - inDeg[v];
        if (diff == 1) startNode = v;
        else if (diff == -1) endNode = v;
    }

    if (startNode.empty()) {                                            // (6)
        return eulerianCycle(graph, graph.begin()->first);
    }

    graph[endNode].push_back(startNode);                                // (7)
    auto cycle = eulerianCycle(graph, startNode);
    graph[endNode].pop_back();

    for (size_t i = 0; i + 1 < cycle.size(); ++i) {                    // (8)
        if (cycle[i] == endNode && cycle[i + 1] == startNode) {
            std::vector<std::string> path;
            for (size_t j = i + 1; j < cycle.size(); ++j)
                path.push_back(cycle[j]);
            for (size_t j = 1; j <= i; ++j)
                path.push_back(cycle[j]);
            return path;
        }
    }

    return cycle;
}

std::string solve(const InputData& input) {
    auto graph = buildPairedDeBruijn(input.reads);                      // (9)
    auto path = eulerianPath(graph);                                    // (10)

    std::string firstStr, secondStr;
    for (size_t i = 0; i < path.size(); ++i) {                         // (11)
        auto [a, b] = splitNode(path[i]);
        if (i == 0) {
            firstStr = a;
            secondStr = b;
        } else {
            firstStr += a.back();
            secondStr += b.back();
        }
    }

    size_t k = input.k;
    size_t d = input.d;
    size_t overlap = firstStr.size() - k - d;                           // (12)

    for (size_t i = 0; i < overlap; ++i) {                              // (13)
        if (firstStr[k + d + i] != secondStr[i]) {
            throw std::runtime_error("Mismatch at overlap position " + std::to_string(i));
        }
    }

    return firstStr + secondStr.substr(overlap);                        // (14)
}
```

## Разбор кода

**1. Построение парных узлов: префикс**
```cpp
std::string prefixNode = makeNode(
    read.first.substr(0, read.first.size() - 1),
    read.second.substr(0, read.second.size() - 1));
```
- `read.first` --- первое чтение из пары (например, `"GAGA"`). `substr(0, size - 1)` отсекает последний символ: `"GAGA"` -> `"GAG"`.
- `makeNode(a, b)` возвращает `a + "|" + b` --- строковое представление парного узла, например `"GAG|TTG"`.
- Каждый парный узел кодирует **два** (k-1)-мера, разделённых `|`. Это компактное представление: одна строка вместо структуры с двумя полями.

**2. Построение парных узлов: суффикс**
```cpp
std::string suffixNode = makeNode(
    read.first.substr(1),
    read.second.substr(1));
```
- `substr(1)` без второго аргумента --- с позиции 1 до конца: `"GAGA"` -> `"AGA"`.
- Итого ребро `"GAG|TTG"` -> `"AGA|TGA"` соответствует (k,d)-меру `(GAGA, TTGA)`.

**3. Добавление ребра в граф**
```cpp
graph[prefixNode].push_back(suffixNode);
```
- `Graph = std::map<std::string, std::vector<std::string>>`. Ключ --- парный узел (строка `"GAG|TTG"`), значение --- вектор целевых узлов.
- `operator[]` у `map`: если ключа нет --- создаёт запись с пустым вектором. Это автоматическая «склейка»: если два разных (k,d)-мера имеют одинаковый префикс, их суффиксы добавляются в один вектор.
- Стоимость: O(k * log n), где n --- число узлов (поиск в красно-чёрном дереве по строковому ключу).

**4. Подсчёт степеней**
```cpp
std::map<std::string, int> inDeg, outDeg;
std::set<std::string> allNodes;
```
- Нужны для нахождения стартовой (out > in) и финишной (in > out) вершин эйлерова пути.
- `allNodes` собирает все вершины: некоторые могут быть только в значениях (target), но не ключами графа.
- `static_cast<int>(neighbors.size())` --- явное приведение `size_t` -> `int`. Необходимо, потому что `outDeg` хранит `int`, а `size()` возвращает `size_t`.

**5. Поиск стартовой и финишной вершин**
```cpp
int diff = outDeg[v] - inDeg[v];
if (diff == 1) startNode = v;
else if (diff == -1) endNode = v;
```
- Для эйлерова пути (не цикла) ровно одна вершина с `out - in = 1` (старт) и ровно одна с `in - out = 1` (финиш). Все остальные --- сбалансированы.
- `outDeg[v]` и `inDeg[v]` --- обращение к `map`, если ключа нет --- вернёт 0 (default для `int`).

**6. Частный случай: граф уже эйлеров**
```cpp
if (startNode.empty()) {
    return eulerianCycle(graph, graph.begin()->first);
}
```
- Если `startNode` не найден --- все вершины сбалансированы, граф имеет эйлеров цикл.
- `graph.begin()->first` --- первая вершина в лексикографическом порядке (свойство `std::map`).

**7. Приведение к эйлерову циклу**
```cpp
graph[endNode].push_back(startNode);
auto cycle = eulerianCycle(graph, startNode);
graph[endNode].pop_back();
```
- Стандартный трюк: добавляем фиктивное ребро endNode -> startNode, делая граф эйлеровым. Ищем цикл, затем удаляем фиктивное ребро.
- `pop_back()` удаляет последний элемент вектора --- именно то ребро, которое мы добавили.
- Восстановление графа (pop_back) --- хорошая практика, даже если граф больше не используется.

**8. Разрезание цикла в путь**
```cpp
for (size_t i = 0; i + 1 < cycle.size(); ++i) {
    if (cycle[i] == endNode && cycle[i + 1] == startNode) {
```
- Ищем позицию фиктивного ребра в цикле и разрезаем: всё после разреза + всё до разреза.
- `i + 1 < cycle.size()` вместо `i < cycle.size() - 1` --- защита от переполнения `size_t` при пустом цикле.
- Результат: путь от startNode до endNode, проходящий по всем рёбрам ровно раз.

**9--10. Основная функция: построение графа и эйлеров путь**
```cpp
auto graph = buildPairedDeBruijn(input.reads);
auto path = eulerianPath(graph);
```
- `auto` выводит тип как `Graph` (для graph) и `std::vector<std::string>` (для path).
- Разделение на два вызова: сначала граф, потом обход --- чистая архитектура, легко тестировать каждый этап отдельно.

**11. Сборка двух строк из пути**
```cpp
auto [a, b] = splitNode(path[i]);
if (i == 0) { firstStr = a; secondStr = b; }
else { firstStr += a.back(); secondStr += b.back(); }
```
- `splitNode` разбивает `"GAG|TTG"` на `("GAG", "TTG")` по разделителю `|`. Возвращает `std::pair<std::string, std::string>`.
- C++17 structured binding: `auto [a, b]` распаковывает пару в две переменные.
- Первый узел копируется целиком, для последующих берётся только последний символ --- `back()` (O(1)).
- В результате `firstStr` = строка, «прочитанная» из первых компонент всех узлов, `secondStr` --- из вторых.

**12. Вычисление длины перекрытия**
```cpp
size_t overlap = firstStr.size() - k - d;
```
- Длина firstStr = |path| + k - 2 = n - k - d + 1 + k - 2 = n - d - 1 (где n --- длина результата).
- Длина secondStr аналогична.
- secondStr начинается с позиции k+d в результирующей строке. Перекрытие = |firstStr| - (k+d) символов.

**13. Верификация перекрытия**
```cpp
if (firstStr[k + d + i] != secondStr[i])
    throw std::runtime_error("Mismatch at overlap position " + std::to_string(i));
```
- Проверяем посимвольно: хвост firstStr (с позиции k+d) должен совпадать с началом secondStr.
- Если не совпадает --- входные данные не порождены одной строкой. `std::runtime_error` --- для ошибок времени выполнения (некорректные данные).
- `std::to_string(i)` --- конвертация числа в строку для диагностики.

**14. Конкатенация результата**
```cpp
return firstStr + secondStr.substr(overlap);
```
- `secondStr.substr(overlap)` --- часть secondStr, которая **не** перекрывается с firstStr.
- `operator+` для строк --- конкатенация с аллокацией новой строки.
- Результат: полная реконструированная строка длины n = |firstStr| + |secondStr| - overlap.

## Запуск

```bash
# Входной файл в формате Rosalind
./ba3j-reconstruct_from_pairs --input rosalind_ba3j.txt --output answer.txt

# Короткая форма
./ba3j-reconstruct_from_pairs rosalind_ba3j.txt

# Без аргументов — sample dataset
./ba3j-reconstruct_from_pairs
```

## Анализ сложности

### Время: O(n * k * log n)

Где n --- число (k,d)-меров (парных чтений).

**Обоснование.** Алгоритм состоит из трёх этапов:

1. **Построение парного графа де Брёйна** --- O(n * k * log n). Для каждого из n чтений: 4 вызова `substr` --- O(k), конкатенация в `makeNode` --- O(k), вставка в `map` --- O(k * log n).

2. **Эйлеров путь** --- O(n). Граф имеет O(n) рёбер. Алгоритм Хирхольцера обходит каждое ребро ровно раз, но сравнение строковых ключей добавляет O(k) на каждую операцию в `map`. Итого O(n * k) для обхода. Однако из-за `edgeIndex` (прямой доступ по индексу) внутренний цикл Хирхольцера --- O(1) на ребро, а O(k) тратится только при `edgeIndex.count(v)` и `graph[v].size()`.

3. **Сборка строк и верификация** --- O(n * k). Проход по пути (n+1 узлов), `splitNode` --- O(k) каждый.

**Узкое место в коде --- построение графа:**
```cpp
graph[prefixNode].push_back(suffixNode);  // O(k * log n) — поиск в map по строковому ключу
```
Каждый вызов `operator[]` у `std::map` выполняет O(log n) сравнений ключей, каждое сравнение строк длины 2(k-1)+1 --- O(k). Итого: O(n * k * log n).

С `std::unordered_map` было бы O(k) амортизированно на вставку (хэширование строки), но `map` обеспечивает детерминированный порядок обхода.

### Память: O(n * k)

**Обоснование.** Граф хранит до n рёбер, каждое ребро --- две строки длины 2(k-1)+1.

**Основные структуры данных:**
```cpp
Graph graph;                        // O(n) записей, ключи и значения — строки O(k)
std::vector<std::string> stack;     // O(n) в худшем случае, строки O(k)
std::vector<std::string> cycle;     // O(n+1) строк O(k)
std::string firstStr, secondStr;    // O(n + k) каждая
```
- Парные узлы --- строки вида `"GAG|TTG"`, длина 2(k-1)+1 = 2k-1.
- Для n чтений граф содержит не более n+1 уникальных узлов.
- Стек и цикл в алгоритме Хирхольцера хранят строки, а не целые числа --- это увеличивает расход памяти в k раз по сравнению с числовым графом (BA3F).
- Итого: O(n * k) на граф + O(n * k) на стек/цикл + O(n) на результирующую строку = **O(n * k)**.
