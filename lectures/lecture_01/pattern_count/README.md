# Pattern Count

## Задача

Посчитать количество всех вхождений паттерна `Pattern` в текст `Text` (включая перекрывающиеся).

**Вход:** строки `Text` и `Pattern`  
**Выход:** число вхождений

**Пример:**
```
Text: GCGCG
Pattern: GCG
Вывод: 2  (позиции 0 и 2)
```

## Алгоритм

### Идея
Используем **скользящее окно** фиксированной длины `|Pattern|`. На каждой позиции сравниваем подстроку текста с паттерном. Если совпадает — увеличиваем счётчик.

### Псевдокод
```
PatternCount(Text, Pattern):
    count ← 0
    for i ← 0 to |Text| - |Pattern|:
        if Text[i..i+|Pattern|-1] = Pattern:
            count ← count + 1
    return count
```

### Почему перекрывающиеся вхождения учитываются?
После проверки позиции `i` сдвигаемся ровно на 1, а не на длину паттерна. Поэтому для `GCGCG` и `GCG`:
- i=0: `GCG`CG → совпадение
- i=1: G`CGC`G → нет
- i=2: GC`GCG` → совпадение

## Реализация на C++

```cpp
int patternCount(const std::string& text, const std::string& pattern) {
    if (pattern.empty() || pattern.size() > text.size())
        return 0;

    int count = 0;
    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {
        if (text.substr(i, pattern.size()) == pattern)
            ++count;
    }
    return count;
}
```

### Особенности реализации

**1. Проверка граничных условий**
```cpp
if (pattern.empty() || pattern.size() > text.size())
    return 0;
```
- Пустой паттерн не может иметь вхождений
- Паттерн длиннее текста — вхождений нет

**2. Тип `size_t` для индексов**
```cpp
for (size_t i = 0; ...
```
- `size_t` — беззнаковый тип, используется для индексов и размеров
- Предотвращает предупреждения компилятора при сравнении с `.size()`

**3. `substr()` для извлечения подстроки**
```cpp
text.substr(i, pattern.size())
```
- Создаёт новую строку длиной `pattern.size()` начиная с позиции `i`
- **O(k)** по времени и памяти, где k = длина паттерна

**4. Сравнение строк**
```cpp
... == pattern
```
- Оператор `==` для `std::string` сравнивает посимвольно
- **O(k)** по времени

## Альтернативы

**`std::string_view` (C++17)** — для избежания копирования:
```cpp
std::string_view view(text);
if (view.substr(i, pattern.size()) == pattern) ...
```

**`text.compare()`** — более явное сравнение:
```cpp
if (text.compare(i, pattern.size(), pattern) == 0) ...
```

## Запуск

```bash
# Без аргументов — тестовые данные
./pattern_count

# Аргументы: <Text> <Pattern>
./pattern_count "GCGCG" "GCG"
# Вывод: 2
```

## Сложность

### Время: O(|Text| × |Pattern|)

**Обоснование:**
- Внешний цикл: **O(|Text| - |Pattern| + 1)** = **O(|Text|)** итераций
- На каждой итерации:
  - `substr()` — **O(|Pattern|)** копирование символов
  - Сравнение строк — **O(|Pattern|)** посимвольное сравнение
- Итого: **O(|Text| × |Pattern|)**

**Пример:** Text = 1000 символов, Pattern = 10 символов  
→ ~1000 × 10 = 10,000 операций

### Память: O(|Pattern|)

**Обоснование:**
- `substr()` создаёт временную строку длиной |Pattern|
- Счётчик `count` — O(1)
- Итого: **O(|Pattern|)** дополнительной памяти

**Пример:** Pattern = 10 символов → ~10 байт временно
