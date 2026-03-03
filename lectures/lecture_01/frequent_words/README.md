# Frequent Words

## Задача

Найти все наиболее часто встречающиеся k-меры в тексте.

**Вход:** строка `Text` и число `k`  
**Выход:** все k-меры с максимальной частотой

**Пример:**
```
Text: ACGTTGCATGTCGCATGATGCATGAGAGCT
k = 4
Вывод: ATGC CATG GCAT
```

## Теория

### k-мер
Подстрока длины `k`. В биоинформатике k-меры используются для:
- Поиска мотивов в ДНК
- Сборки генома (de Bruijn graphs)
- Выравнивания последовательностей

Для строки длины `n` количество k-меров равно `n - k + 1`.

## Алгоритм

### Идея
Однопроходный алгоритм с хэш-таблицей:
1. Извлекаем каждый k-мер из текста
2. Считаем частоту в хэш-таблице
3. Находим максимальную частоту
4. Собираем все k-меры с этой частотой

### Псевдокод
```
FrequentWords(Text, k):
    freq ← пустая хэш-таблица
    maxCount ← 0
    
    for i ← 0 to |Text| - k:
        pattern ← Text[i..i+k-1]
        freq[pattern] ← freq[pattern] + 1
        maxCount ← max(maxCount, freq[pattern])
    
    frequentPatterns ← пустое множество
    for each pattern in freq:
        if freq[pattern] = maxCount:
            добавить pattern в frequentPatterns
    
    return frequentPatterns
```

## Реализация на C++

```cpp
std::set<std::string> frequentWords(const std::string& text, size_t k) {
    if (k == 0 || k > text.size())
        return {};

    std::unordered_map<std::string, int> freq;  // (1)
    int maxCount = 0;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        int count = ++freq[pattern];  // (2)
        maxCount = std::max(maxCount, count);
    }

    std::set<std::string> result;  // (3)
    for (const auto& [pattern, count] : freq) {  // (4)
        if (count == maxCount)
            result.insert(pattern);
    }
    return result;
}
```

### Особенности реализации

**1. `std::unordered_map` — хэш-таблица**
```cpp
std::unordered_map<std::string, int> freq;
```
- Ключ = k-мер, значение = частота
- Вставка и доступ в среднем **O(1)**
- В худшем случае **O(n)** при коллизиях

**2. Инкремент в хэш-таблице**
```cpp
int count = ++freq[pattern];
```
- `freq[pattern]` создаёт запись со значением 0, если ключа нет
- Префиксный `++` сначала увеличивает, потом возвращает значение
- Возвращаемое значение сразу используется для обновления `maxCount`

**3. `std::set` для результата**
```cpp
std::set<std::string> result;
```
- Автоматически удаляет дубликаты
- Сортирует элементы лексикографически
- Основан на красно-чёрном дереве (**O(log n)** вставка)

**4. Structured binding (C++17)**
```cpp
for (const auto& [pattern, count] : freq)
```
- Деструктуризация пары ключ-значение
- Эквивалентно: `for (const auto& pair : freq)` с `pair.first`, `pair.second`

## Сравнение с наивным подходом

**Наивный (как в оригинале из Pattern Count):**
```cpp
// Для каждого k-мера вызываем pattern_count — O(n) каждый
for i ← 0 to |Text|-k:
    count ← pattern_count(Text, Text[i..i+k-1])  // O(n)
```
**Итого: O(n² × k)**

**Оптимизированный (наш):**
```cpp
// Один проход с хэш-таблицей
for i ← 0 to |Text|-k:
    freq[Text[i..i+k-1]]++  // O(1) в среднем
```
**Итого: O(n × k)**

## Запуск

```bash
# Без аргументов — тестовые данные
./frequent_words

# Аргументы: <Text> <k>
./frequent_words "ACGTTGCATGTCGCATGATGCATGAGAGCT" 4
# Вывод: ATGC CATG GCAT
```

## Сложность

### Время: O(n × k)

**Обоснование:**
- Внешний цикл: **n - k + 1** ≈ **O(n)** итераций (где n = |Text|)
- На каждой итерации:
  - `substr()` — **O(k)** копирование k-мера
  - Хэширование строки — **O(k)** (вычисление хэша по всем символам)
  - Вставка в `unordered_map` — **O(1)** в среднем (после хэширования)
  - `std::max()` — **O(1)**
- Второй цикл по хэш-таблице: **O(число_уникальных_k-меров)** ≤ **O(n)**
- Итого: **O(n × k)**

**Пример:** Text = 1000 символов, k = 4  
→ ~1000 × 4 = 4,000 операций

### Память: O(u × k)

**Обоснование:**
- `unordered_map` хранит уникальные k-меры
- u = количество уникальных k-меров (u ≤ n - k + 1)
- Каждый k-мер: **O(k)** байт
- `std::set` результата: **O(m × k)** где m = число частых k-меров
- Итого: **O(u × k)** ≤ **O(n × k)**

**Пример:** Text = 1000, k = 4, все k-меры уникальны  
→ ~1000 × 4 = 4,000 байт в хэш-таблице
