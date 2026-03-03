# Number To Pattern (BA1M)

## Задача

Преобразовать число в ДНК-строку заданной длины.

**Вход:** целое число index и длина k
**Выход:** NumberToPattern(index, k)

**Пример:**
```
index: 45
k: 4
Вывод: AGTC
```

## Теория

### Обратное преобразование

NumberToPattern — обратная операция к PatternToNumber. Число переводится из десятичной в четверичную систему, цифры заменяются нуклеотидами:
```
0 = A, 1 = C, 2 = G, 3 = T
```

**Пример для index=45, k=4:**
```
45 / 4 = 11, остаток 1 → C  (младший разряд)
11 / 4 = 2,  остаток 3 → T
2  / 4 = 0,  остаток 2 → G
0  / 4 = 0,  остаток 0 → A  (старший разряд)

Результат (справа налево): AGTC
```

**Применение в биоинформатике:**
- Восстановление k-мера по индексу в массиве частот
- Перечисление всех k-меров в лексикографическом порядке
- Парная функция к PatternToNumber

## Алгоритм

### Идея
Последовательно делим число на 4, записывая остатки как нуклеотиды от младшего разряда к старшему.

### Псевдокод
```
NumberToPattern(index, k):
    pattern ← строка из k символов 'A'
    for i ← k-1 downto 0:
        pattern[i] ← NumberToNucleotide(index mod 4)
        index ← index / 4
    return pattern
```

## Реализация на C++

```cpp
char numberToNucleotide(int n) {
    switch (n) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return '?';
    }
}

std::string numberToPattern(int64_t index, size_t k) {
    std::string pattern(k, 'A');  // (1)
    for (size_t i = k; i > 0; --i) {  // (2)
        pattern[i - 1] = numberToNucleotide(index % 4);  // (3)
        index /= 4;  // (4)
    }
    return pattern;
}
```

### Особенности реализации

**1. Конструктор заполнения `std::string`**
```cpp
std::string pattern(k, 'A');
```
- Конструктор `string(count, char)` создаёт строку из `count` копий символа
- Выделяет память один раз — нет реаллокаций при `push_back`
- Далее символы перезаписываются через `pattern[i] = ...`
- Альтернатива: `std::string pattern; pattern.resize(k);` — но символы не инициализированы

**2. Обратный цикл с `size_t` — ловушка underflow**
```cpp
for (size_t i = k; i > 0; --i)
    pattern[i - 1] = ...;
```
- `size_t` — **беззнаковый** тип (обычно `unsigned long`)
- Написать `for (size_t i = k-1; i >= 0; --i)` **НЕЛЬЗЯ**: при `i == 0` декремент даёт `SIZE_MAX` (≈ 1.8×10¹⁹), цикл не завершится
- Безопасный паттерн: итерируем `i` от `k` до `1`, индексируем через `i - 1`
- Другой вариант: `for (int i = k-1; i >= 0; --i)` — но смешивать `int` и `size_t` опасно

**3. Операторы `%` и `/` для целых чисел**
```cpp
index % 4   // остаток от деления
index /= 4  // целочисленное деление (отбрасывает дробную часть)
```
- `%` — оператор остатка (modulo): `45 % 4 = 1`, `11 % 4 = 3`
- `/` для целых чисел усекает к нулю: `45 / 4 = 11` (не 11.25)
- `/=` — составной оператор: `index /= 4` эквивалентно `index = index / 4`
- Для степеней двойки: `% 4` эквивалентно `& 3`, а `/ 4` эквивалентно `>> 2`

**4. Парсинг аргументов CLI**
```cpp
index = std::stoll(argv[1]);  // string to long long
k = std::stoul(argv[2]);      // string to unsigned long
```
- `std::stoll` — преобразует строку в `long long` (≥ 64 бит), подходит для `int64_t`
- `std::stoul` — преобразует строку в `unsigned long`, подходит для `size_t`
- При ошибке бросают `std::invalid_argument` или `std::out_of_range`
- `argv[i]` имеет тип `char*` — неявно конвертируется в `std::string`

## Альтернативные реализации

**Рекурсивная версия:**
```cpp
std::string numberToPattern(int64_t index, size_t k) {
    if (k == 0) return "";
    char last = numberToNucleotide(index % 4);
    return numberToPattern(index / 4, k - 1) + last;
}
```
- Соответствует математическому определению
- Менее эффективна: конкатенация строк на каждом шаге — O(k²) суммарно

## Запуск

```bash
# Без аргументов — тестовые данные
./number_to_pattern

# Аргументы: <index> <k>
./number_to_pattern 45 4
# Вывод: AGTC
```

## Сложность

### Время: O(k)
- Один проход по длине k
- На каждой итерации: деление, остаток, маппинг — O(1)

### Память: O(k)
- Строка результата длины k
