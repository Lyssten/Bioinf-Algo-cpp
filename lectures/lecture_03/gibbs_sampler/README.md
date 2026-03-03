# Gibbs Sampler (BA2G)

## Задача

Найти набор мотивов алгоритмом Гиббса.

**Вход:** целые числа k, t, N и набор строк Dna
**Выход:** лучший набор мотивов за 20 запусков GibbsSampler по N итераций каждый

**Пример:**
```
k: 8, t: 5, N: 100
Dna: CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
     GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
     ...
Вывод:
TCTCGGGG
CCAAGGTG
TACAGGCG
TTCAGGTG
TCCACGTG
```

## Теория

### Отличие от RandomizedMotifSearch

RandomizedMotifSearch заменяет **все** мотивы одновременно на каждой итерации. Gibbs Sampler заменяет **один случайный** мотив за итерацию — более плавное изменение.

### Profile-randomly generated k-mer

В отличие от Profile-**most** probable, здесь k-мер выбирается **случайно** пропорционально вероятности:
```
P(выбрать k-мер_i) = Pr(k-мер_i | Profile) / Σ Pr(k-мер_j | Profile)
```

Это позволяет «исследовать» пространство поиска и избегать ранней сходимости к плохому локальному минимуму.

### Сравнение подходов

| Свойство | Greedy | Randomized | Gibbs |
|----------|--------|------------|-------|
| Начало | детерминированное | случайное | случайное |
| Обновление | последовательное | все сразу | один случайный |
| Выбор k-мера | most probable | most probable | randomly generated |
| Сходимость | быстрая, но к ближайшему минимуму | быстрая, много рестартов | медленная, но более точная |

## Алгоритм

### Псевдокод
```
GibbsSampler(Dna, k, t, N):
    случайно выбрать Motifs
    BestMotifs ← Motifs
    for j ← 1 to N:
        i ← Random(0, t-1)  // выбираем строку для замены
        Profile ← ProfileWithPseudocounts(Motifs без Motifs[i])
        Motifs[i] ← ProfileRandomlyGenerated(Dna[i], k, Profile)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs ← Motifs
    return BestMotifs

// Основной цикл: 20 рестартов
best ← GibbsSampler с худшим Score
repeat 20 times:
    result ← GibbsSampler(Dna, k, t, N)
    if Score(result) < Score(best):
        best ← result
```

## Реализация на C++

```cpp
std::string profileRandomlyGenerated(const std::string& text, size_t k,
                                      const std::vector<std::vector<double>>& profile,
                                      std::mt19937& rng) {
    std::vector<double> probs;
    double total = 0.0;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        double p = probability(text.substr(i, k), profile);
        probs.push_back(p);
        total += p;  // (1)
    }

    std::uniform_real_distribution<double> dist(0.0, total);
    double r = dist(rng);
    double cumulative = 0.0;

    for (size_t i = 0; i < probs.size(); ++i) {
        cumulative += probs[i];
        if (cumulative >= r)  // (2)
            return text.substr(i, k);
    }
    return text.substr(text.size() - k, k);
}
```

### Особенности реализации

**1. `std::uniform_real_distribution<double>` — непрерывное распределение**
```cpp
std::uniform_real_distribution<double> dist(0.0, total);
double r = dist(rng);
```
- В отличие от `uniform_int_distribution`, генерирует вещественные числа
- Диапазон [0.0, total) — включая 0, не включая total
- `double` — тип генерируемых чисел (можно `float` для экономии)
- Используется для взвешенного случайного выбора

**2. Взвешенная случайная выборка (кумулятивный метод)**
```cpp
std::vector<double> probs;   // ненормализованные веса
double total = 0.0;
for (...) {
    double p = probability(kmer, profile);
    probs.push_back(p);
    total += p;
}

double r = dist(rng);        // случайная точка на [0, total)
double cumulative = 0.0;
for (size_t i = 0; i < probs.size(); ++i) {
    cumulative += probs[i];
    if (cumulative >= r)      // нашли интервал
        return text.substr(i, k);
}
```
- Алгоритм: каждый k-мер «занимает» отрезок на числовой прямой [0, total)
- Длина отрезка = вероятность k-мера (чем вероятнее — тем длиннее отрезок)
- Генерируем случайную точку `r` → находим, в чей отрезок она попала
- **Не нужна нормализация:** `dist(0, total)` уже учитывает ненормализованную сумму
- Время выборки: O(n) (один проход по массиву вероятностей)

**3. Построение профиля с исключением — паттерн `continue`**
```cpp
for (size_t i = 0; i < motifs.size(); ++i) {
    if (i == exclude) continue;  // пропускаем i-ю строку
    for (size_t j = 0; j < k; ++j)
        ++profile[nucleotideIndex(motifs[i][j])][j];
}
```
- `continue` — переходит к следующей итерации цикла, пропуская тело
- Здесь: строим профиль из всех мотивов, кроме исключённого
- **Альтернатива:** два цикла `[0, exclude)` и `(exclude, size)` — менее читаемо
- `motifs.size() - 1` в знаменателе нормализации (t-1 строк вместо t)

**4. Защита от граничного случая**
```cpp
return text.substr(text.size() - k, k);  // fallback
```
- Из-за особенностей `double` кумулятивная сумма может не достичь `r` точно
- Последний k-мер возвращается как запасной вариант (на практике редко срабатывает)
- Defensive programming: код корректен даже при ошибках округления

## Запуск

```bash
./gibbs_sampler 8 5 100 "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA" "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG" ...
```

## Сложность

### Время: O(20 × N × (t × k + n × k))
- 20 рестартов × N итераций
- На каждой итерации: построение профиля O(t×k) + выбор k-мера O(n×k)

### Память: O(t × k + n)
- Мотивы + массив вероятностей для одной строки
