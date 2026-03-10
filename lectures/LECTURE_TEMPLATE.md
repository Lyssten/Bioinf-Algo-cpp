# Lecture XX: Topic Name

## Задачи

| № | Задача | Описание |
|---|--------|----------|
| 1 | task_name/ | Краткое описание |

## Добавление задачи

1. Создайте папку `lectures/lecture_XX/<task_name>/`
2. Добавьте `main.cpp` с решением
3. Сразу реализуйте единый CLI-контракт:
   - входной файл Rosalind через `--input <path>` или просто `<path>` как единственный аргумент
   - positional-аргументы в порядке Rosalind как запасной режим
   - вывод ответа в `stdout` и в файл через `--output <path>` или дефолтный `<input>.out`
4. Создайте `CMakeLists.txt`:
   ```cmake
   add_executable(<task_name> main.cpp)
   ```
   После общей сборки бинарник должен появляться прямо в `lectures/lecture_XX/<task_name>/`
5. Добавьте `add_subdirectory(lectures/lecture_XX/<task_name>)` в корневой `CMakeLists.txt`
6. Добавьте пример входного файла `rosalind_*.txt`, если задача уже есть в Rosalind
7. Напишите `README.md` с описанием алгоритма, расшифровкой параметров (`k`, `t`, `N`, `d` и т.д.) и примерами запуска через файл
