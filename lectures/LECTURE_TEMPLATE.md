# Lecture XX: Topic Name

## Задачи

| № | Задача | Описание |
|---|--------|----------|
| 1 | task_name/ | Краткое описание |

## Добавление задачи

1. Создайте папку `lectures/lecture_XX/<task_name>/`
2. Добавьте `main.cpp` с решением
3. Создайте `CMakeLists.txt`:
   ```cmake
   add_executable(<task_name> main.cpp)
   ```
4. Добавьте `add_subdirectory(lectures/lecture_XX/<task_name>)` в корневой `CMakeLists.txt`
5. Напишите `README.md` с описанием алгоритма
