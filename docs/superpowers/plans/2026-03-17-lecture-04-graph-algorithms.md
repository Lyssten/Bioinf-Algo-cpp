# Lecture 04: Graph Algorithms — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add lecture_04 with 6 Rosalind problems (BA3A–BA3F) about graph-based genome assembly algorithms.

**Architecture:** Each task is a self-contained `main.cpp` following the lecture_03 CLI pattern (CliOptions, parseInputFile, parseCli, sampleInput, validateInput, solve, writeOutputFile). No shared headers. READMEs follow the project's 8-section style (Задача, Теория, Алгоритм, Реализация на C++, Особенности, Альтернативные реализации, Запуск, Сложность).

**Tech Stack:** C++17, CMake 3.10+

**Delegation Strategy:**
- **Gemini CLI** — generates README.md files (template-driven, ~100-200 lines each, verifiable against existing READMEs)
- **Claude agents** — implement main.cpp files (require understanding of the full CLI boilerplate pattern)

---

## Chunk 1: Task Implementations (main.cpp + CMakeLists.txt)

### Task 1: string_composition (BA3A)

**Files:**
- Create: `lectures/lecture_04/string_composition/main.cpp`
- Create: `lectures/lecture_04/string_composition/CMakeLists.txt`

**Problem:** Given integer k and string Text, return all k-mers of Text (one per line, sorted).

**Input format:** Line 1: k. Line 2: Text.
**Sample input:** k=5, Text="CAATCCAAC" → AATCC ATCCA CAATC CCAAC TCCAA

**Algorithm:** Slide window of size k, collect substrings, sort lexicographically.

- [ ] **Step 1:** Create `CMakeLists.txt` — `add_executable(string_composition main.cpp)`
- [ ] **Step 2:** Write `main.cpp` with full CLI boilerplate + algorithm
  - InputData: `{size_t k; std::string text;}`
  - solve: extract all k-mers, sort, join with newlines
  - sampleInput: k=5, text="CAATCCAAC"
- [ ] **Step 3:** Build and verify with sample data

---

### Task 2: genome_path_string (BA3B)

**Files:**
- Create: `lectures/lecture_04/genome_path_string/main.cpp`
- Create: `lectures/lecture_04/genome_path_string/CMakeLists.txt`

**Problem:** Given a sequence of k-mers forming a genome path, reconstruct the string.

**Input format:** One k-mer per line.
**Sample input:** ACCGA CCGAA CGAAG GAAGC AAGCT → ACCGAAGCT

**Algorithm:** Take first k-mer fully, then append last char of each subsequent k-mer.

- [ ] **Step 1:** Create `CMakeLists.txt` — `add_executable(genome_path_string main.cpp)`
- [ ] **Step 2:** Write `main.cpp` with full CLI boilerplate + algorithm
  - InputData: `{std::vector<std::string> kmers;}`
  - solve: result = kmers[0]; for i=1..n: result += kmers[i].back()
  - sampleInput: {"ACCGA","CCGAA","CGAAG","GAAGC","AAGCT"}
- [ ] **Step 3:** Build and verify with sample data

---

### Task 3: overlap_graph (BA3C)

**Files:**
- Create: `lectures/lecture_04/overlap_graph/main.cpp`
- Create: `lectures/lecture_04/overlap_graph/CMakeLists.txt`

**Problem:** Given collection of k-mers, build overlap graph (adjacency list). Edge from A→B if suffix(A) == prefix(B).

**Input format:** One k-mer per line.
**Sample:** ATGCG GCATG CATGC AGGCA GGCAT → AGGCA->GGCAT, CATGC->ATGCG, GCATG->CATGC, GGCAT->GCATG

**Algorithm:** For each pair (a,b): if a.substr(1) == b.substr(0,k-1), add edge a→b. Sort output.

- [ ] **Step 1:** Create `CMakeLists.txt`
- [ ] **Step 2:** Write `main.cpp` — use map<string, vector<string>> for adjacency list, sort output
- [ ] **Step 3:** Build and verify

---

### Task 4: de_bruijn_string (BA3D)

**Files:**
- Create: `lectures/lecture_04/de_bruijn_string/main.cpp`
- Create: `lectures/lecture_04/de_bruijn_string/CMakeLists.txt`

**Problem:** Given k and Text, construct DeBruijn_k(Text) as adjacency list.

**Input format:** Line 1: k. Line 2: Text.
**Sample:** k=4, Text="AAGATTCTCTAC" → AAG->AGA, AGA->GAT, ATT->TTC, CTA->TAC, CTC->TCT, GAT->ATT, TCT->CTA,CTC, TTC->TCT

**Algorithm:** Extract (k-1)-mer nodes from each k-mer edge. For each k-mer: prefix→suffix. Group by prefix using map<string, vector<string>>.

- [ ] **Step 1:** Create `CMakeLists.txt`
- [ ] **Step 2:** Write `main.cpp` — edges from k-mers of Text, group by (k-1)-mer prefix
- [ ] **Step 3:** Build and verify

---

### Task 5: de_bruijn_kmers (BA3E)

**Files:**
- Create: `lectures/lecture_04/de_bruijn_kmers/main.cpp`
- Create: `lectures/lecture_04/de_bruijn_kmers/CMakeLists.txt`

**Problem:** Given collection of k-mers, construct De Bruijn graph.

**Input format:** One k-mer per line.
**Sample:** GAGG CAGG GGGG GGGA CAGG AGGG GGAG → AGG->GGG, CAG->AGG,AGG, GAG->AGG, GGA->GAG, GGG->GGA,GGG

**Algorithm:** For each k-mer: add edge prefix→suffix. Group, preserve duplicates.

- [ ] **Step 1:** Create `CMakeLists.txt`
- [ ] **Step 2:** Write `main.cpp` — similar to Task 4, but input is k-mers directly
- [ ] **Step 3:** Build and verify

---

### Task 6: eulerian_cycle (BA3F)

**Files:**
- Create: `lectures/lecture_04/eulerian_cycle/main.cpp`
- Create: `lectures/lecture_04/eulerian_cycle/CMakeLists.txt`

**Problem:** Find Eulerian cycle in directed graph (adjacency list input).

**Input format:** Lines of "node -> neighbor1,neighbor2,..."
**Sample:** 0->3, 1->0, ... → 6->8->7->9->6->5->4->2->1->0->3->2->6

**Algorithm:** Hierholzer's algorithm:
1. Start from any node, walk until stuck (forms initial cycle)
2. Find node in cycle with unused edges
3. Start new walk from that node, splice into cycle
4. Repeat until all edges used

- [ ] **Step 1:** Create `CMakeLists.txt`
- [ ] **Step 2:** Write `main.cpp` — Hierholzer's algorithm with adjacency list I/O parsing
  - InputData: `{std::map<int, std::vector<int>> graph;}`
  - Graph nodes are integers, parse "node -> n1,n2,..." format
  - Output: cycle as "n1->n2->...->n1"
- [ ] **Step 3:** Build and verify

---

## Chunk 2: README.md files (delegated to Gemini)

Each README follows the project style (see `lectures/lecture_03/motif_enumeration/README.md` as reference).

**Delegation:** Use `gemini -y -p "..."` for each README, providing:
- The Rosalind problem description
- The implemented main.cpp as context
- An existing README as style template
- Instructions to write in Russian, following the 8-section structure

- [ ] **Step 7:** Generate README.md for each of the 6 tasks via Gemini
- [ ] **Step 8:** Review generated READMEs for correctness and style consistency

---

## Chunk 3: Project Configuration Updates

### Task 9: Update root CMakeLists.txt

**Files:**
- Modify: `CMakeLists.txt`

- [ ] **Step 1:** Add 6 `add_bioinf_task(...)` lines under new "# Lecture 04: Graph algorithms" comment

---

### Task 10: Update root README.md

**Files:**
- Modify: `README.md`

- [ ] **Step 1:** Add lecture_04 to directory tree
- [ ] **Step 2:** Add Lecture 04 table with all 6 tasks
- [ ] **Step 3:** Add Lecture 04 run examples

---

### Task 11: Update CLAUDE.md

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1:** Add lecture_04 to structure tree
- [ ] **Step 2:** Add Lecture 04 task table
- [ ] **Step 3:** Add lecture_04 utilities list if any new ones

---

## Chunk 4: Build & Verify

- [ ] **Step 1:** Full build: `cmake -S . -B .cmake-build && cmake --build .cmake-build -j$(nproc)`
- [ ] **Step 2:** Run each task with sample data, verify output matches Rosalind samples
- [ ] **Step 3:** Commit all changes
