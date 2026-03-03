#include <iostream>
#include <string>
#include <unordered_map>
#include <set>

// Обратный комплемент нуклеотида
char complement(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return c;
    }
}

// Обратный комплемент строки
std::string reverseComplement(const std::string& pattern) {
    std::string result;
    result.reserve(pattern.size());
    for (auto it = pattern.rbegin(); it != pattern.rend(); ++it) {
        result += complement(*it);
    }
    return result;
}

// Генерируем все строки с расстоянием Хэмминга ≤ d от pattern
void generateNeighbors(const std::string& pattern, size_t d, std::set<std::string>& result) {
    if (d == 0) {
        result.insert(pattern);
        return;
    }

    auto generate = [&](auto&& self, size_t pos, std::string& current, size_t mismatches) -> void {
        if (mismatches > d) return;
        if (pos == current.size()) {
            result.insert(current);
            return;
        }

        char orig = pattern[pos];
        for (char c : {'A', 'C', 'G', 'T'}) {
            current[pos] = c;
            size_t newMismatches = mismatches + (c != orig ? 1 : 0);
            if (newMismatches <= d)
                self(self, pos + 1, current, newMismatches);
        }
        current[pos] = orig;
    };

    std::string mutablePattern = pattern;
    generate(generate, 0, mutablePattern, 0);
}

// Наиболее частые k-меры с учётом mismatches И обратных комплементов
std::set<std::string> frequentWordsWithMismatchesAndComplements(const std::string& text, size_t k, size_t d) {
    // Шаг 1: считаем Count_d(Text, Pattern) для каждого Pattern
    std::unordered_map<std::string, int> freq;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        std::set<std::string> neighbors;
        generateNeighbors(pattern, d, neighbors);

        for (const auto& neighbor : neighbors) {
            ++freq[neighbor];
        }
    }

    // Шаг 2: для каждого кандидата суммируем freq[P] + freq[RC(P)]
    int maxCount = 0;
    std::unordered_map<std::string, int> combined;

    for (const auto& [p, count] : freq) {
        if (combined.count(p)) continue;
        std::string rc = reverseComplement(p);
        int total = count + (rc != p ? freq[rc] : 0);
        combined[p] = total;
        combined[rc] = total;
        maxCount = std::max(maxCount, total);
    }

    std::set<std::string> result;
    for (const auto& [p, count] : combined) {
        if (count == maxCount)
            result.insert(p);
    }
    return result;
}

int main(int argc, char* argv[]) {
    std::string text = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
    size_t k = 4, d = 1;

    if (argc == 4) {
        text = argv[1];
        k = std::stoul(argv[2]);
        d = std::stoul(argv[3]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto patterns = frequentWordsWithMismatchesAndComplements(text, k, d);
    bool first = true;
    for (const auto& p : patterns) {
        if (!first) std::cout << " ";
        std::cout << p;
        first = false;
    }
    std::cout << std::endl;

    return 0;
}
