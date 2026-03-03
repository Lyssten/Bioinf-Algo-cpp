#include <iostream>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>

int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

// Генерируем все k-меры с ≤ d mismatches
void generateNeighbors(const std::string& pattern, size_t d, std::set<std::string>& result) {
    if (d == 0) {
        result.insert(pattern);
        return;
    }
    
    // Добавляем сам паттерн
    result.insert(pattern);
    
    if (d >= pattern.size()) {
        // Если d >= k, то все возможные k-меры
        std::string current(pattern.size(), 'A');
        std::vector<std::string> all;
        all.reserve(4 * pattern.size());
        
        auto gen = [&](auto&& self, size_t pos, std::string& s) -> void {
            if (pos == s.size()) {
                result.insert(s);
                return;
            }
            for (char c : {'A', 'C', 'G', 'T'}) {
                s[pos] = c;
                self(self, pos + 1, s);
            }
        };
        gen(gen, 0, current);
        return;
    }
    
    // Рекурсивная генерация
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

std::set<std::string> frequentWordsWithMismatches(const std::string& text, size_t k, size_t d) {
    std::unordered_map<std::string, int> freq;
    int maxCount = 0;

    // Для каждого k-мера в тексте генерируем все k-меры с ≤ d mismatches
    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        std::set<std::string> neighbors;
        generateNeighbors(pattern, d, neighbors);
        
        for (const auto& neighbor : neighbors) {
            int count = ++freq[neighbor];
            maxCount = std::max(maxCount, count);
        }
    }

    std::set<std::string> result;
    for (const auto& [p, count] : freq) {
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

    auto patterns = frequentWordsWithMismatches(text, k, d);
    bool first = true;
    for (const auto& p : patterns) {
        if (!first) std::cout << " ";
        std::cout << p;
        first = false;
    }
    std::cout << std::endl;

    return 0;
}
