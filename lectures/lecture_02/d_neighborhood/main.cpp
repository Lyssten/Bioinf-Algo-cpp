#include <iostream>
#include <string>
#include <set>

// Генерируем все строки с расстоянием Хэмминга ≤ d от pattern
std::set<std::string> neighbors(const std::string& pattern, size_t d) {
    std::set<std::string> result;

    if (d == 0) {
        result.insert(pattern);
        return result;
    }

    // Рекурсивная генерация всех соседей
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
    return result;
}

int main(int argc, char* argv[]) {
    std::string pattern = "ACG";
    size_t d = 1;

    if (argc == 3) {
        pattern = argv[1];
        d = std::stoul(argv[2]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto result = neighbors(pattern, d);
    for (const auto& s : result) {
        std::cout << s << std::endl;
    }

    return 0;
}
