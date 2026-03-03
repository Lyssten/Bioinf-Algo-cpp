#include <iostream>
#include <string>
#include <vector>

int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

// Находим все позиции, где pattern встречается в text с ≤ d несовпадениями
std::vector<int> approximatePatternMatching(const std::string& text, const std::string& pattern, size_t d) {
    std::vector<int> positions;
    if (pattern.empty() || pattern.size() > text.size())
        return positions;

    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {
        std::string window = text.substr(i, pattern.size());
        if (hammingDistance(window, pattern) <= static_cast<int>(d))
            positions.push_back(static_cast<int>(i));
    }
    return positions;
}

int main(int argc, char* argv[]) {
    std::string pattern = "ATTCTGGA";
    std::string text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC";
    size_t d = 3;

    if (argc == 4) {
        pattern = argv[1];
        text = argv[2];
        d = std::stoul(argv[3]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto positions = approximatePatternMatching(text, pattern, d);
    for (size_t i = 0; i < positions.size(); ++i) {
        if (i > 0) std::cout << " ";
        std::cout << positions[i];
    }
    std::cout << std::endl;

    return 0;
}
