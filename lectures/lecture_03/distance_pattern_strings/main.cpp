#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <climits>

int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

// Сумма минимальных расстояний Хэмминга между pattern и каждой строкой из dna
int distanceBetweenPatternAndStrings(const std::string& pattern, const std::vector<std::string>& dna) {
    size_t k = pattern.size();
    int totalDist = 0;

    for (const auto& text : dna) {
        int minDist = INT_MAX;
        for (size_t i = 0; i <= text.size() - k; ++i) {
            int dist = hammingDistance(pattern, text.substr(i, k));
            minDist = std::min(minDist, dist);
        }
        totalDist += minDist;
    }
    return totalDist;
}

int main(int argc, char* argv[]) {
    std::string pattern = "AAA";
    std::vector<std::string> dna = {"TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"};

    if (argc >= 3) {
        pattern = argv[1];
        dna.clear();
        // Вторым аргументом передаётся строка с ДНК через пробелы
        std::istringstream iss(argv[2]);
        std::string s;
        while (iss >> s)
            dna.push_back(s);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << distanceBetweenPatternAndStrings(pattern, dna) << std::endl;

    return 0;
}
