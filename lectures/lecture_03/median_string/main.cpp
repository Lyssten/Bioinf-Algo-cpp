#include <iostream>
#include <string>
#include <vector>
#include <climits>
#include <cstdint>

int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

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
    std::string pattern(k, 'A');
    for (size_t i = k; i > 0; --i) {
        pattern[i - 1] = numberToNucleotide(index % 4);
        index /= 4;
    }
    return pattern;
}

// Перебираем все 4^k k-меров и находим тот, что минимизирует d(Pattern, Dna)
std::string medianString(const std::vector<std::string>& dna, size_t k) {
    int minDist = INT_MAX;
    std::string median;

    int64_t total = 1;
    for (size_t i = 0; i < k; ++i) total *= 4;

    for (int64_t i = 0; i < total; ++i) {
        std::string pattern = numberToPattern(i, k);
        int dist = distanceBetweenPatternAndStrings(pattern, dna);
        if (dist < minDist) {
            minDist = dist;
            median = pattern;
        }
    }
    return median;
}

int main(int argc, char* argv[]) {
    size_t k = 3;
    std::vector<std::string> dna = {
        "AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG",
        "GCTGAGCACCGG", "AGTACGGGACAG"
    };

    if (argc >= 3) {
        k = std::stoul(argv[1]);
        dna.clear();
        for (int i = 2; i < argc; ++i)
            dna.push_back(argv[i]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << medianString(dna, k) << std::endl;

    return 0;
}
