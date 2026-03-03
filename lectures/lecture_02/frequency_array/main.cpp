#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdint>

int nucleotideToNumber(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

int64_t patternToNumber(const std::string& pattern) {
    int64_t number = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        number = number * 4 + nucleotideToNumber(pattern[i]);
    }
    return number;
}

// Генерируем массив частот всех k-меров в тексте
// Индекс = PatternToNumber(k-мер), значение = количество вхождений
std::vector<int> computeFrequencyArray(const std::string& text, size_t k) {
    size_t arraySize = static_cast<size_t>(std::pow(4, k));
    std::vector<int> freqArray(arraySize, 0);

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        int64_t index = patternToNumber(pattern);
        ++freqArray[index];
    }
    return freqArray;
}

int main(int argc, char* argv[]) {
    std::string text = "ACGCGGCTCTGAAA";
    size_t k = 2;

    if (argc == 3) {
        text = argv[1];
        k = std::stoul(argv[2]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto freqArray = computeFrequencyArray(text, k);
    for (size_t i = 0; i < freqArray.size(); ++i) {
        if (i > 0) std::cout << " ";
        std::cout << freqArray[i];
    }
    std::cout << std::endl;

    return 0;
}
