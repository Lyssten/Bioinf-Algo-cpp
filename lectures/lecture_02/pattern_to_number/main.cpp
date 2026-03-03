#include <iostream>
#include <string>
#include <cstdint>

// Преобразуем нуклеотид в число: A=0, C=1, G=2, T=3
int nucleotideToNumber(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Преобразуем ДНК-строку в число (позиционная система счисления с основанием 4)
int64_t patternToNumber(const std::string& pattern) {
    int64_t number = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        number = number * 4 + nucleotideToNumber(pattern[i]);
    }
    return number;
}

int main(int argc, char* argv[]) {
    std::string pattern = "AGT";

    if (argc == 2) {
        pattern = argv[1];
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << patternToNumber(pattern) << std::endl;

    return 0;
}
