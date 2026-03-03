#include <iostream>
#include <string>
#include <cstdint>

// Преобразуем число в нуклеотид: 0=A, 1=C, 2=G, 3=T
char numberToNucleotide(int n) {
    switch (n) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return '?';
    }
}

// Преобразуем число обратно в ДНК-строку длины k
std::string numberToPattern(int64_t index, size_t k) {
    std::string pattern(k, 'A');
    for (size_t i = k; i > 0; --i) {
        pattern[i - 1] = numberToNucleotide(index % 4);
        index /= 4;
    }
    return pattern;
}

int main(int argc, char* argv[]) {
    int64_t index = 45;
    size_t k = 4;

    if (argc == 3) {
        index = std::stoll(argv[1]);
        k = std::stoul(argv[2]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << numberToPattern(index, k) << std::endl;

    return 0;
}
