#include <iostream>
#include <string>
#include <vector>
#include <sstream>

// Индекс нуклеотида в профильной матрице: A=0, C=1, G=2, T=3
int nucleotideIndex(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// Вычисляем вероятность k-мера по профильной матрице
double probability(const std::string& kmer, const std::vector<std::vector<double>>& profile) {
    double prob = 1.0;
    for (size_t i = 0; i < kmer.size(); ++i) {
        prob *= profile[nucleotideIndex(kmer[i])][i];
    }
    return prob;
}

// Находим k-мер в text с максимальной вероятностью по profile
std::string profileMostProbable(const std::string& text, size_t k,
                                const std::vector<std::vector<double>>& profile) {
    double maxProb = -1.0;
    std::string bestKmer = text.substr(0, k);

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string kmer = text.substr(i, k);
        double prob = probability(kmer, profile);
        if (prob > maxProb) {
            maxProb = prob;
            bestKmer = kmer;
        }
    }
    return bestKmer;
}

int main(int argc, char* argv[]) {
    std::string text = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT";
    size_t k = 5;
    std::vector<std::vector<double>> profile = {
        {0.2, 0.2, 0.3, 0.2, 0.3},  // A
        {0.4, 0.3, 0.1, 0.5, 0.1},  // C
        {0.3, 0.3, 0.5, 0.2, 0.4},  // G
        {0.1, 0.2, 0.1, 0.1, 0.2}   // T
    };

    if (argc >= 3) {
        text = argv[1];
        k = std::stoul(argv[2]);
        profile.assign(4, std::vector<double>());
        // Профиль передаётся как 4 строки: каждая строка — вероятности через пробел
        for (int row = 0; row < 4 && (3 + row) < argc; ++row) {
            std::istringstream iss(argv[3 + row]);
            double val;
            while (iss >> val)
                profile[row].push_back(val);
        }
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << profileMostProbable(text, k, profile) << std::endl;

    return 0;
}
