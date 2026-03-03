#include <iostream>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>

std::set<std::string> clumpFinding(const std::string& genome, size_t k, size_t L, size_t t) {
    std::set<std::string> result;
    
    if (k > genome.size() || L > genome.size())
        return result;

    for (size_t i = 0; i <= genome.size() - L; ++i) {
        std::string window = genome.substr(i, L);
        std::unordered_map<std::string, int> freq;
        
        for (size_t j = 0; j <= L - k; ++j) {
            std::string pattern = window.substr(j, k);
            if (++freq[pattern] >= t)
                result.insert(pattern);
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    std::string genome = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC";
    size_t k = 5, L = 75, t = 4;

    if (argc == 5) {
        genome = argv[1];
        k = std::stoul(argv[2]);
        L = std::stoul(argv[3]);
        t = std::stoul(argv[4]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto patterns = clumpFinding(genome, k, L, t);
    bool first = true;
    for (const auto& p : patterns) {
        if (!first) std::cout << " ";
        std::cout << p;
        first = false;
    }
    std::cout << std::endl;

    return 0;
}
