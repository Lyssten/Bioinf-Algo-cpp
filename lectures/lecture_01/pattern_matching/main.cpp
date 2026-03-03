#include <iostream>
#include <string>
#include <vector>

std::vector<int> findPatternPositions(const std::string& genome, const std::string& pattern) {
    std::vector<int> positions;
    
    if (pattern.empty() || pattern.size() > genome.size())
        return positions;

    for (size_t i = 0; i <= genome.size() - pattern.size(); ++i) {
        if (genome.substr(i, pattern.size()) == pattern)
            positions.push_back(i);
    }
    return positions;
}

int main(int argc, char* argv[]) {
    std::string pattern = "ATAT";
    std::string genome = "GATATATGCATATACTT";

    if (argc == 3) {
        pattern = argv[1];
        genome = argv[2];
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto positions = findPatternPositions(genome, pattern);
    for (size_t i = 0; i < positions.size(); ++i) {
        if (i > 0) std::cout << " ";
        std::cout << positions[i];
    }
    std::cout << std::endl;

    return 0;
}
