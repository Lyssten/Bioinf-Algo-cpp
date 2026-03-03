#include <iostream>
#include <string>
#include <vector>

std::vector<int> minimumSkew(const std::string& genome) {
    std::vector<int> positions;
    int skew = 0, minSkew = 0;

    positions.push_back(0); // Skew(0) = 0

    for (size_t i = 0; i < genome.size(); ++i) {
        if (genome[i] == 'G') ++skew;
        else if (genome[i] == 'C') --skew;

        if (skew < minSkew) {
            minSkew = skew;
            positions.clear();
            positions.push_back(i + 1);
        } else if (skew == minSkew) {
            positions.push_back(i + 1);
        }
    }
    return positions;
}

int main(int argc, char* argv[]) {
    std::string genome = "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG";

    if (argc == 2) {
        genome = argv[1];
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto positions = minimumSkew(genome);
    for (size_t i = 0; i < positions.size(); ++i) {
        if (i > 0) std::cout << " ";
        std::cout << positions[i];
    }
    std::cout << std::endl;

    return 0;
}
