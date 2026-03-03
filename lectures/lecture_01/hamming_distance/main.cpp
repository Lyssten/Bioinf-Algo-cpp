#include <iostream>
#include <string>

int hammingDistance(const std::string& s1, const std::string& s2) {
    if (s1.size() != s2.size())
        return -1;
    
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

int main(int argc, char* argv[]) {
    std::string s1 = "GGGCCGTTGGT", s2 = "GGACCGTTGAC";

    if (argc == 3) {
        s1 = argv[1];
        s2 = argv[2];
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    int dist = hammingDistance(s1, s2);
    if (dist >= 0)
        std::cout << dist << std::endl;
    else
        std::cerr << "Error: strings must have equal length" << std::endl;

    return 0;
}
