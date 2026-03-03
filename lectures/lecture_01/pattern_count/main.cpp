#include <iostream>
#include <string>

int patternCount(const std::string& text, const std::string& pattern) {
    if (pattern.empty() || pattern.size() > text.size())
        return 0;

    int count = 0;
    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {
        if (text.substr(i, pattern.size()) == pattern)
            ++count;
    }
    return count;
}

int main(int argc, char* argv[]) {
    std::string text = "GCGCG", pattern = "GCG";

    if (argc == 3) {
        text = argv[1];
        pattern = argv[2];
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << patternCount(text, pattern) << std::endl;
    return 0;
}
