#include <iostream>
#include <string>
#include <unordered_map>
#include <set>

std::set<std::string> frequentWords(const std::string& text, size_t k) {
    if (k == 0 || k > text.size())
        return {};

    std::unordered_map<std::string, int> freq;
    int maxCount = 0;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        std::string pattern = text.substr(i, k);
        int count = ++freq[pattern];
        maxCount = std::max(maxCount, count);
    }

    std::set<std::string> result;
    for (const auto& [pattern, count] : freq) {
        if (count == maxCount)
            result.insert(pattern);
    }
    return result;
}

int main(int argc, char* argv[]) {
    std::string text = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
    size_t k = 4;

    if (argc == 3) {
        text = argv[1];
        k = std::stoul(argv[2]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto patterns = frequentWords(text, k);
    for (auto it = patterns.begin(); it != patterns.end(); ++it) {
        if (it != patterns.begin()) std::cout << " ";
        std::cout << *it;
    }
    std::cout << std::endl;

    return 0;
}
