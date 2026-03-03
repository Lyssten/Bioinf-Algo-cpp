#include <iostream>
#include <string>
#include <algorithm>

std::string reverseComplement(const std::string& text) {
    std::string result;
    result.reserve(text.size());

    for (auto it = text.rbegin(); it != text.rend(); ++it) {
        switch (*it) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default: result += *it;
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    std::string text = "AAAACCCGGT";

    if (argc == 2) {
        text = argv[1];
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    std::cout << reverseComplement(text) << std::endl;
    return 0;
}
