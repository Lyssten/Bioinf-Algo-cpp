#include <iostream>
#include <string>
#include <vector>
#include <set>

int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

// Генерируем все k-меры с расстоянием Хэмминга ≤ d от pattern
void generateNeighbors(const std::string& pattern, size_t d, std::set<std::string>& result) {
    auto generate = [&](auto&& self, size_t pos, std::string& current, size_t mismatches) -> void {
        if (mismatches > d) return;
        if (pos == current.size()) {
            result.insert(current);
            return;
        }
        char orig = pattern[pos];
        for (char c : {'A', 'C', 'G', 'T'}) {
            current[pos] = c;
            size_t newMismatches = mismatches + (c != orig ? 1 : 0);
            if (newMismatches <= d)
                self(self, pos + 1, current, newMismatches);
        }
        current[pos] = orig;
    };

    std::string mutable_pattern = pattern;
    generate(generate, 0, mutable_pattern, 0);
}

// Проверяем, встречается ли pattern в text с ≤ d несовпадениями
bool appearsWithMismatches(const std::string& text, const std::string& pattern, size_t d) {
    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {
        if (hammingDistance(text.substr(i, pattern.size()), pattern) <= static_cast<int>(d))
            return true;
    }
    return false;
}

// Находим все (k,d)-мотивы — k-меры, встречающиеся в каждой строке с ≤ d mismatches
std::set<std::string> motifEnumeration(const std::vector<std::string>& dna, size_t k, size_t d) {
    std::set<std::string> motifs;

    // Генерируем кандидатов из первой строки
    for (size_t i = 0; i <= dna[0].size() - k; ++i) {
        std::string kmer = dna[0].substr(i, k);
        std::set<std::string> neighbors;
        generateNeighbors(kmer, d, neighbors);

        for (const auto& candidate : neighbors) {
            bool inAll = true;
            for (size_t j = 1; j < dna.size(); ++j) {
                if (!appearsWithMismatches(dna[j], candidate, d)) {
                    inAll = false;
                    break;
                }
            }
            if (inAll)
                motifs.insert(candidate);
        }
    }
    return motifs;
}

int main(int argc, char* argv[]) {
    size_t k = 3, d = 1;
    std::vector<std::string> dna = {"ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"};

    if (argc >= 3) {
        k = std::stoul(argv[1]);
        d = std::stoul(argv[2]);
        dna.clear();
        for (int i = 3; i < argc; ++i)
            dna.push_back(argv[i]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto motifs = motifEnumeration(dna, k, d);
    bool first = true;
    for (const auto& m : motifs) {
        if (!first) std::cout << " ";
        std::cout << m;
        first = false;
    }
    std::cout << std::endl;

    return 0;
}
