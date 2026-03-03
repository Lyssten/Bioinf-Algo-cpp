#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <climits>

int nucleotideIndex(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

std::vector<std::vector<double>> buildProfileWithPseudocounts(const std::vector<std::string>& motifs) {
    size_t k = motifs[0].size();
    size_t t = motifs.size();
    std::vector<std::vector<double>> profile(4, std::vector<double>(k, 1.0));

    for (const auto& motif : motifs) {
        for (size_t j = 0; j < k; ++j)
            ++profile[nucleotideIndex(motif[j])][j];
    }

    double total = static_cast<double>(t) + 4.0;
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < k; ++j)
            profile[i][j] /= total;

    return profile;
}

double probability(const std::string& kmer, const std::vector<std::vector<double>>& profile) {
    double prob = 1.0;
    for (size_t i = 0; i < kmer.size(); ++i)
        prob *= profile[nucleotideIndex(kmer[i])][i];
    return prob;
}

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

int score(const std::vector<std::string>& motifs) {
    size_t k = motifs[0].size();
    int totalScore = 0;
    for (size_t j = 0; j < k; ++j) {
        int count[4] = {0, 0, 0, 0};
        for (const auto& motif : motifs)
            ++count[nucleotideIndex(motif[j])];
        int maxCount = 0;
        for (int i = 0; i < 4; ++i)
            maxCount = std::max(maxCount, count[i]);
        totalScore += static_cast<int>(motifs.size()) - maxCount;
    }
    return totalScore;
}

// Одна итерация RandomizedMotifSearch
std::vector<std::string> randomizedMotifSearchOnce(const std::vector<std::string>& dna,
                                                    size_t k, size_t t, std::mt19937& rng) {
    // Случайно выбираем начальные k-меры
    std::vector<std::string> motifs;
    for (size_t i = 0; i < t; ++i) {
        std::uniform_int_distribution<size_t> dist(0, dna[i].size() - k);
        size_t pos = dist(rng);
        motifs.push_back(dna[i].substr(pos, k));
    }

    std::vector<std::string> bestMotifs = motifs;

    while (true) {
        auto profile = buildProfileWithPseudocounts(motifs);
        for (size_t i = 0; i < t; ++i)
            motifs[i] = profileMostProbable(dna[i], k, profile);

        if (score(motifs) < score(bestMotifs))
            bestMotifs = motifs;
        else
            return bestMotifs;
    }
}

// Запускаем 1000 раз и берём лучший результат
std::vector<std::string> randomizedMotifSearch(const std::vector<std::string>& dna,
                                                size_t k, size_t t) {
    std::mt19937 rng(42);
    std::vector<std::string> bestMotifs;
    int bestScore = INT_MAX;

    for (int iter = 0; iter < 1000; ++iter) {
        auto motifs = randomizedMotifSearchOnce(dna, k, t, rng);
        int s = score(motifs);
        if (s < bestScore) {
            bestScore = s;
            bestMotifs = motifs;
        }
    }
    return bestMotifs;
}

int main(int argc, char* argv[]) {
    size_t k = 8, t = 5;
    std::vector<std::string> dna = {
        "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
        "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
        "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
        "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
        "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
    };

    if (argc >= 3) {
        k = std::stoul(argv[1]);
        t = std::stoul(argv[2]);
        dna.clear();
        for (int i = 3; i < argc; ++i)
            dna.push_back(argv[i]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto motifs = randomizedMotifSearch(dna, k, t);
    for (const auto& m : motifs)
        std::cout << m << std::endl;

    return 0;
}
