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

// Профиль с псевдосчётчиками, исключая motifs[exclude]
std::vector<std::vector<double>> buildProfileExcluding(const std::vector<std::string>& motifs,
                                                        size_t exclude) {
    size_t k = motifs[0].size();
    size_t t = motifs.size() - 1;  // без исключённой строки
    std::vector<std::vector<double>> profile(4, std::vector<double>(k, 1.0));

    for (size_t i = 0; i < motifs.size(); ++i) {
        if (i == exclude) continue;
        for (size_t j = 0; j < k; ++j)
            ++profile[nucleotideIndex(motifs[i][j])][j];
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

// Случайно выбираем k-мер из text пропорционально вероятностям по profile
std::string profileRandomlyGenerated(const std::string& text, size_t k,
                                      const std::vector<std::vector<double>>& profile,
                                      std::mt19937& rng) {
    std::vector<double> probs;
    double total = 0.0;

    for (size_t i = 0; i <= text.size() - k; ++i) {
        double p = probability(text.substr(i, k), profile);
        probs.push_back(p);
        total += p;
    }

    // Нормализуем и выбираем случайно
    std::uniform_real_distribution<double> dist(0.0, total);
    double r = dist(rng);
    double cumulative = 0.0;

    for (size_t i = 0; i < probs.size(); ++i) {
        cumulative += probs[i];
        if (cumulative >= r)
            return text.substr(i, k);
    }
    return text.substr(text.size() - k, k);
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

// Одна итерация Gibbs Sampler
std::vector<std::string> gibbsSamplerOnce(const std::vector<std::string>& dna,
                                           size_t k, size_t t, size_t N,
                                           std::mt19937& rng) {
    // Случайные начальные мотивы
    std::vector<std::string> motifs;
    for (size_t i = 0; i < t; ++i) {
        std::uniform_int_distribution<size_t> dist(0, dna[i].size() - k);
        motifs.push_back(dna[i].substr(dist(rng), k));
    }

    std::vector<std::string> bestMotifs = motifs;

    for (size_t iter = 0; iter < N; ++iter) {
        // Случайно выбираем строку для замены
        std::uniform_int_distribution<size_t> dist(0, t - 1);
        size_t i = dist(rng);

        // Строим профиль без i-й строки
        auto profile = buildProfileExcluding(motifs, i);

        // Выбираем новый k-мер из i-й строки пропорционально вероятностям
        motifs[i] = profileRandomlyGenerated(dna[i], k, profile, rng);

        if (score(motifs) < score(bestMotifs))
            bestMotifs = motifs;
    }
    return bestMotifs;
}

// Запускаем 20 раз и берём лучший результат
std::vector<std::string> gibbsSampler(const std::vector<std::string>& dna,
                                       size_t k, size_t t, size_t N) {
    std::mt19937 rng(42);
    std::vector<std::string> bestMotifs;
    int bestScore = INT_MAX;

    for (int run = 0; run < 20; ++run) {
        auto motifs = gibbsSamplerOnce(dna, k, t, N, rng);
        int s = score(motifs);
        if (s < bestScore) {
            bestScore = s;
            bestMotifs = motifs;
        }
    }
    return bestMotifs;
}

int main(int argc, char* argv[]) {
    size_t k = 8, t = 5, N = 100;
    std::vector<std::string> dna = {
        "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
        "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
        "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
        "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
        "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
    };

    if (argc >= 4) {
        k = std::stoul(argv[1]);
        t = std::stoul(argv[2]);
        N = std::stoul(argv[3]);
        dna.clear();
        for (int i = 4; i < argc; ++i)
            dna.push_back(argv[i]);
    } else {
        std::cout << "No arguments provided. Using sample dataset...\n";
    }

    auto motifs = gibbsSampler(dna, k, t, N);
    for (const auto& m : motifs)
        std::cout << m << std::endl;

    return 0;
}
