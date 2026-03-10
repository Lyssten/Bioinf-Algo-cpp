#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct GibbsInput {
    size_t k;
    size_t t;
    size_t N;
    std::vector<std::string> dna;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

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

GibbsInput sampleInput() {
    return {
        8,
        5,
        100,
        {
            "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
            "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
            "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
            "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
            "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
        }
    };
}

bool isReadableFile(const std::string& path) {
    std::ifstream input(path);
    return input.good();
}

size_t parseSizeT(const std::string& value, const std::string& name) {
    size_t pos = 0;
    unsigned long parsed = 0;

    try {
        parsed = std::stoul(value, &pos);
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid numeric value for " + name + ": " + value);
    }

    if (pos != value.size())
        throw std::invalid_argument("Invalid numeric value for " + name + ": " + value);

    return static_cast<size_t>(parsed);
}

CliOptions parseCli(int argc, char* argv[]) {
    CliOptions options;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--input" || arg == "-i") {
            if (i + 1 >= argc)
                throw std::invalid_argument("Missing value after --input");
            options.inputPath = argv[++i];
            continue;
        }

        if (arg == "--output" || arg == "-o") {
            if (i + 1 >= argc)
                throw std::invalid_argument("Missing value after --output");
            options.outputPath = argv[++i];
            continue;
        }

        options.positional.push_back(arg);
    }

    if (options.inputPath.empty() && options.positional.size() == 1 &&
        isReadableFile(options.positional[0])) {
        options.inputPath = options.positional[0];
        options.positional.clear();
    }

    if (!options.inputPath.empty() && !options.positional.empty())
        throw std::invalid_argument("Use either Rosalind-style positional arguments or an input file, not both");

    return options;
}

GibbsInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    GibbsInput data{};
    if (!(input >> data.k >> data.t >> data.N))
        throw std::runtime_error("Expected first line in Rosalind format: k t N");

    std::string dnaString;
    while (input >> dnaString)
        data.dna.push_back(dnaString);

    return data;
}

GibbsInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() < 4)
        throw std::invalid_argument("Expected Rosalind-style arguments: k t N dna1 dna2 ...");

    GibbsInput data{};
    data.k = parseSizeT(positional[0], "k");
    data.t = parseSizeT(positional[1], "t");
    data.N = parseSizeT(positional[2], "N");
    data.dna.assign(positional.begin() + 3, positional.end());
    return data;
}

void validateInput(const GibbsInput& data) {
    if (data.k == 0)
        throw std::invalid_argument("k must be greater than 0");
    if (data.t == 0)
        throw std::invalid_argument("t must be greater than 0");
    if (data.N == 0)
        throw std::invalid_argument("N must be greater than 0");
    if (data.dna.empty())
        throw std::invalid_argument("At least one DNA string is required");
    if (data.dna.size() != data.t) {
        throw std::invalid_argument(
            "The number of DNA strings does not match t: expected " + std::to_string(data.t) +
            ", got " + std::to_string(data.dna.size())
        );
    }

    for (size_t i = 0; i < data.dna.size(); ++i) {
        if (data.dna[i].size() < data.k) {
            throw std::invalid_argument(
                "DNA string #" + std::to_string(i + 1) + " is shorter than k"
            );
        }
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "gibbs_sampler_output.txt";
}

void writeOutputFile(const std::vector<std::string>& motifs, const std::string& outputPath) {
    std::ofstream output(outputPath);
    if (!output)
        throw std::runtime_error("Cannot open output file: " + outputPath);

    for (const auto& motif : motifs)
        output << motif << '\n';
}

int main(int argc, char* argv[]) {
    try {
        CliOptions options = parseCli(argc, argv);
        GibbsInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        auto motifs = gibbsSampler(input.dna, input.k, input.t, input.N);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(motifs, outputPath);

        for (const auto& motif : motifs)
            std::cout << motif << '\n';

        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./gibbs_sampler [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./gibbs_sampler [--output output.txt] k t N dna1 dna2 ... dna_t\n";
        std::cerr << "  ./gibbs_sampler input.txt [--output output.txt]\n";
        return 1;
    }
}
