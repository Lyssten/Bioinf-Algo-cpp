#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

struct RandomizedMotifSearchInput {
    size_t k;
    size_t t;
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

std::vector<std::vector<double>> buildProfileWithPseudocounts(const std::vector<std::string>& motifs) {
    size_t k = motifs[0].size();
    size_t t = motifs.size();
    std::vector<std::vector<double>> profile(4, std::vector<double>(k, 1.0));

    for (const auto& motif : motifs) {
        for (size_t j = 0; j < k; ++j)
            ++profile[nucleotideIndex(motif[j])][j];
    }

    double total = static_cast<double>(t) + 4.0;
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < k; ++j)
            profile[i][j] /= total;
    }

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
        for (int value : count)
            maxCount = std::max(maxCount, value);

        totalScore += static_cast<int>(motifs.size()) - maxCount;
    }

    return totalScore;
}

std::vector<std::string> randomizedMotifSearchOnce(const std::vector<std::string>& dna,
                                                   size_t k, size_t t, std::mt19937& rng) {
    std::vector<std::string> motifs;
    for (size_t i = 0; i < t; ++i) {
        std::uniform_int_distribution<size_t> dist(0, dna[i].size() - k);
        motifs.push_back(dna[i].substr(dist(rng), k));
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

std::vector<std::string> randomizedMotifSearch(const std::vector<std::string>& dna,
                                               size_t k, size_t t) {
    std::mt19937 rng(42);
    std::vector<std::string> bestMotifs;
    int bestScore = INT_MAX;

    for (int iter = 0; iter < 1000; ++iter) {
        auto motifs = randomizedMotifSearchOnce(dna, k, t, rng);
        int currentScore = score(motifs);
        if (currentScore < bestScore) {
            bestScore = currentScore;
            bestMotifs = motifs;
        }
    }

    return bestMotifs;
}

RandomizedMotifSearchInput sampleInput() {
    return {
        8,
        5,
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

RandomizedMotifSearchInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    RandomizedMotifSearchInput data{};
    if (!(input >> data.k >> data.t))
        throw std::runtime_error("Expected first line in Rosalind format: k t");

    std::string dnaString;
    while (input >> dnaString)
        data.dna.push_back(dnaString);

    return data;
}

RandomizedMotifSearchInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() < 3)
        throw std::invalid_argument("Expected Rosalind-style arguments: k t dna1 dna2 ...");

    RandomizedMotifSearchInput data{};
    data.k = parseSizeT(positional[0], "k");
    data.t = parseSizeT(positional[1], "t");
    data.dna.assign(positional.begin() + 2, positional.end());
    return data;
}

void validateInput(const RandomizedMotifSearchInput& data) {
    if (data.k == 0)
        throw std::invalid_argument("k must be greater than 0");
    if (data.t == 0)
        throw std::invalid_argument("t must be greater than 0");
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
    return "randomized_motif_search_output.txt";
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
        RandomizedMotifSearchInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        auto motifs = randomizedMotifSearch(input.dna, input.k, input.t);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(motifs, outputPath);

        for (const auto& motif : motifs)
            std::cout << motif << '\n';

        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./randomized_motif_search [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./randomized_motif_search [--output output.txt] k t dna1 dna2 ... dna_t\n";
        std::cerr << "  ./randomized_motif_search input.txt [--output output.txt]\n";
        return 1;
    }
}
