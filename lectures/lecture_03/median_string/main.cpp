#include <algorithm>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct MedianStringInput {
    size_t k;
    std::vector<std::string> dna;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

int hammingDistance(const std::string& s1, const std::string& s2) {
    int dist = 0;
    for (size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != s2[i])
            ++dist;
    }
    return dist;
}

int distanceBetweenPatternAndStrings(const std::string& pattern, const std::vector<std::string>& dna) {
    size_t k = pattern.size();
    int totalDist = 0;

    for (const auto& text : dna) {
        int minDist = INT_MAX;
        for (size_t i = 0; i <= text.size() - k; ++i) {
            int dist = hammingDistance(pattern, text.substr(i, k));
            minDist = std::min(minDist, dist);
        }
        totalDist += minDist;
    }

    return totalDist;
}

char numberToNucleotide(int n) {
    switch (n) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return '?';
    }
}

std::string numberToPattern(int64_t index, size_t k) {
    std::string pattern(k, 'A');
    for (size_t i = k; i > 0; --i) {
        pattern[i - 1] = numberToNucleotide(static_cast<int>(index % 4));
        index /= 4;
    }
    return pattern;
}

std::string medianString(const std::vector<std::string>& dna, size_t k) {
    int minDist = INT_MAX;
    std::string median;

    int64_t total = 1;
    for (size_t i = 0; i < k; ++i)
        total *= 4;

    for (int64_t i = 0; i < total; ++i) {
        std::string pattern = numberToPattern(i, k);
        int dist = distanceBetweenPatternAndStrings(pattern, dna);
        if (dist < minDist) {
            minDist = dist;
            median = pattern;
        }
    }

    return median;
}

MedianStringInput sampleInput() {
    return {
        3,
        {
            "AAATTGACGCAT",
            "GACGACCACGTT",
            "CGTCAGCGCCTG",
            "GCTGAGCACCGG",
            "AGTACGGGACAG"
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

MedianStringInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    MedianStringInput data{};
    if (!(input >> data.k))
        throw std::runtime_error("Expected first line in Rosalind format: k");

    std::string dnaString;
    while (input >> dnaString)
        data.dna.push_back(dnaString);

    return data;
}

MedianStringInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() < 2)
        throw std::invalid_argument("Expected Rosalind-style arguments: k dna1 dna2 ...");

    MedianStringInput data{};
    data.k = parseSizeT(positional[0], "k");
    data.dna.assign(positional.begin() + 1, positional.end());
    return data;
}

void validateInput(const MedianStringInput& data) {
    if (data.k == 0)
        throw std::invalid_argument("k must be greater than 0");
    if (data.dna.empty())
        throw std::invalid_argument("At least one DNA string is required");

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
    return "median_string_output.txt";
}

void writeOutputFile(const std::string& content, const std::string& outputPath) {
    std::ofstream output(outputPath);
    if (!output)
        throw std::runtime_error("Cannot open output file: " + outputPath);

    output << content << '\n';
}

int main(int argc, char* argv[]) {
    try {
        CliOptions options = parseCli(argc, argv);
        MedianStringInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        std::string result = medianString(input.dna, input.k);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./median_string [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./median_string [--output output.txt] k dna1 dna2 ...\n";
        std::cerr << "  ./median_string input.txt [--output output.txt]\n";
        return 1;
    }
}
