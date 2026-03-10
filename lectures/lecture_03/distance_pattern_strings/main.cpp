#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct DistanceInput {
    std::string pattern;
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

DistanceInput sampleInput() {
    return {
        "AAA",
        {"TTACCTTAAC", "GATATCTGTC", "ACGGCGTTCG", "CCCTAAAGAG", "CGTCAGAGGT"}
    };
}

bool isReadableFile(const std::string& path) {
    std::ifstream input(path);
    return input.good();
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

std::vector<std::string> splitWhitespaceTokens(const std::vector<std::string>& tokens) {
    std::vector<std::string> result;

    for (const auto& token : tokens) {
        std::istringstream stream(token);
        std::string value;
        while (stream >> value)
            result.push_back(value);
    }

    return result;
}

DistanceInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    DistanceInput data{};
    if (!(input >> data.pattern))
        throw std::runtime_error("Expected Rosalind format: pattern followed by DNA strings");

    std::string dnaString;
    while (input >> dnaString)
        data.dna.push_back(dnaString);

    return data;
}

DistanceInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() < 2)
        throw std::invalid_argument("Expected Rosalind-style arguments: pattern dna1 dna2 ...");

    DistanceInput data{};
    data.pattern = positional[0];
    data.dna = splitWhitespaceTokens(std::vector<std::string>(positional.begin() + 1, positional.end()));
    return data;
}

void validateInput(const DistanceInput& data) {
    if (data.pattern.empty())
        throw std::invalid_argument("pattern must not be empty");
    if (data.dna.empty())
        throw std::invalid_argument("At least one DNA string is required");

    for (size_t i = 0; i < data.dna.size(); ++i) {
        if (data.dna[i].size() < data.pattern.size()) {
            throw std::invalid_argument(
                "DNA string #" + std::to_string(i + 1) + " is shorter than the pattern"
            );
        }
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "distance_pattern_strings_output.txt";
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
        DistanceInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        std::string result = std::to_string(distanceBetweenPatternAndStrings(input.pattern, input.dna));
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./distance_pattern_strings [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./distance_pattern_strings [--output output.txt] pattern dna1 dna2 ...\n";
        std::cerr << "  ./distance_pattern_strings input.txt [--output output.txt]\n";
        return 1;
    }
}
