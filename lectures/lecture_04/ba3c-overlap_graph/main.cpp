#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

struct OverlapGraphInput {
    std::vector<std::string> kmers;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

std::map<std::string, std::vector<std::string>> solve(const std::vector<std::string>& kmers) {
    std::map<std::string, std::vector<std::string>> adjacency;

    if (kmers.empty())
        return adjacency;

    size_t k = kmers[0].size();

    for (const auto& a : kmers) {
        std::string suffix = a.substr(1);
        for (const auto& b : kmers) {
            std::string prefix = b.substr(0, k - 1);
            if (suffix == prefix) {
                adjacency[a].push_back(b);
            }
        }
    }

    return adjacency;
}

OverlapGraphInput sampleInput() {
    return {
        {"ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT"}
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

OverlapGraphInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    OverlapGraphInput data{};
    std::string kmer;
    while (input >> kmer)
        data.kmers.push_back(kmer);

    return data;
}

OverlapGraphInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.empty())
        throw std::invalid_argument("Expected at least one k-mer as argument");

    OverlapGraphInput data{};
    data.kmers = positional;
    return data;
}

void validateInput(const OverlapGraphInput& data) {
    if (data.kmers.empty())
        throw std::invalid_argument("At least one k-mer is required");

    size_t k = data.kmers[0].size();
    if (k < 2)
        throw std::invalid_argument("k-mers must have length >= 2");

    for (size_t i = 1; i < data.kmers.size(); ++i) {
        if (data.kmers[i].size() != k) {
            throw std::invalid_argument(
                "All k-mers must have the same length, but k-mer #" +
                std::to_string(i + 1) + " has length " +
                std::to_string(data.kmers[i].size()) + " (expected " +
                std::to_string(k) + ")"
            );
        }
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "overlap_graph_output.txt";
}

std::string formatResult(const std::map<std::string, std::vector<std::string>>& adjacency) {
    std::string output;
    for (const auto& [node, targets] : adjacency) {
        output += node + " -> ";
        for (size_t i = 0; i < targets.size(); ++i) {
            if (i > 0)
                output += ",";
            output += targets[i];
        }
        output += '\n';
    }
    if (!output.empty() && output.back() == '\n')
        output.pop_back();
    return output;
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
        OverlapGraphInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        auto result = solve(input.kmers);
        std::string formatted = formatResult(result);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(formatted, outputPath);

        std::cout << formatted << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3c-overlap_graph [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3c-overlap_graph [--output output.txt] kmer1 kmer2 ...\n";
        std::cerr << "  ./ba3c-overlap_graph input.txt [--output output.txt]\n";
        return 1;
    }
}
