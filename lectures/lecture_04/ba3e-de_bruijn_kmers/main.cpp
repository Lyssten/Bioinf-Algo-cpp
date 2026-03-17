#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

struct InputData {
    std::vector<std::string> kmers;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

using AdjList = std::map<std::string, std::vector<std::string>>;

AdjList solve(const InputData& input) {
    AdjList graph;

    for (const auto& kmer : input.kmers) {
        std::string prefix = kmer.substr(0, kmer.size() - 1);
        std::string suffix = kmer.substr(1);
        graph[prefix].push_back(suffix);
    }

    for (auto& [node, targets] : graph) {
        std::sort(targets.begin(), targets.end());
    }

    return graph;
}

InputData sampleInput() {
    return {{"GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG"}};
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

InputData parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    InputData data{};
    std::string kmer;
    while (input >> kmer)
        data.kmers.push_back(kmer);

    return data;
}

InputData parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.empty())
        throw std::invalid_argument("Expected Rosalind-style arguments: kmer1 kmer2 ...");

    InputData data{};
    data.kmers = positional;
    return data;
}

void validateInput(const InputData& data) {
    if (data.kmers.empty())
        throw std::invalid_argument("At least one k-mer is required");

    size_t k = data.kmers[0].size();
    if (k < 2)
        throw std::invalid_argument("k-mers must have length >= 2");

    for (size_t i = 1; i < data.kmers.size(); ++i) {
        if (data.kmers[i].size() != k)
            throw std::invalid_argument("All k-mers must have the same length");
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "de_bruijn_kmers_output.txt";
}

std::string formatResult(const AdjList& graph) {
    std::string result;
    for (const auto& [node, targets] : graph) {
        result += node + " -> ";
        for (size_t i = 0; i < targets.size(); ++i) {
            if (i > 0)
                result += ',';
            result += targets[i];
        }
        result += '\n';
    }
    if (!result.empty() && result.back() == '\n')
        result.pop_back();
    return result;
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
        InputData input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        auto graph = solve(input);
        std::string result = formatResult(graph);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3e-de_bruijn_kmers [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3e-de_bruijn_kmers [--output output.txt] kmer1 kmer2 ...\n";
        std::cerr << "  ./ba3e-de_bruijn_kmers input.txt [--output output.txt]\n";
        return 1;
    }
}
