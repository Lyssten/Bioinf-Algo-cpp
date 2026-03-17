#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct InputData {
    size_t k;
    std::string text;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

std::vector<std::string> solve(const InputData& input) {
    std::vector<std::string> kmers;

    for (size_t i = 0; i + input.k <= input.text.size(); ++i) {
        kmers.push_back(input.text.substr(i, input.k));
    }

    std::sort(kmers.begin(), kmers.end());
    return kmers;
}

InputData sampleInput() {
    return {5, "CAATCCAAC"};
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

InputData parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    InputData data{};
    if (!(input >> data.k))
        throw std::runtime_error("Expected first line in Rosalind format: k");

    if (!(input >> data.text))
        throw std::runtime_error("Expected second line in Rosalind format: Text");

    return data;
}

InputData parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() != 2)
        throw std::invalid_argument("Expected Rosalind-style arguments: k Text");

    InputData data{};
    data.k = parseSizeT(positional[0], "k");
    data.text = positional[1];
    return data;
}

void validateInput(const InputData& data) {
    if (data.k == 0)
        throw std::invalid_argument("k must be greater than 0");
    if (data.text.empty())
        throw std::invalid_argument("Text must not be empty");
    if (data.text.size() < data.k)
        throw std::invalid_argument("Text is shorter than k");
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "string_composition_output.txt";
}

void writeOutputFile(const std::vector<std::string>& result, const std::string& outputPath) {
    std::ofstream output(outputPath);
    if (!output)
        throw std::runtime_error("Cannot open output file: " + outputPath);

    for (const auto& kmer : result)
        output << kmer << '\n';
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

        auto result = solve(input);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        for (const auto& kmer : result)
            std::cout << kmer << '\n';

        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3a-string_composition [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3a-string_composition [--output output.txt] k Text\n";
        std::cerr << "  ./ba3a-string_composition input.txt [--output output.txt]\n";
        return 1;
    }
}
