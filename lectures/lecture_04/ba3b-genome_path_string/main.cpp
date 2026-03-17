#include <fstream>
#include <iostream>
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

std::string solve(const InputData& data) {
    std::string result = data.kmers[0];
    for (size_t i = 1; i < data.kmers.size(); ++i) {
        result += data.kmers[i].back();
    }
    return result;
}

InputData sampleInput() {
    return {
        {"ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"}
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

InputData parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    InputData data{};
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty()) {
            // Strip trailing \r for Windows line endings
            if (line.back() == '\r')
                line.pop_back();
            if (!line.empty())
                data.kmers.push_back(line);
        }
    }

    return data;
}

InputData parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.empty())
        throw std::invalid_argument("Expected at least one k-mer as argument");

    InputData data{};
    data.kmers = positional;
    return data;
}

void validateInput(const InputData& data) {
    if (data.kmers.empty())
        throw std::invalid_argument("At least one k-mer is required");

    size_t k = data.kmers[0].size();
    if (k == 0)
        throw std::invalid_argument("k-mers must be non-empty");

    for (size_t i = 1; i < data.kmers.size(); ++i) {
        if (data.kmers[i].size() != k) {
            throw std::invalid_argument(
                "k-mer #" + std::to_string(i + 1) + " has length " +
                std::to_string(data.kmers[i].size()) + ", expected " + std::to_string(k)
            );
        }
    }

    for (size_t i = 0; i + 1 < data.kmers.size(); ++i) {
        std::string suffix = data.kmers[i].substr(1);
        std::string prefix = data.kmers[i + 1].substr(0, k - 1);
        if (suffix != prefix) {
            throw std::invalid_argument(
                "Overlap mismatch between k-mer #" + std::to_string(i + 1) +
                " and #" + std::to_string(i + 2) +
                ": suffix \"" + suffix + "\" != prefix \"" + prefix + "\""
            );
        }
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "genome_path_string_output.txt";
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

        std::string result = solve(input);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3b-genome_path_string [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3b-genome_path_string [--output output.txt] kmer1 kmer2 ...\n";
        std::cerr << "  ./ba3b-genome_path_string input.txt [--output output.txt]\n";
        return 1;
    }
}
