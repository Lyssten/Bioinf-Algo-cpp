#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct PairedRead {
    std::string first;
    std::string second;
};

struct InputData {
    size_t k;
    size_t d;
    std::vector<PairedRead> reads;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

std::string solve(const InputData& input) {
    std::string firstStr = input.reads[0].first;
    std::string secondStr = input.reads[0].second;

    for (size_t i = 1; i < input.reads.size(); ++i) {
        firstStr += input.reads[i].first.back();
        secondStr += input.reads[i].second.back();
    }

    size_t k = input.k;
    size_t d = input.d;
    size_t overlap = firstStr.size() - k - d;

    for (size_t i = 0; i < overlap; ++i) {
        if (firstStr[k + d + i] != secondStr[i]) {
            throw std::runtime_error(
                "Mismatch at overlap position " + std::to_string(i) +
                ": first='" + std::string(1, firstStr[k + d + i]) +
                "' second='" + std::string(1, secondStr[i]) + "'");
        }
    }

    return firstStr + secondStr.substr(overlap);
}

InputData sampleInput() {
    return {4, 2, {
        {"GACC", "GCGC"},
        {"ACCG", "CGCC"},
        {"CCGA", "GCCG"},
        {"CGAG", "CCGG"},
        {"GAGC", "CGGA"}
    }};
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
        throw std::invalid_argument("Use either an input file or positional arguments, not both");

    return options;
}

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

InputData parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    InputData data{};
    if (!(input >> data.k >> data.d))
        throw std::runtime_error("Expected first line: k d");

    std::string line;
    std::getline(input, line); // consume rest of first line
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty()) continue;

        size_t sep = line.find('|');
        if (sep == std::string::npos)
            throw std::runtime_error("Expected format 'first|second': " + line);

        data.reads.push_back({line.substr(0, sep), line.substr(sep + 1)});
    }

    return data;
}

void validateInput(const InputData& data) {
    if (data.reads.empty())
        throw std::invalid_argument("At least one paired read is required");
    if (data.k == 0)
        throw std::invalid_argument("k must be greater than 0");

    for (const auto& read : data.reads) {
        if (read.first.size() != data.k || read.second.size() != data.k)
            throw std::invalid_argument("All reads must have length k");
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "gapped_genome_path_output.txt";
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
            throw std::invalid_argument("Positional arguments not supported; use --input <file>");
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
        std::cerr << "  ./ba3l-gapped_genome_path [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3l-gapped_genome_path input.txt\n";
        return 1;
    }
}
