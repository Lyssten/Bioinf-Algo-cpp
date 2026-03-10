#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct ProfileMostProbableInput {
    std::string text;
    size_t k;
    std::vector<std::vector<double>> profile;
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

ProfileMostProbableInput sampleInput() {
    return {
        "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",
        5,
        {
            {0.2, 0.2, 0.3, 0.2, 0.3},
            {0.4, 0.3, 0.1, 0.5, 0.1},
            {0.3, 0.3, 0.5, 0.2, 0.4},
            {0.1, 0.2, 0.1, 0.1, 0.2}
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

double parseDouble(const std::string& value, const std::string& name) {
    size_t pos = 0;
    double parsed = 0.0;

    try {
        parsed = std::stod(value, &pos);
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid numeric value for " + name + ": " + value);
    }

    if (pos != value.size())
        throw std::invalid_argument("Invalid numeric value for " + name + ": " + value);

    return parsed;
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

ProfileMostProbableInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    ProfileMostProbableInput data{};
    if (!(input >> data.text >> data.k))
        throw std::runtime_error("Expected Rosalind format: text, k, then 4 profile rows");

    data.profile.assign(4, std::vector<double>(data.k, 0.0));
    for (size_t row = 0; row < 4; ++row) {
        for (size_t col = 0; col < data.k; ++col) {
            if (!(input >> data.profile[row][col])) {
                throw std::runtime_error(
                    "Expected 4 profile rows with exactly " + std::to_string(data.k) + " columns"
                );
            }
        }
    }

    return data;
}

std::vector<double> parseProfileValues(const std::vector<std::string>& args) {
    std::vector<double> values;

    for (const auto& arg : args) {
        std::istringstream stream(arg);
        std::string token;
        while (stream >> token)
            values.push_back(parseDouble(token, "profile"));
    }

    return values;
}

ProfileMostProbableInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() < 3)
        throw std::invalid_argument("Expected Rosalind-style arguments: text k profile");

    ProfileMostProbableInput data{};
    data.text = positional[0];
    data.k = parseSizeT(positional[1], "k");

    std::vector<double> values = parseProfileValues(
        std::vector<std::string>(positional.begin() + 2, positional.end())
    );

    if (values.size() != 4 * data.k) {
        throw std::invalid_argument(
            "Expected exactly " + std::to_string(4 * data.k) + " profile values after text and k"
        );
    }

    data.profile.assign(4, std::vector<double>(data.k, 0.0));
    size_t index = 0;
    for (size_t row = 0; row < 4; ++row) {
        for (size_t col = 0; col < data.k; ++col)
            data.profile[row][col] = values[index++];
    }

    return data;
}

void validateInput(const ProfileMostProbableInput& data) {
    if (data.k == 0)
        throw std::invalid_argument("k must be greater than 0");
    if (data.text.size() < data.k)
        throw std::invalid_argument("text is shorter than k");
    if (data.profile.size() != 4)
        throw std::invalid_argument("profile must contain exactly 4 rows");

    for (size_t row = 0; row < data.profile.size(); ++row) {
        if (data.profile[row].size() != data.k) {
            throw std::invalid_argument(
                "Profile row #" + std::to_string(row + 1) + " must contain exactly " +
                std::to_string(data.k) + " values"
            );
        }
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "profile_most_probable_output.txt";
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
        ProfileMostProbableInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        std::string result = profileMostProbable(input.text, input.k, input.profile);
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./profile_most_probable [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./profile_most_probable [--output output.txt] text k \"rowA\" \"rowC\" \"rowG\" \"rowT\"\n";
        std::cerr << "  ./profile_most_probable input.txt [--output output.txt]\n";
        return 1;
    }
}
