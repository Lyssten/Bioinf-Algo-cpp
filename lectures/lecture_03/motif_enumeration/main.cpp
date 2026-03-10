#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct MotifEnumerationInput {
    size_t k;
    size_t d;
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

// Генерируем все k-меры с расстоянием Хэмминга <= d от pattern
void generateNeighbors(const std::string& pattern, size_t d, std::set<std::string>& result) {
    auto generate = [&](auto&& self, size_t pos, std::string& current, size_t mismatches) -> void {
        if (mismatches > d)
            return;
        if (pos == current.size()) {
            result.insert(current);
            return;
        }

        char original = pattern[pos];
        for (char nucleotide : {'A', 'C', 'G', 'T'}) {
            current[pos] = nucleotide;
            size_t newMismatches = mismatches + (nucleotide != original ? 1 : 0);
            if (newMismatches <= d)
                self(self, pos + 1, current, newMismatches);
        }
        current[pos] = original;
    };

    std::string mutablePattern = pattern;
    generate(generate, 0, mutablePattern, 0);
}

bool appearsWithMismatches(const std::string& text, const std::string& pattern, size_t d) {
    for (size_t i = 0; i <= text.size() - pattern.size(); ++i) {
        if (hammingDistance(text.substr(i, pattern.size()), pattern) <= static_cast<int>(d))
            return true;
    }
    return false;
}

std::set<std::string> motifEnumeration(const std::vector<std::string>& dna, size_t k, size_t d) {
    std::set<std::string> motifs;

    for (size_t i = 0; i <= dna[0].size() - k; ++i) {
        std::string kmer = dna[0].substr(i, k);
        std::set<std::string> neighbors;
        generateNeighbors(kmer, d, neighbors);

        for (const auto& candidate : neighbors) {
            bool inAll = true;
            for (size_t j = 1; j < dna.size(); ++j) {
                if (!appearsWithMismatches(dna[j], candidate, d)) {
                    inAll = false;
                    break;
                }
            }

            if (inAll)
                motifs.insert(candidate);
        }
    }

    return motifs;
}

MotifEnumerationInput sampleInput() {
    return {
        3,
        1,
        {"ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"}
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

MotifEnumerationInput parseInputFile(const std::string& path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open input file: " + path);

    MotifEnumerationInput data{};
    if (!(input >> data.k >> data.d))
        throw std::runtime_error("Expected first line in Rosalind format: k d");

    std::string dnaString;
    while (input >> dnaString)
        data.dna.push_back(dnaString);

    return data;
}

MotifEnumerationInput parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() < 3)
        throw std::invalid_argument("Expected Rosalind-style arguments: k d dna1 dna2 ...");

    MotifEnumerationInput data{};
    data.k = parseSizeT(positional[0], "k");
    data.d = parseSizeT(positional[1], "d");
    data.dna.assign(positional.begin() + 2, positional.end());
    return data;
}

void validateInput(const MotifEnumerationInput& data) {
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

std::string joinMotifs(const std::set<std::string>& motifs) {
    std::ostringstream output;
    bool first = true;

    for (const auto& motif : motifs) {
        if (!first)
            output << ' ';
        output << motif;
        first = false;
    }

    return output.str();
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "motif_enumeration_output.txt";
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
        MotifEnumerationInput input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            input = parsePositionalArgs(options.positional);
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        std::string result = joinMotifs(motifEnumeration(input.dna, input.k, input.d));
        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./motif_enumeration [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./motif_enumeration [--output output.txt] k d dna1 dna2 ...\n";
        std::cerr << "  ./motif_enumeration input.txt [--output output.txt]\n";
        return 1;
    }
}
