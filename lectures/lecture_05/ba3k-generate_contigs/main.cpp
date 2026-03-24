#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
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

using Graph = std::map<std::string, std::vector<std::string>>;

Graph buildDeBruijn(const std::vector<std::string>& kmers) {
    Graph graph;
    for (const auto& kmer : kmers) {
        std::string prefix = kmer.substr(0, kmer.size() - 1);
        std::string suffix = kmer.substr(1);
        graph[prefix].push_back(suffix);
    }
    return graph;
}

std::string spellPath(const std::vector<std::string>& path) {
    if (path.empty()) return "";
    std::string result = path[0];
    for (size_t i = 1; i < path.size(); ++i)
        result += path[i].back();
    return result;
}

std::vector<std::string> generateContigs(const std::vector<std::string>& kmers) {
    auto graph = buildDeBruijn(kmers);

    // Compute degrees
    std::map<std::string, int> inDeg, outDeg;
    std::set<std::string> allNodes;

    for (auto& [v, neighbors] : graph) {
        allNodes.insert(v);
        outDeg[v] = static_cast<int>(neighbors.size());
        for (auto& u : neighbors) {
            inDeg[u]++;
            allNodes.insert(u);
        }
    }

    for (auto& v : allNodes) {
        if (!inDeg.count(v)) inDeg[v] = 0;
        if (!outDeg.count(v)) outDeg[v] = 0;
    }

    auto is1in1out = [&](const std::string& v) {
        return inDeg[v] == 1 && outDeg[v] == 1;
    };

    std::vector<std::string> contigs;
    std::set<std::string> visited;

    // Phase 1: paths from non-1-in-1-out nodes
    for (auto& v : allNodes) {
        if (!is1in1out(v) && outDeg[v] > 0) {
            auto it = graph.find(v);
            if (it == graph.end()) continue;

            for (auto& w : it->second) {
                std::vector<std::string> path = {v, w};
                visited.insert(v);
                visited.insert(w);

                std::string current = w;
                while (is1in1out(current)) {
                    auto cIt = graph.find(current);
                    if (cIt == graph.end() || cIt->second.empty()) break;
                    current = cIt->second[0];
                    path.push_back(current);
                    visited.insert(current);
                }

                contigs.push_back(spellPath(path));
            }
        }
    }

    // Phase 2: isolated cycles
    for (auto& v : allNodes) {
        if (is1in1out(v) && !visited.count(v)) {
            auto it = graph.find(v);
            if (it == graph.end()) continue;

            std::vector<std::string> cycle = {v};
            visited.insert(v);
            std::string current = it->second[0];

            while (current != v) {
                cycle.push_back(current);
                visited.insert(current);
                auto cIt = graph.find(current);
                if (cIt == graph.end()) break;
                current = cIt->second[0];
            }

            cycle.push_back(v);
            contigs.push_back(spellPath(cycle));
        }
    }

    return contigs;
}

InputData sampleInput() {
    return {{"ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA"}};
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
    std::string kmer;
    while (input >> kmer) {
        kmer = trim(kmer);
        if (!kmer.empty())
            data.kmers.push_back(kmer);
    }

    return data;
}

InputData parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.empty())
        throw std::invalid_argument("Expected k-mers as arguments");
    return {positional};
}

void validateInput(const InputData& data) {
    if (data.kmers.empty())
        throw std::invalid_argument("At least one k-mer is required");
    size_t k = data.kmers[0].size();
    if (k < 2)
        throw std::invalid_argument("k-mers must have length >= 2");
    for (const auto& kmer : data.kmers) {
        if (kmer.size() != k)
            throw std::invalid_argument("All k-mers must have the same length");
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "generate_contigs_output.txt";
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

        auto contigs = generateContigs(input.kmers);

        std::string result;
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (i > 0) result += ' ';
            result += contigs[i];
        }

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3k-generate_contigs [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3k-generate_contigs input.txt\n";
        std::cerr << "  ./ba3k-generate_contigs kmer1 kmer2 ...\n";
        return 1;
    }
}
