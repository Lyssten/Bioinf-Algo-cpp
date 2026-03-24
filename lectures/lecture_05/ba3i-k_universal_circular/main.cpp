#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

struct InputData {
    size_t k;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

std::string intToBinary(int val, size_t len) {
    std::string s(len, '0');
    for (size_t i = len; i > 0; --i) {
        s[i - 1] = '0' + (val & 1);
        val >>= 1;
    }
    return s;
}

std::string solve(const InputData& input) {
    size_t k = input.k;
    int numNodes = 1 << (k - 1); // 2^(k-1)

    // Build de Bruijn graph: node i has edges to (i*2) % numNodes and (i*2+1) % numNodes
    std::map<int, std::vector<int>> graph;
    for (int i = 0; i < numNodes; ++i) {
        graph[i].push_back((i * 2) % numNodes);
        graph[i].push_back((i * 2 + 1) % numNodes);
    }

    // Eulerian cycle via Hierholzer
    std::map<int, size_t> edgeIndex;
    for (auto& [v, _] : graph)
        edgeIndex[v] = 0;

    std::vector<int> stack = {0};
    std::vector<int> cycle;

    while (!stack.empty()) {
        int v = stack.back();
        if (edgeIndex[v] < graph[v].size()) {
            int u = graph[v][edgeIndex[v]++];
            stack.push_back(u);
        } else {
            cycle.push_back(v);
            stack.pop_back();
        }
    }

    std::reverse(cycle.begin(), cycle.end());

    // Spell circular string: each edge contributes one bit
    // Edge from node v to node u represents the k-mer whose binary is (v shifted left + last bit of u)
    // For circular string, take last bit of each target node
    std::string result;
    for (size_t i = 0; i + 1 < cycle.size(); ++i) {
        result += ('0' + (cycle[i + 1] & 1));
    }

    return result;
}

InputData sampleInput() {
    return {4};
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

    if (options.inputPath.empty() && options.positional.size() == 1) {
        if (isReadableFile(options.positional[0])) {
            options.inputPath = options.positional[0];
            options.positional.clear();
        }
    }

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
    std::string line;
    std::getline(input, line);
    data.k = parseSizeT(trim(line), "k");
    return data;
}

InputData parsePositionalArgs(const std::vector<std::string>& positional) {
    if (positional.size() != 1)
        throw std::invalid_argument("Expected one argument: k");
    return {parseSizeT(positional[0], "k")};
}

void validateInput(const InputData& data) {
    if (data.k < 1)
        throw std::invalid_argument("k must be at least 1");
    if (data.k > 20)
        throw std::invalid_argument("k too large (max 20, would produce 2^20 characters)");
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "k_universal_circular_output.txt";
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
        std::cerr << "  ./ba3i-k_universal_circular [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3i-k_universal_circular k\n";
        return 1;
    }
}
