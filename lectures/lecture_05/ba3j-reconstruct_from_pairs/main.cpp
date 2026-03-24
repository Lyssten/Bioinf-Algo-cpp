#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
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

// Node in paired de Bruijn graph: "prefix1|prefix2"
using Graph = std::map<std::string, std::vector<std::string>>;

std::string makeNode(const std::string& a, const std::string& b) {
    return a + "|" + b;
}

std::pair<std::string, std::string> splitNode(const std::string& node) {
    size_t sep = node.find('|');
    return {node.substr(0, sep), node.substr(sep + 1)};
}

Graph buildPairedDeBruijn(const std::vector<PairedRead>& reads) {
    Graph graph;
    for (const auto& read : reads) {
        std::string prefixNode = makeNode(
            read.first.substr(0, read.first.size() - 1),
            read.second.substr(0, read.second.size() - 1));
        std::string suffixNode = makeNode(
            read.first.substr(1),
            read.second.substr(1));
        graph[prefixNode].push_back(suffixNode);
    }
    return graph;
}

std::vector<std::string> eulerianCycle(Graph& graph, const std::string& startNode) {
    std::map<std::string, size_t> edgeIndex;
    for (auto& [v, _] : graph)
        edgeIndex[v] = 0;

    std::vector<std::string> stack = {startNode};
    std::vector<std::string> cycle;

    while (!stack.empty()) {
        auto& v = stack.back();
        if (edgeIndex.count(v) && edgeIndex[v] < graph[v].size()) {
            std::string next = graph[v][edgeIndex[v]++];
            stack.push_back(next);
        } else {
            cycle.push_back(v);
            stack.pop_back();
        }
    }

    std::reverse(cycle.begin(), cycle.end());
    return cycle;
}

std::vector<std::string> eulerianPath(Graph& graph) {
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

    std::string startNode, endNode;
    for (auto& v : allNodes) {
        int diff = outDeg[v] - inDeg[v];
        if (diff == 1) startNode = v;
        else if (diff == -1) endNode = v;
    }

    if (startNode.empty()) {
        return eulerianCycle(graph, graph.begin()->first);
    }

    graph[endNode].push_back(startNode);
    auto cycle = eulerianCycle(graph, startNode);
    graph[endNode].pop_back();

    for (size_t i = 0; i + 1 < cycle.size(); ++i) {
        if (cycle[i] == endNode && cycle[i + 1] == startNode) {
            std::vector<std::string> path;
            for (size_t j = i + 1; j < cycle.size(); ++j)
                path.push_back(cycle[j]);
            for (size_t j = 1; j <= i; ++j)
                path.push_back(cycle[j]);
            return path;
        }
    }

    return cycle;
}

std::string solve(const InputData& input) {
    auto graph = buildPairedDeBruijn(input.reads);
    auto path = eulerianPath(graph);

    // Spell strings from first and second components
    std::string firstStr, secondStr;
    for (size_t i = 0; i < path.size(); ++i) {
        auto [a, b] = splitNode(path[i]);
        if (i == 0) {
            firstStr = a;
            secondStr = b;
        } else {
            firstStr += a.back();
            secondStr += b.back();
        }
    }

    // Combine: secondStr starts at position k+d in the result
    size_t k = input.k;
    size_t d = input.d;
    size_t overlap = firstStr.size() - k - d;

    for (size_t i = 0; i < overlap; ++i) {
        if (firstStr[k + d + i] != secondStr[i]) {
            throw std::runtime_error("Mismatch at overlap position " + std::to_string(i));
        }
    }

    return firstStr + secondStr.substr(overlap);
}

InputData sampleInput() {
    return {4, 2, {
        {"GAGA", "TTGA"},
        {"TCGT", "GATG"},
        {"CGTG", "ATGT"},
        {"TGGT", "TGAG"},
        {"GTGA", "TGTT"},
        {"GTGG", "GTGA"},
        {"TGAG", "GTTG"},
        {"GGTC", "GAGA"},
        {"GTCG", "AGAT"}
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
    if (data.k < 2)
        throw std::invalid_argument("k must be at least 2");
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
    return "reconstruct_from_pairs_output.txt";
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
        std::cerr << "  ./ba3j-reconstruct_from_pairs [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3j-reconstruct_from_pairs input.txt\n";
        return 1;
    }
}
