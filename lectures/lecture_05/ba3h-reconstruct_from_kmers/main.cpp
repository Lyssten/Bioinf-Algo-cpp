#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

struct InputData {
    size_t k;
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

std::string spellString(const std::vector<std::string>& path) {
    if (path.empty()) return "";
    std::string result = path[0];
    for (size_t i = 1; i < path.size(); ++i) {
        result += path[i].back();
    }
    return result;
}

std::string solve(const InputData& input) {
    auto graph = buildDeBruijn(input.kmers);
    auto path = eulerianPath(graph);
    return spellString(path);
}

InputData sampleInput() {
    return {4, {"CTTA", "ACCA", "TACC", "GGCT", "GCTT", "TTAC"}};
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
    std::string line;
    std::getline(input, line);
    data.k = std::stoul(trim(line));

    while (std::getline(input, line)) {
        line = trim(line);
        if (!line.empty())
            data.kmers.push_back(line);
    }

    return data;
}

void validateInput(const InputData& data) {
    if (data.kmers.empty())
        throw std::invalid_argument("At least one k-mer is required");
    if (data.k < 2)
        throw std::invalid_argument("k must be at least 2");
    for (const auto& kmer : data.kmers) {
        if (kmer.size() != data.k)
            throw std::invalid_argument("All k-mers must have length k");
    }
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "reconstruct_from_kmers_output.txt";
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
        std::cerr << "  ./ba3h-reconstruct_from_kmers [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3h-reconstruct_from_kmers input.txt\n";
        return 1;
    }
}
