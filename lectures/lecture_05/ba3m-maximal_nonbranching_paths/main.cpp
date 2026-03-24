#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct InputData {
    std::map<int, std::vector<int>> graph;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

std::vector<std::vector<int>> maximalNonBranchingPaths(
        const std::map<int, std::vector<int>>& graph) {

    std::map<int, int> inDeg, outDeg;
    std::set<int> allNodes;

    for (auto& [v, neighbors] : graph) {
        allNodes.insert(v);
        outDeg[v] = static_cast<int>(neighbors.size());
        for (int u : neighbors) {
            inDeg[u]++;
            allNodes.insert(u);
        }
    }

    // Ensure all nodes have degree entries
    for (int v : allNodes) {
        if (!inDeg.count(v)) inDeg[v] = 0;
        if (!outDeg.count(v)) outDeg[v] = 0;
    }

    auto is1in1out = [&](int v) {
        return inDeg[v] == 1 && outDeg[v] == 1;
    };

    std::vector<std::vector<int>> paths;
    std::set<int> visited;

    // Phase 1: paths starting from non-1-in-1-out nodes
    for (int v : allNodes) {
        if (!is1in1out(v) && outDeg[v] > 0) {
            auto it = graph.find(v);
            if (it == graph.end()) continue;

            for (int w : it->second) {
                std::vector<int> path = {v, w};
                visited.insert(v);
                visited.insert(w);

                while (is1in1out(w)) {
                    auto wIt = graph.find(w);
                    if (wIt == graph.end() || wIt->second.empty()) break;
                    w = wIt->second[0];
                    path.push_back(w);
                    visited.insert(w);
                }

                paths.push_back(path);
            }
        }
    }

    // Phase 2: isolated cycles (all nodes are 1-in-1-out)
    for (int v : allNodes) {
        if (is1in1out(v) && !visited.count(v)) {
            auto it = graph.find(v);
            if (it == graph.end()) continue;

            std::vector<int> cycle = {v};
            visited.insert(v);
            int w = it->second[0];

            while (w != v) {
                cycle.push_back(w);
                visited.insert(w);
                auto wIt = graph.find(w);
                if (wIt == graph.end()) break;
                w = wIt->second[0];
            }

            cycle.push_back(v); // close the cycle
            paths.push_back(cycle);
        }
    }

    return paths;
}

InputData sampleInput() {
    InputData data;
    data.graph[1] = {2};
    data.graph[2] = {3};
    data.graph[3] = {4, 5};
    data.graph[6] = {7};
    data.graph[7] = {6};
    return data;
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

    InputData data;
    std::string line;

    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty()) continue;

        size_t arrowPos = line.find(" -> ");
        if (arrowPos == std::string::npos)
            throw std::runtime_error("Invalid line format: " + line);

        int source = std::stoi(trim(line.substr(0, arrowPos)));
        std::string targetsStr = line.substr(arrowPos + 4);

        std::istringstream ss(targetsStr);
        std::string token;
        while (std::getline(ss, token, ',')) {
            data.graph[source].push_back(std::stoi(trim(token)));
        }
    }

    return data;
}

void validateInput(const InputData& data) {
    if (data.graph.empty())
        throw std::invalid_argument("Graph must not be empty");
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty())
        return options.outputPath;
    if (!options.inputPath.empty())
        return options.inputPath + ".out";
    return "maximal_nonbranching_paths_output.txt";
}

std::string formatResult(const std::vector<std::vector<int>>& paths) {
    std::string result;
    for (const auto& path : paths) {
        if (!result.empty()) result += '\n';
        for (size_t i = 0; i < path.size(); ++i) {
            if (i > 0) result += " -> ";
            result += std::to_string(path[i]);
        }
    }
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
            throw std::invalid_argument("Positional arguments not supported; use --input <file>");
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);

        auto paths = maximalNonBranchingPaths(input.graph);
        std::string result = formatResult(paths);

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3m-maximal_nonbranching_paths [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3m-maximal_nonbranching_paths input.txt\n";
        return 1;
    }
}
