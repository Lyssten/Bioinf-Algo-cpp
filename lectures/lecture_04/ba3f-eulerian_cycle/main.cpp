#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
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

std::vector<int> eulerianCycle(std::map<int, std::vector<int>>& graph) {
    std::map<int, size_t> edgeIndex;
    for (auto& [node, neighbors] : graph) {
        edgeIndex[node] = 0;
    }

    int startNode = graph.begin()->first;
    std::vector<int> stack;
    std::vector<int> cycle;

    stack.push_back(startNode);

    while (!stack.empty()) {
        int v = stack.back();
        if (edgeIndex[v] < graph[v].size()) {
            int u = graph[v][edgeIndex[v]];
            edgeIndex[v]++;
            stack.push_back(u);
        } else {
            cycle.push_back(v);
            stack.pop_back();
        }
    }

    std::reverse(cycle.begin(), cycle.end());
    return cycle;
}

InputData sampleInput() {
    InputData data;
    data.graph[0] = {3};
    data.graph[1] = {0};
    data.graph[2] = {1, 6};
    data.graph[3] = {2};
    data.graph[4] = {2};
    data.graph[5] = {4};
    data.graph[6] = {5, 8};
    data.graph[7] = {9};
    data.graph[8] = {7};
    data.graph[9] = {6};
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
            throw std::runtime_error("Invalid line format (expected ' -> '): " + line);

        int source = std::stoi(trim(line.substr(0, arrowPos)));
        std::string targetsStr = line.substr(arrowPos + 4);

        std::vector<int> targets;
        std::istringstream ss(targetsStr);
        std::string token;
        while (std::getline(ss, token, ',')) {
            targets.push_back(std::stoi(trim(token)));
        }

        data.graph[source] = targets;
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
    return "eulerian_cycle_output.txt";
}

std::string formatResult(const std::vector<int>& cycle) {
    std::string result;
    for (size_t i = 0; i < cycle.size(); ++i) {
        if (i > 0) result += "->";
        result += std::to_string(cycle[i]);
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

        std::vector<int> cycle = eulerianCycle(input.graph);
        std::string result = formatResult(cycle);

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(result, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        std::cerr << "Usage:\n";
        std::cerr << "  ./ba3f-eulerian_cycle [--input input.txt] [--output output.txt]\n";
        std::cerr << "  ./ba3f-eulerian_cycle input.txt [--output output.txt]\n";
        return 1;
    }
}
