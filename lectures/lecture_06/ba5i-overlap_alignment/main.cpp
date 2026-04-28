#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// match=+1, mismatch=-2, indel=-2
static const int MATCH = 1;
static const int MIS   = -2;
static const int SIGMA = 2;

static int score(char a, char b) {
    return (a == b) ? MATCH : MIS;
}

struct Alignment {
    int score;
    std::string alignedV;
    std::string alignedW;
};

Alignment overlapAlignment(const std::string& v, const std::string& w) {
    int n = static_cast<int>(v.size());
    int m = static_cast<int>(w.size());

    // dp[i][j] = best score aligning v[0..i-1] with w[0..j-1]
    // dp[i][0] = 0 for all i: free to skip any prefix of v
    // dp[0][j] = -j*SIGMA: must pay for leading gaps in v for w's prefix
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    for (int j = 1; j <= m; ++j)
        dp[0][j] = -j * SIGMA;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int diag = dp[i-1][j-1] + score(v[i-1], w[j-1]);
            int up   = dp[i-1][j]   - SIGMA;
            int left = dp[i][j-1]   - SIGMA;
            dp[i][j] = std::max({diag, up, left});
        }
    }

    // Find last j with maximum score in last row (free to skip suffix of w)
    int bestScore = dp[n][0];
    int bestJ = 0;
    for (int j = 1; j <= m; ++j) {
        if (dp[n][j] >= bestScore) {
            bestScore = dp[n][j];
            bestJ = j;
        }
    }

    // Backtrack from (n, bestJ)
    std::string alignV, alignW;
    int i = n, j = bestJ;
    while (i > 0 && j > 0) {
        if (dp[i][j] == dp[i-1][j-1] + score(v[i-1], w[j-1])) {
            alignV += v[i-1]; alignW += w[j-1]; --i; --j;
        } else if (dp[i][j] == dp[i-1][j] - SIGMA) {
            alignV += v[i-1]; alignW += '-'; --i;
        } else {
            alignV += '-'; alignW += w[j-1]; --j;
        }
    }
    // Remaining w consumed as leading gaps in v (from dp[0][j] = -j*SIGMA)
    while (j > 0) { alignV += '-'; alignW += w[j-1]; --j; }
    // Remaining i: free prefix of v — NOT added to output

    std::reverse(alignV.begin(), alignV.end());
    std::reverse(alignW.begin(), alignW.end());

    return {bestScore, alignV, alignW};
}

struct InputData {
    std::string v;
    std::string w;
};

struct CliOptions {
    std::string inputPath;
    std::string outputPath;
    std::vector<std::string> positional;
};

InputData sampleInput() {
    return {"PAWHEAE", "HEAGAWGHEE"};
}

bool isReadableFile(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

CliOptions parseCli(int argc, char* argv[]) {
    CliOptions options;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i + 1 >= argc) throw std::invalid_argument("Missing value after --input");
            options.inputPath = argv[++i];
        } else if (arg == "--output" || arg == "-o") {
            if (i + 1 >= argc) throw std::invalid_argument("Missing value after --output");
            options.outputPath = argv[++i];
        } else {
            options.positional.push_back(arg);
        }
    }
    if (options.inputPath.empty() && options.positional.size() == 1 &&
        isReadableFile(options.positional[0])) {
        options.inputPath = options.positional[0];
        options.positional.clear();
    }
    return options;
}

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

InputData parseInputFile(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    InputData data;
    std::getline(f, data.v); data.v = trim(data.v);
    std::getline(f, data.w); data.w = trim(data.w);
    if (data.v.empty() || data.w.empty())
        throw std::runtime_error("Expected two non-empty lines");
    return data;
}

void validateInput(const InputData& data) {
    if (data.v.empty() || data.w.empty())
        throw std::invalid_argument("Both strings must be non-empty");
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty()) return options.outputPath;
    if (!options.inputPath.empty())  return options.inputPath + ".out";
    return "ba5i-overlap_alignment_output.txt";
}

void writeOutputFile(const Alignment& aln, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output: " + path);
    f << aln.score << '\n' << aln.alignedV << '\n' << aln.alignedW << '\n';
}

int main(int argc, char* argv[]) {
    try {
        CliOptions options = parseCli(argc, argv);
        InputData input = sampleInput();

        if (!options.inputPath.empty()) {
            input = parseInputFile(options.inputPath);
        } else if (!options.positional.empty()) {
            throw std::invalid_argument("Positional args not supported; use --input <file>");
        } else {
            std::cerr << "No input provided. Using sample dataset.\n";
        }

        validateInput(input);
        Alignment aln = overlapAlignment(input.v, input.w);

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(aln, outputPath);

        std::cout << aln.score << '\n' << aln.alignedV << '\n' << aln.alignedW << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "Usage: ./ba5i-overlap_alignment [--input in.txt] [--output out.txt]\n";
        return 1;
    }
}
