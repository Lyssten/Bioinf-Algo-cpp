#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// PAM250 substitution matrix, amino acids in alphabetical order: ACDEFGHIKLMNPQRSTVWY
static const int PAM250[20][20] = {
    //A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    { 2, -2,  0,  0, -3,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3}, // A
    {-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0}, // C
    { 0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4}, // D
    { 0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4}, // E
    {-3, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -3, -5, -5, -4, -3, -3, -1,  0,  7}, // F
    { 1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0,  0, -1, -3,  1,  0, -1, -7, -5}, // G
    {-1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0}, // H
    {-1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1}, // I
    {-1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4}, // K
    {-2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1}, // L
    {-1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1, -1, -2, -1,  2, -4, -2}, // M
    { 0, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,  0,  1,  0,  1,  0, -2, -4, -2}, // N
    { 1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,  0,  0,  1,  0, -1, -6, -5}, // P
    { 0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1,  0, -1, -2, -5, -4}, // Q
    {-2, -4, -1, -1, -4, -3,  2, -2,  3, -3, -1,  0,  0,  1,  6,  0, -1, -2,  2, -4}, // R
    { 1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1,  0,  0,  2,  1, -1, -2, -3}, // S
    { 1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3}, // T
    { 0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2}, // V
    {-6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0}, // W
    {-3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10}, // Y
};

static const std::string AA_ALPHA = "ACDEFGHIKLMNPQRSTVWY";

static int aaIndex(char c) {
    auto pos = AA_ALPHA.find(c);
    if (pos == std::string::npos)
        throw std::runtime_error(std::string("Unknown amino acid: ") + c);
    return static_cast<int>(pos);
}

static int score(char a, char b) {
    return PAM250[aaIndex(a)][aaIndex(b)];
}

struct Alignment {
    int score;
    std::string alignedV;
    std::string alignedW;
};

Alignment localAlignment(const std::string& v, const std::string& w, int sigma) {
    int n = static_cast<int>(v.size());
    int m = static_cast<int>(w.size());

    // dp[i][j] = best local alignment score ending at v[i-1], w[j-1]
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

    int bestScore = 0;
    int bestI = 0, bestJ = 0;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int diag = dp[i-1][j-1] + score(v[i-1], w[j-1]);
            int up   = dp[i-1][j]   - sigma;
            int left = dp[i][j-1]   - sigma;
            dp[i][j] = std::max({0, diag, up, left});
            if (dp[i][j] > bestScore) {
                bestScore = dp[i][j];
                bestI = i;
                bestJ = j;
            }
        }
    }

    // Backtrack from (bestI, bestJ) to first 0
    std::string alignV, alignW;
    int i = bestI, j = bestJ;
    while (i > 0 && j > 0 && dp[i][j] != 0) {
        if (dp[i][j] == dp[i-1][j-1] + score(v[i-1], w[j-1])) {
            alignV += v[i-1];
            alignW += w[j-1];
            --i; --j;
        } else if (dp[i][j] == dp[i-1][j] - sigma) {
            alignV += v[i-1];
            alignW += '-';
            --i;
        } else {
            alignV += '-';
            alignW += w[j-1];
            --j;
        }
    }

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
    return {"MEANLY", "PENALTY"};
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
    for (char c : data.v) aaIndex(c);
    for (char c : data.w) aaIndex(c);
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty()) return options.outputPath;
    if (!options.inputPath.empty())  return options.inputPath + ".out";
    return "ba5f-local_alignment_output.txt";
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
        Alignment aln = localAlignment(input.v, input.w, 5);

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(aln, outputPath);

        std::cout << aln.score << '\n' << aln.alignedV << '\n' << aln.alignedW << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "Usage: ./ba5f-local_alignment [--input in.txt] [--output out.txt]\n";
        return 1;
    }
}
