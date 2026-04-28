#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// BLOSUM62, amino acids in alphabetical order: ACDEFGHIKLMNPQRSTVWY
static const int BLOSUM62[20][20] = {
    //A   C   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   T   V   W   Y
    { 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2}, // A
    { 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2}, // C
    {-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3}, // D
    {-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2}, // E
    {-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3}, // F
    { 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3}, // G
    {-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2}, // H
    {-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1}, // I
    {-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2}, // K
    {-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1}, // L
    {-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1}, // M
    {-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2}, // N
    {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3}, // P
    {-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1}, // Q
    {-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2}, // R
    { 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2}, // S
    { 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2}, // T
    { 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1}, // V
    {-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2}, // W
    {-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7}, // Y
};

static const std::string AA_ALPHA = "ACDEFGHIKLMNPQRSTVWY";
static const int NEG_INF = std::numeric_limits<int>::min() / 2;

static int aaIndex(char c) {
    auto pos = AA_ALPHA.find(c);
    if (pos == std::string::npos)
        throw std::runtime_error(std::string("Unknown amino acid: ") + c);
    return static_cast<int>(pos);
}

static int blosum(char a, char b) {
    return BLOSUM62[aaIndex(a)][aaIndex(b)];
}

// Global alignment score column at position j (linear gap penalty)
// Returns vector of length n+1 with scores dp[0..n][j]
static std::vector<int> alignColumn(const std::string& v, const std::string& w,
                                    int sigma, int targetCol) {
    int n = static_cast<int>(v.size());
    std::vector<int> prev(n + 1), curr(n + 1);
    for (int i = 0; i <= n; ++i) prev[i] = -i * sigma;

    for (int j = 1; j <= targetCol; ++j) {
        curr[0] = -j * sigma;
        for (int i = 1; i <= n; ++i) {
            int diag = prev[i-1] + blosum(v[i-1], w[j-1]);
            int up   = prev[i]   - sigma;
            int left = curr[i-1] - sigma;
            curr[i] = std::max({diag, up, left});
        }
        std::swap(prev, curr);
    }
    return prev; // after swap, prev holds the result at targetCol
}

struct MiddleEdge {
    int r1, c1, r2, c2; // edge from (r1,c1) to (r2,c2)
};

MiddleEdge findMiddleEdge(const std::string& v, const std::string& w, int sigma) {
    int n = static_cast<int>(v.size());
    int m = static_cast<int>(w.size());
    int mid = m / 2;

    // from_source[i] = best score to align v[0..i-1] with w[0..mid-1]
    std::vector<int> fromSrc = alignColumn(v, w, sigma, mid);

    // to_sink[i] = best score to align v[i..n-1] with w[mid+1..m-1]
    // Computed by reversing: v_rev and w_suf_rev = reverse(w[mid+1..m-1])
    std::string vRev(v.rbegin(), v.rend());
    std::string wSufRev(w.rbegin(), w.rbegin() + (m - mid - 1));

    std::vector<int> revCol = alignColumn(vRev, wSufRev, sigma,
                                          static_cast<int>(wSufRev.size()));
    // to_sink[i] = revCol[n - i]
    auto toSink = [&](int i) -> int { return revCol[n - i]; };

    // Find best middle edge: (i, mid) → (i, mid+1) or (i, mid) → (i+1, mid+1)
    int bestScore = NEG_INF;
    int bestI = 0;
    bool bestDiag = false;

    for (int i = 0; i <= n; ++i) {
        // Horizontal edge: (i, mid) → (i, mid+1), score -sigma
        if (mid + 1 <= m) {
            int total = fromSrc[i] - sigma + toSink(i);
            if (total > bestScore) {
                bestScore = total;
                bestI = i;
                bestDiag = false;
            }
        }
        // Diagonal edge: (i, mid) → (i+1, mid+1)
        if (i < n && mid < m) {
            int total = fromSrc[i] + blosum(v[i], w[mid]) + toSink(i + 1);
            if (total > bestScore) {
                bestScore = total;
                bestI = i;
                bestDiag = true;
            }
        }
    }

    if (bestDiag)
        return {bestI, mid, bestI + 1, mid + 1};
    else
        return {bestI, mid, bestI, mid + 1};
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
    return {"PLEASANTLY", "MEANLY"};
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
    return "ba5k-middle_edge_output.txt";
}

void writeOutputFile(const MiddleEdge& e, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output: " + path);
    f << '(' << e.r1 << ", " << e.c1 << ") (" << e.r2 << ", " << e.c2 << ")\n";
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
        MiddleEdge e = findMiddleEdge(input.v, input.w, 5);

        std::string result = "(" + std::to_string(e.r1) + ", " + std::to_string(e.c1) +
                             ") (" + std::to_string(e.r2) + ", " + std::to_string(e.c2) + ")";

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(e, outputPath);

        std::cout << result << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "Usage: ./ba5k-middle_edge [--input in.txt] [--output out.txt]\n";
        return 1;
    }
}
