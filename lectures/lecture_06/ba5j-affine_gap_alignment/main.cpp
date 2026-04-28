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

static int aaIndex(char c) {
    auto pos = AA_ALPHA.find(c);
    if (pos == std::string::npos)
        throw std::runtime_error(std::string("Unknown amino acid: ") + c);
    return static_cast<int>(pos);
}

static int blosum(char a, char b) {
    return BLOSUM62[aaIndex(a)][aaIndex(b)];
}

static const int NEG_INF = std::numeric_limits<int>::min() / 2;

struct Alignment {
    int score;
    std::string alignedV;
    std::string alignedW;
};

// State constants for backtracking
static const int8_t FROM_LOWER  = 0; // gap in v continues/opens
static const int8_t FROM_UPPER  = 1; // gap in w continues/opens
static const int8_t FROM_DIAG   = 2; // match/mismatch
static const int8_t OPEN_LOWER  = 3; // opened gap in v (from middle)
static const int8_t OPEN_UPPER  = 4; // opened gap in w (from middle)

Alignment affineGapAlignment(const std::string& v, const std::string& w,
                              int sigma, int epsilon) {
    int n = static_cast<int>(v.size());
    int m = static_cast<int>(w.size());

    // Three DP matrices
    std::vector<std::vector<int>> lower(n+1, std::vector<int>(m+1, NEG_INF));
    std::vector<std::vector<int>> middle(n+1, std::vector<int>(m+1, NEG_INF));
    std::vector<std::vector<int>> upper(n+1, std::vector<int>(m+1, NEG_INF));

    // Backtrack arrays: for each matrix, store what transition was taken
    std::vector<std::vector<int8_t>> backL(n+1, std::vector<int8_t>(m+1, -1));
    std::vector<std::vector<int8_t>> backU(n+1, std::vector<int8_t>(m+1, -1));
    std::vector<std::vector<int8_t>> backM(n+1, std::vector<int8_t>(m+1, -1));

    middle[0][0] = 0;

    // Init first row: gap in v (horizontal)
    for (int j = 1; j <= m; ++j) {
        lower[0][j]  = -sigma - (j-1) * epsilon;
        middle[0][j] = lower[0][j];
        backL[0][j]  = (j == 1) ? OPEN_LOWER : FROM_LOWER;
        backM[0][j]  = FROM_LOWER;
    }
    // Init first column: gap in w (vertical)
    for (int i = 1; i <= n; ++i) {
        upper[i][0]  = -sigma - (i-1) * epsilon;
        middle[i][0] = upper[i][0];
        backU[i][0]  = (i == 1) ? OPEN_UPPER : FROM_UPPER;
        backM[i][0]  = FROM_UPPER;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            // lower: gap in v (came from left)
            int lFromL = (lower[i][j-1]  != NEG_INF) ? lower[i][j-1]  - epsilon : NEG_INF;
            int lFromM = (middle[i][j-1] != NEG_INF) ? middle[i][j-1] - sigma   : NEG_INF;
            if (lFromL >= lFromM) { lower[i][j] = lFromL; backL[i][j] = FROM_LOWER; }
            else                  { lower[i][j] = lFromM; backL[i][j] = OPEN_LOWER; }

            // upper: gap in w (came from above)
            int uFromU = (upper[i-1][j]  != NEG_INF) ? upper[i-1][j]  - epsilon : NEG_INF;
            int uFromM = (middle[i-1][j] != NEG_INF) ? middle[i-1][j] - sigma   : NEG_INF;
            if (uFromU >= uFromM) { upper[i][j] = uFromU; backU[i][j] = FROM_UPPER; }
            else                  { upper[i][j] = uFromM; backU[i][j] = OPEN_UPPER; }

            // middle
            int mFromL = lower[i][j];
            int mFromU = upper[i][j];
            int mFromD = (middle[i-1][j-1] != NEG_INF)
                           ? middle[i-1][j-1] + blosum(v[i-1], w[j-1])
                           : NEG_INF;
            middle[i][j] = std::max({mFromL, mFromU, mFromD});
            if      (mFromD >= mFromL && mFromD >= mFromU) backM[i][j] = FROM_DIAG;
            else if (mFromL >= mFromU)                     backM[i][j] = FROM_LOWER;
            else                                           backM[i][j] = FROM_UPPER;
        }
    }

    // Find best final state
    int bestScore;
    int8_t startState;
    if (middle[n][m] >= lower[n][m] && middle[n][m] >= upper[n][m]) {
        bestScore = middle[n][m]; startState = FROM_DIAG;
    } else if (lower[n][m] >= upper[n][m]) {
        bestScore = lower[n][m]; startState = FROM_LOWER;
    } else {
        bestScore = upper[n][m]; startState = FROM_UPPER;
    }

    // Backtrack
    std::string alignV, alignW;
    int i = n, j = m;
    int8_t state = startState;
    // Recode: 0=in lower, 1=in upper, 2=in middle
    // When start is FROM_DIAG, we're in middle
    // When start is FROM_LOWER, we're in lower
    // When start is FROM_UPPER, we're in upper
    int curState = (startState == FROM_LOWER) ? 0 : (startState == FROM_UPPER) ? 1 : 2;

    while (i > 0 || j > 0) {
        if (curState == 0) { // in lower matrix (gap in v)
            alignV += '-';
            alignW += w[j-1];
            --j;
            // Did we extend or open?
            int8_t bt = backL[i][j+1]; // we came to (i, j+1)
            curState = (bt == FROM_LOWER) ? 0 : 2; // extend stays in lower, open goes to middle
        } else if (curState == 1) { // in upper matrix (gap in w)
            alignV += v[i-1];
            alignW += '-';
            --i;
            int8_t bt = backU[i+1][j]; // we came to (i+1, j)
            curState = (bt == FROM_UPPER) ? 1 : 2;
        } else { // in middle matrix
            int8_t bt = backM[i][j];
            if (bt == FROM_DIAG) {
                alignV += v[i-1];
                alignW += w[j-1];
                --i; --j;
                curState = 2;
            } else if (bt == FROM_LOWER) {
                curState = 0;
                // Don't move i or j yet — process in lower next iteration
            } else { // FROM_UPPER
                curState = 1;
            }
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
    return {"PRTEINS", "PRTWPSEIN"};
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
    return "ba5j-affine_gap_alignment_output.txt";
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
        Alignment aln = affineGapAlignment(input.v, input.w, 11, 1);

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(aln, outputPath);

        std::cout << aln.score << '\n' << aln.alignedV << '\n' << aln.alignedW << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "Usage: ./ba5j-affine_gap_alignment [--input in.txt] [--output out.txt]\n";
        return 1;
    }
}
