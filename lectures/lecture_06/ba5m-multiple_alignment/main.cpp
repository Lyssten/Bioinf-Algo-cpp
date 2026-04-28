#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Scoring: +1 when all three characters match (no gaps), 0 otherwise (no gap penalty)
struct Alignment {
    int score;
    std::string a1, a2, a3;
};

Alignment multipleAlignment(const std::string& v1, const std::string& v2,
                             const std::string& v3) {
    int n1 = static_cast<int>(v1.size());
    int n2 = static_cast<int>(v2.size());
    int n3 = static_cast<int>(v3.size());

    // dp[i][j][k] = best score aligning v1[0..i-1], v2[0..j-1], v3[0..k-1]
    std::vector<std::vector<std::vector<int>>> dp(
        n1+1, std::vector<std::vector<int>>(n2+1, std::vector<int>(n3+1, 0)));
    // back[i][j][k]: bits d1d2d3 — which sequences advanced into (i,j,k)
    std::vector<std::vector<std::vector<int8_t>>> back(
        n1+1, std::vector<std::vector<int8_t>>(n2+1, std::vector<int8_t>(n3+1, 0)));

    const int NEG_INF = std::numeric_limits<int>::min() / 2;

    for (int i = 0; i <= n1; ++i) {
        for (int j = 0; j <= n2; ++j) {
            for (int k = 0; k <= n3; ++k) {
                if (i == 0 && j == 0 && k == 0) continue;

                int best = NEG_INF;
                int8_t bestDir = 0;

                auto tryUpdate = [&](int pi, int pj, int pk, int8_t dir, int add) {
                    int val = dp[pi][pj][pk] + add;
                    if (val > best) { best = val; bestDir = dir; }
                };

                // All 7 possible predecessor transitions (d1d2d3 encoded as bits)
                if (i>0 && j>0 && k>0) {
                    int s = (v1[i-1]==v2[j-1] && v2[j-1]==v3[k-1]) ? 1 : 0;
                    tryUpdate(i-1, j-1, k-1, 7, s); // 111
                }
                if (i>0 && j>0)         tryUpdate(i-1, j-1, k,   6, 0); // 110
                if (i>0 && k>0)         tryUpdate(i-1, j,   k-1, 5, 0); // 101
                if (j>0 && k>0)         tryUpdate(i,   j-1, k-1, 3, 0); // 011
                if (i>0)                tryUpdate(i-1, j,   k,   4, 0); // 100
                if (j>0)                tryUpdate(i,   j-1, k,   2, 0); // 010
                if (k>0)                tryUpdate(i,   j,   k-1, 1, 0); // 001

                dp[i][j][k] = best;
                back[i][j][k] = bestDir;
            }
        }
    }

    // Backtrack
    std::string a1, a2, a3;
    int i = n1, j = n2, k = n3;
    while (i > 0 || j > 0 || k > 0) {
        int8_t dir = back[i][j][k];
        int d1 = (dir >> 2) & 1;
        int d2 = (dir >> 1) & 1;
        int d3 = (dir >> 0) & 1;
        a1 += (d1 ? v1[i-1] : '-');
        a2 += (d2 ? v2[j-1] : '-');
        a3 += (d3 ? v3[k-1] : '-');
        i -= d1; j -= d2; k -= d3;
    }

    std::reverse(a1.begin(), a1.end());
    std::reverse(a2.begin(), a2.end());
    std::reverse(a3.begin(), a3.end());

    return {dp[n1][n2][n3], a1, a2, a3};
}

struct InputData { std::string v1, v2, v3; };
struct CliOptions {
    std::string inputPath, outputPath;
    std::vector<std::string> positional;
};

InputData sampleInput() {
    return {"ATATCCG", "TCCGA", "ATGTACTG"};
}

bool isReadableFile(const std::string& path) {
    std::ifstream f(path); return f.good();
}

CliOptions parseCli(int argc, char* argv[]) {
    CliOptions options;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i+1 >= argc) throw std::invalid_argument("Missing value after --input");
            options.inputPath = argv[++i];
        } else if (arg == "--output" || arg == "-o") {
            if (i+1 >= argc) throw std::invalid_argument("Missing value after --output");
            options.outputPath = argv[++i];
        } else {
            options.positional.push_back(arg);
        }
    }
    if (options.inputPath.empty() && options.positional.size()==1 &&
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
    std::getline(f, data.v1); data.v1 = trim(data.v1);
    std::getline(f, data.v2); data.v2 = trim(data.v2);
    std::getline(f, data.v3); data.v3 = trim(data.v3);
    return data;
}

void validateInput(const InputData& data) {
    if (data.v1.empty() || data.v2.empty() || data.v3.empty())
        throw std::invalid_argument("All three strings must be non-empty");
}

std::string resolveOutputPath(const CliOptions& options) {
    if (!options.outputPath.empty()) return options.outputPath;
    if (!options.inputPath.empty())  return options.inputPath + ".out";
    return "ba5m-multiple_alignment_output.txt";
}

void writeOutputFile(const Alignment& aln, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output: " + path);
    f << aln.score << '\n' << aln.a1 << '\n' << aln.a2 << '\n' << aln.a3 << '\n';
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
        Alignment aln = multipleAlignment(input.v1, input.v2, input.v3);

        std::string outputPath = resolveOutputPath(options);
        writeOutputFile(aln, outputPath);

        std::cout << aln.score << '\n' << aln.a1 << '\n' << aln.a2 << '\n' << aln.a3 << '\n';
        std::cerr << "Saved output to " << outputPath << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        std::cerr << "Usage: ./ba5m-multiple_alignment [--input in.txt] [--output out.txt]\n";
        return 1;
    }
}
