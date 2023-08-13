#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <vector>

//--------------------------------------------------------------------------------------------------
constexpr int k30ArcSecsPerDegree = 2 * 60;
constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
constexpr int kNumRows = 180 * k30ArcSecsPerDegree;

//--------------------------------------------------------------------------------------------------
void print(const std::string &line) { std::cout << line; }

void printLine(const std::string &line) {
    print(line);
    std::cout << "\n";
}

//--------------------------------------------------------------------------------------------------
struct Location {
    double lat;
    double lon;
};

// 0 indexed, x from west to east, y from south to north
struct Pixel {
    int x;
    int y;
};

struct LatLonRect {
    // TODO
};

// (x ∈ [W, E], y ∈ [S, N])
struct PixelRect {
    int W;
    int E;
    int S;
    int N;

    // TODO: Test timings with and without & (both here and everywhere else. How small is enough?)
    bool contains(const Pixel &p) const { return W <= p.x && p.x <= E && N <= p.y && p.y <= S; }
};

// struct OffsetPixelRect {
//     int offsetW;
//     int offsetE;
//     int offsetS;
//     int offsetN;

//     // TODO: dims aren't a Pixel
//     std::vector<PixelRect> ToPixelRect(const Pixel &center, const Pixel &gridSize) {
//         std::vector<PixelRect> rects;

//         int W = center.x + offsetW;
//         int E = center.x + offsetE;
//         int S = center.y + offsetS;
//         int N = center.y + offsetN;

//         if (W < 0) {
//             rects.emplace_back(0, E, S, N);
//             rects.emplace_back(gridSize.x + W, gridSize.x, S, N);
//         } else if (E >= gridSize.x) {
//             rects.emplace_back(W, gridSize.x - 1, S, N);
//             rects.emplace_back(0, E - gridSize.x, S, N);
//         } else {
//             rects.emplace_back(W, E, S, N);
//         }

//         return rects;
//     }
// };

// Summary: Setting up Karol's BCBE for trading, tradenotificationlistener subscribe by portfolio,
// was working on stuff related to moving the mdrc configs,
// 1. BC UI to OMS staging
// 3. Futures order python client to BCBE to staging reflector  verify boltx looks as expected

// 2. Python client to OMS staging

struct PixelCircle {
    Pixel center;
    // Each rect stretches from west to east side of circle
    std::vector<PixelRect> rects;
};

//--------------------------------------------------------------------------------------------------
// Returns distance in kilometers between two points on Earth's surface.
// TODO: Find source and give credit
static double distance(Location loc1, Location loc2) {
    double req = 6378137.0;

    // alexmijo added code for dealing with equatorial points or identical points
    if (loc1.lat == loc2.lat && loc1.lon == loc2.lon) {
        return 0.0;
    } else if (loc1.lat == 0.0 && loc2.lat == 0.0) {
        double lonDiff = abs(loc1.lon - loc2.lon);
        double fractionOfEquator = lonDiff / 360.0;
        double equatorLength = 2.0 * 3.14159265358979323 * req;
        return fractionOfEquator * equatorLength;
    }

    const double latitude_01 = loc1.lat * M_PI / 180.0;
    const double longitude_01 = loc1.lon * M_PI / 180.0;

    const double latitude_02 = loc2.lat * M_PI / 180.0;
    const double longitude_02 = loc2.lon * M_PI / 180.0;

    const double a = 6378137.0;
    const double b = 6356752.314245;

    // Flattening
    const double f = (a - b) / a;

    // tan U1
    const double tan_U1 = (1 - f) * std::tan(latitude_01);
    const double tan_U2 = (1 - f) * std::tan(latitude_02);

    // Longitudinal Distance
    const double cos_U1 = 1 / std::sqrt(1 + tan_U1 * tan_U1);
    const double cos_U2 = 1 / std::sqrt(1 + tan_U2 * tan_U2);
    const double sin_U1 = tan_U1 * cos_U1;
    const double sin_U2 = tan_U2 * cos_U2;

    // Iterate until complete
    const double L = longitude_02 - longitude_01;
    double lambda = L;
    double diff, sigma;
    double cos_alpha_sq, cos_2sigma_m;
    double cos_sigma, sin_sigma;

    while (true) {

        double sin_lambda = std::sin(lambda);
        double cos_lambda = std::cos(lambda);

        double c1 = (cos_U2 * sin_lambda) * (cos_U2 * sin_lambda);
        double c2 = (cos_U1 * sin_U2);
        double c3 = (sin_U1 * cos_U2 * cos_lambda);

        //  sin sigma
        sin_sigma = std::sqrt(c1 + (c2 - c3) * (c2 - c3));

        // cos sigma
        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda;

        // sigma
        sigma = std::atan2(sin_sigma, cos_sigma);

        // sin alpha
        double sin_alpha = (cos_U1 * cos_U2 * sin_lambda) / (sin_sigma);

        // cos^2 alpha
        cos_alpha_sq = 1 - (sin_alpha * sin_alpha);

        // cos^2 2sigmam
        cos_2sigma_m = cos_sigma - (2 * sin_U1 * sin_U2) / (cos_alpha_sq);

        // C
        double C = (f / 16.0) * cos_alpha_sq * (4 + f * (4 - 3 * cos_alpha_sq));

        // Update Lambda
        diff = lambda;
        lambda = L + (1 - C) * f * sin_alpha *
                         (sigma + C * sin_sigma *
                                      (cos_2sigma_m +
                                       C * cos_sigma * (-1 + 2 * cos_2sigma_m * cos_2sigma_m)));
        diff = lambda - diff;
        if (std::fabs(diff) < 0.00001) {
            break;
        }
    }

    // U2
    double u_sq = cos_alpha_sq * (a * a - b * b) / (b * b);

    // Compute A, B
    double A = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));

    double B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))));

    // Sigma
    double cos_2sigma_m_sq = cos_2sigma_m * cos_2sigma_m;
    double delta_sigma =
        B * sin_sigma *
        (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq) -
                                     (B / 6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma * sin_sigma) *
                                         (-3 + 4 * cos_2sigma_m_sq)));

    // Distance
    double s = b * A * (sigma - delta_sigma);

    // Convert from meters to kilometers before returning
    return s / 1000.;
}

//-----------------------------------------------------------------------------SumTable.cpp
template <typename T> class SumTable {
  public:
    SumTable(const std::string &sumTableFilename) {
        std::ifstream sumTableFile(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&height), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&width), sizeof(int));
        // The 0s will remain only as a column of 0s to the west and a row of 0s to the south.
        // This is padding for the data which will replace the rest of the zeros, and explains why
        // sumTable has dimensions one greater than both height and width.
        sumTable.resize(height + 1, std::vector<T>(width + 1, 0));
        for (int r = 1; r < height + 1; r++) {
            // Read in an entire row at once
            sumTableFile.read(reinterpret_cast<char *>(&sumTable[r][1]), sizeof(T) * width);
        }
    }

    // TODO: Something with the padding?
    // TODO: Move?
    SumTable(const std::vector<std::vector<T>> &sumTableOrValues, const bool isSumTable) {
        if (isSumTable) {
            sumTable = sumTableOrValues;
            height = sumTableOrValues.size() - 1;
            width = sumTableOrValues[0].size();
        } else {
            height = sumTableOrValues.size();
            width = height > 0 ? sumTableOrValues[0].size() : 0;

            // Initialize the sumTable with an extra row and column for padding
            sumTable.resize(height + 1, std::vector<T>(width + 1, 0));

            // Calculate the sumTable
            for (int r = 1; r <= height; ++r) {
                for (int c = 1; c <= width; ++c) {
                    sumTable[r][c] = sumTableOrValues[r - 1][c - 1] + sumTable[r - 1][c] +
                                     sumTable[r][c - 1] - sumTable[r - 1][c - 1];
                }
            }
        }
    }

    T sumWithinRectangle(const PixelRect &rect) {
        return sumTable[rect.N + 1][rect.E + 1] - sumTable[rect.S][rect.E + 1] -
               sumTable[rect.N + 1][rect.W] + sumTable[rect.S][rect.W];
    }

    T sumWithinRectangle(const int W, const int E, const int S, const int N) {
        return sumWithinRectangle({W, E, S, N});
    }

    int width;
    int height;

  private:
    std::vector<std::vector<T>> sumTable;
};

//-----------------------------------------------------------------------------------------Tests.cpp
class TestRunner {
  public:
    void addTest(const std::function<bool()> &testFunction, const std::string &testName = "") {
        tests.emplace_back(testName, testFunction);
    }

    void runTests() {
        for (const auto &test : tests) {
            bool result = test.second();
            std::cout << "Test '" << test.first << "' ";
            if (result) {
                std::cout << "Passed\n";
            } else {
                std::cout << "Failed\n";
            }
        }
    }

  private:
    std::vector<std::pair<std::string, std::function<bool()>>> tests;
};

//-----------------------------------------------------------------------------SumTableTests.cpp
bool testInitializeSumTableFromValues() {
    // 1 2 3 4
    // 0 1 2 3
    std::vector<std::vector<double>> values = {{0, 1, 2, 3}, {1, 2, 3, 4}};
    SumTable<double> sumTable(values, false);

    // 1 4 9 16
    // 0 1 3 6
    // Manually computed expected sums
    std::vector<std::vector<double>> expectedSums = {{0, 1, 3, 6}, {1, 4, 9, 16}};

    for (int r = 0; r < values.size(); ++r) {
        for (int c = 0; c < values[0].size(); ++c) {
            if (sumTable.sumWithinRectangle({0, c, 0, r}) != expectedSums[r][c]) {
                std::cout << "Actual: " << sumTable.sumWithinRectangle({0, c, 0, r})
                          << ", Expected: " << expectedSums[r][c] << std::endl;
                return false;
            }

            if (sumTable.sumWithinRectangle({c, c, r, r}) != values[r][c]) {
                std::cout << "Actual: " << sumTable.sumWithinRectangle({c, c, r, r})
                          << ", Expected: " << values[r][c] << std::endl;
                return false;
            }

            if (sumTable.sumWithinRectangle(0, c, 0, r) != expectedSums[r][c]) {
                std::cout << "Actual: " << sumTable.sumWithinRectangle(0, c, 0, r)
                          << ", Expected: " << expectedSums[r][c] << std::endl;
                return false;
            }

            if (sumTable.sumWithinRectangle(c, c, r, r) != values[r][c]) {
                std::cout << "Actual: " << sumTable.sumWithinRectangle(c, c, r, r)
                          << ", Expected: " << values[r][c] << std::endl;
                return false;
            }
        }
    }

    return true;
}

bool testInitializeSumTableFromFile() { return true; }

void main1() {
    TestRunner runner;
    runner.addTest(testInitializeSumTableFromValues, "testInitializeSumTableFromValues");
    runner.addTest(testInitializeSumTableFromFile, "testInitializeSumTableFromFile");
    runner.runTests();
}

//------------------------------------------------------------------------------------ImageMaker.cpp
struct RGB {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

RGB adjustBrightness(const RGB &color, double brightness) {
    RGB newColor;

    newColor.r = std::clamp(static_cast<int>(color.r * brightness), 0, 255);
    newColor.g = std::clamp(static_cast<int>(color.g * brightness), 0, 255);
    newColor.b = std::clamp(static_cast<int>(color.b * brightness), 0, 255);

    return newColor;
}

const std::unordered_map<std::string, const RGB> kColors = {
    {"r", {255, 0, 0}}, {"g", {0, 255, 0}}, {"b", {0, 0, 255}}};

void makePpm(const std::vector<std::vector<std::string>> &pixels) {
    std::cout << "P3\n" << pixels[0].size() << " " << pixels.size() << "\n255\n";
    for (const auto &row : pixels) {
        for (const auto &pixel : row) {
            if (kColors.count(pixel) > 0) {
                const auto &color = kColors.at(pixel);
                std::cout << (int)color.r << " " << (int)color.g << " " << (int)color.b << "\n";
            }
        }
    }
}

void makePpm(const std::vector<std::vector<RGB>> &pixels) {
    std::cout << "P3\n" << pixels[0].size() << " " << pixels.size() << "\n255\n";
    for (const auto &row : pixels) {
        for (const auto &pixel : row) {
            std::cout << (int)pixel.r << " " << (int)pixel.g << " " << (int)pixel.b << "\n";
        }
    }
}

void main2() {
    std::vector<std::vector<std::string>> image1 = {
        {"r", "r", "r", "r", "g", "b"}, {"r", "r", "r", "r", "g", "b"},
        {"r", "r", "r", "r", "g", "b"}, {"r", "r", "r", "r", "g", "b"},
        {"r", "r", "r", "r", "g", "b"}, {"r", "r", "r", "r", "g", "b"}};

    makePpm(image1);
}


//--------------------------------------------------------------------------------------------------
int main() {
    main1();
}