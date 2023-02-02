#include "RasterData.cpp"
#include <iostream>

constexpr int k30ArcSecsPerDegree = 2 * 60;
constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
constexpr int kNumRows = 180 * k30ArcSecsPerDegree;

int main() {
    std::string sumTableFilename = "popSumTable2020.bin";
    RasterData rd(sumTableFilename);

    for (int x = 0; x < kNumCols; x += 600) {
        std::cout << "x: " << x << "\n";
        const auto lon = rd.lon(x);
        std::cout << "lon: " << lon << "\n";
        const auto lonToXResult = rd.lonToX(lon);
        std::cout << "lonToXResult: " << lonToXResult << "\n";
    }

    for (int y = 0; y < kNumRows; y += 600) {
        std::cout << "y: " << y << "\n";
        const auto lat = rd.lat(y);
        std::cout << "lat: " << lat << "\n";
        const auto latToYResult = rd.latToY(lat);
        std::cout << "latToYResult: " << latToYResult << "\n";
    }
}