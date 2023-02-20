#include "RasterData.cpp"
#include <chrono>
#include <iostream>

constexpr int k30ArcSecsPerDegree = 2 * 60;
constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
constexpr int kNumRows = 180 * k30ArcSecsPerDegree;

class RasterDataTests {
  public:
    RasterDataTests(const std::string &sumTableFilename) : rd(sumTableFilename) {}

    bool inversesTestXLon() {
        static const std::string kFunctionName{"inversesTestXLon"};
        for (int x = 0; x < kNumCols; x += 600) {
            const double lon = rd.lon(x);
            const int lonToXResult = rd.lonToX(lon);
            if (lonToXResult != x) {
                std::cout << "x: " << x << "\n";
                std::cout << "lon: " << lon << "\n";
                std::cout << "lonToXResult: " << lonToXResult << "\n";
            }
        }
    }

    bool inversesTestYLat() {
        for (int y = 0; y < kNumRows; y += 600) {
            std::cout << "y: " << y << "\n";
            const auto lat = rd.lat(y);
            std::cout << "lat: " << lat << "\n";
            const auto latToYResult = rd.latToY(lat);
            std::cout << "latToYResult: " << latToYResult << "\n";
        }
    }

  private:
    RasterData rd;
};

int main() {
    std::string sumTableFilename = "popSumTable2020.bin";
    RasterData rd(sumTableFilename);
}