#include "RasterData.cpp"
#include "Utility.h"
#include "fmt-9.1.0/include/fmt/format.h"
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
                utils::printLine(fmt::format("{} - x: {}, lon: {}, lonToXResult: {}", kFunctionName,
                                             x, lon, lonToXResult));
                return false;
            }
        }
        return true;
    }

    bool inversesTestYLat() {
        static const std::string kFunctionName{"inversesTestYLat"};
        for (int y = 0; y < kNumRows; y += 600) {
            const double lat = rd.lat(y);
            const int latToYResult = rd.latToY(lat);
            if (latToYResult != y) {
                utils::printLine(fmt::format("{} - y: {}, lat: {}, latToYResult: {}", kFunctionName,
                                             y, lat, latToYResult));
                return false;
            }
        }
        return true;
    }

  private:
    RasterData rd;
};

int main() {
    const std::string sumTableFilename = "popSumTable2020.bin";
    RasterDataTests rdTests(sumTableFilename);
    rdTests.inversesTestXLon();
    rdTests.inversesTestYLat();
}