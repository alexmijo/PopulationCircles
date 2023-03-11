#include "RasterData.cpp"
#include "Utility.h"
#include "fmt/include/fmt/core.h"
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

// constexpr int k30ArcSecsPerDegree = 2 * 60;
// constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
// constexpr int kNumRows = 180 * k30ArcSecsPerDegree;
constexpr bool kUse2020Data = true;
constexpr double kWorldPop = kUse2020Data ? 7757982599.3135586 : 7346242908.863955;

class RasterDataTests {
  public:
    RasterDataTests(const std::string &sumTableFilename) { loadData(sumTableFilename); }
    RasterDataTests() {}

    void loadData(const std::string &sumTableFilename) { rd.loadData(sumTableFilename); }

    bool inversesTestXLon() {
        static const std::string kFunctionName{"inversesTestXLon"};
        for (int x = 0; x < kNumCols; x += 600) {
            const double lon = rd.lon(x);
            const int lonToXResult = rd.lonToX(lon);
            if (lonToXResult != x) {
                utils::printLine(kFunctionName + " - x: " + std::to_string(x) +
                                 ", lon: " + std::to_string(lon) +
                                 ", lonToXResult: " + std::to_string(lonToXResult));
                // utils::printLine(fmt::format("{} - x: {}, lon: {}, lonToXResult: {}",
                // kFunctionName, x, lon, lonToXResult));
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
                utils::printLine(kFunctionName + " - y: " + std::to_string(y) +
                                 ", lat: " + std::to_string(lat) +
                                 ", latToYResult: " + std::to_string(latToYResult));
                // utils::printLine(fmt::format("{} - y: {}, lat: {}, latToYResult: {}",
                // kFunctionName, y, lat, latToYResult));
                return false;
            }
        }
        return true;
    }

    bool totalPopulationRectangle() {
        static const std::string kFunctionName{"totalPopulationRectangle"};
        if (!rd.hasData()) {
            return false;
        }
        PixelRectangle pr{0, kNumCols, 0, kNumRows};
        const double result = rd.sumWithinRectangle(pr);
        if (result != kWorldPop) {
            utils::printLine(kFunctionName + " - result: " + std::to_string(result) +
                             ", kWorldPop: " + std::to_string(kWorldPop));
            return false;
        }
        return true;
    }

    int numTestsRun() { return m_NumTestsRun; }
    int numTestsFailed() { return m_NumTestsFailed; }

    bool runFastTests(bool resetNumTests = true) {
        if (resetNumTests) {
            m_NumTestsRun = 0;
            m_NumTestsFailed = 0;
        }
        const int startingNumTestsFailed = m_NumTestsFailed;

        m_NumTestsRun++;
        m_NumTestsFailed += inversesTestXLon() ? 0 : 1;

        m_NumTestsRun++;
        m_NumTestsFailed += inversesTestYLat() ? 0 : 1;

        return m_NumTestsFailed == startingNumTestsFailed;
    }

    bool runSlowTests(bool resetNumTests = true) {
        if (resetNumTests) {
            m_NumTestsRun = 0;
            m_NumTestsFailed = 0;
        }
        int startingNumTestsFailed = m_NumTestsFailed;
        if (!rd.hasData()) {
            return false;
        }

        m_NumTestsRun++;
        m_NumTestsFailed += totalPopulationRectangle() ? 0 : 1;

        return m_NumTestsFailed == startingNumTestsFailed;
    }

  private:
    RasterData rd;
    int m_NumTestsRun{0};
    int m_NumTestsFailed{0};
};

int main() {
    const std::string sumTableFilename = "popSumTable2020.bin";
    // RasterDataTests rdTests(sumTableFilename);
    // if (rdTests.runTests()) {
    //     utils::printLine("Passed " + std::to_string(rdTests.numTestsRun()) + " tests");
    //     // utils::printLine(fmt::format("Passed {} tests", rdTests.numTestsRun()));
    // } else {
    //     utils::printLine("Failed " + std::to_string(rdTests.numTestsFailed()) + "/" +
    //                      std::to_string(rdTests.numTestsRun()) + " tests");
    //     // utils::printLine(
    //     //     fmt::format("Failed {}/{} tests", rdTests.numTestsFailed() /
    //     rdTests.numTestsRun()));
    // }

    RasterDataTests rdTests;
    if (rdTests.runFastTests()) {
        utils::printLine("Passed " + std::to_string(rdTests.numTestsRun()) + " tests");
    } else {
        utils::printLine("Failed " + std::to_string(rdTests.numTestsFailed()) + "/" +
                         std::to_string(rdTests.numTestsRun()) + " tests");
    }
    
    rdTests.loadData(sumTableFilename);
    if (rdTests.runSlowTests()) {
        utils::printLine("Passed " + std::to_string(rdTests.numTestsRun()) + " tests");
    } else {
        utils::printLine("Failed " + std::to_string(rdTests.numTestsFailed()) + "/" +
                         std::to_string(rdTests.numTestsRun()) + " tests");
    }
}