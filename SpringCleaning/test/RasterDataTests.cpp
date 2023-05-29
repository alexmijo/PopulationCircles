#include "../src/RasterData.cpp"
#include "../src/Utility.cpp"

#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

constexpr bool kUse2020Data = true;
// TODO: Constants JSON
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
                printLine(kFunctionName + " - x: " + std::to_string(x) + ", lon: " +
                          std::to_string(lon) + ", lonToXResult: " + std::to_string(lonToXResult));
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
                printLine(kFunctionName + " - y: " + std::to_string(y) + ", lat: " +
                          std::to_string(lat) + ", latToYResult: " + std::to_string(latToYResult));
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
        PixelRectangle pr{0, kNumCols - 1, kNumRows - 1, 0};
        const double result = rd.sumWithinRectangle(pr);
        if (result != kWorldPop) {
            printLine(kFunctionName + " - result: " + std::to_string(result) +
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
    const std::string sumTableFilename =
        kUse2020Data ? "PopData/popSumTable2020.bin" : "PopData/popSumTable2015.bin";

    RasterDataTests rdTests;
    if (rdTests.runFastTests()) {
        printLine("Passed " + std::to_string(rdTests.numTestsRun()) + " tests");
    } else {
        printLine("Failed " + std::to_string(rdTests.numTestsFailed()) + "/" +
                  std::to_string(rdTests.numTestsRun()) + " tests");
    }

    rdTests.loadData(sumTableFilename);
    if (rdTests.runSlowTests()) {
        printLine("Passed " + std::to_string(rdTests.numTestsRun()) + " tests");
    } else {
        printLine("Failed " + std::to_string(rdTests.numTestsFailed()) + "/" +
                  std::to_string(rdTests.numTestsRun()) + " tests");
    }
}