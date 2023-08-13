#include "../src/SumTable.cpp"
#include "Tests.cpp"
#include <vector>

//--------------------------------------------------------------------------------------------------
bool testInitializeSumTable() {
    std::vector<std::vector<double>> values = {{0, 1, 2, 3}, {1, 2, 3, 4}};
    SumTable<double> sumTable(values, false);

    // Manually computed expected sums
    std::vector<std::vector<double>> expectedSums = {
        {0, 0, 0, 0, 0}, {0, 0, 1, 3, 6}, {0, 1, 4, 9, 16}};
sumTable.sumWithinRectangle(PixelRect{0,0,0,0});
    for (int r = 1; r <= 1; ++r) {
        for (int c = 1; c <= 3; ++c) {
            if (sumTable.sumWithinRectangle(PixelRect{0,0,r,c}) != expectedSums[r][c]) { // Assuming sumTable member is made
                                                                 // public for testing purposes
                std::cout << "Actual: " << sumTable.sumTable[r][c]
                          << ", Expected: " << expectedSums[r][c] << std::endl;
                return false;
            }
        }
    }
    return true;
}

// TODO: M
bool test2() { return false; }

int main() {
    TestRunner runner;
    runner.addTest(testInitializeSumTable, "Test 1");
    runner.addTest(test2, "Test 2");
    runner.runTests();
}