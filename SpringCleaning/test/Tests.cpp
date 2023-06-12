#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

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