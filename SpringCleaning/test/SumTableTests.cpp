#include "Tests.cpp"
#include "../src/SumTable.cpp"

void createSumTable

//--------------------------------------------------------------------------------------------------
bool testInitializeSumTable() {
    return true;
}

bool test2() {
    return false;
}

int main() {
    TestRunner runner;
    runner.addTest(test1, "Test 1");
    runner.addTest(test2, "Test 2");
    runner.runTests();
}