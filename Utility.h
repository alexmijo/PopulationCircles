#include <iostream>


namespace utils {
void print(const std::string &line) { std::cout << line; }

void printLine(const std::string &line) {
    print(line);
    std::cout << "\n";
}
} // namespace utils