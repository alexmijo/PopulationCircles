#include <iostream>
#include <unordered_map>
#include <vector>

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

int main() {
    std::vector<std::vector<std::string>> image1 = {
        {"r", "r", "r", "r", "g", "b"}, {"r", "r", "r", "r", "g", "b"},
        {"r", "r", "r", "r", "g", "b"}, {"r", "r", "r", "r", "g", "b"},
        {"r", "r", "r", "r", "g", "b"}, {"r", "r", "r", "r", "g", "b"}};

    makePpm(image1);

    return 0;
}
