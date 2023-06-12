#include "PixelShapes.cpp"
#include "Utility.cpp"
#include <vector>

template <typename T> class SumTable {
  public:
    SumTable(const std::string &sumTableFilename) {
        std::ifstream sumTableFile(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&height), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&width), sizeof(int));
        // The 0s will remain only as a column of 0s to the west and a row of 0s to the south.
        // This is padding for the data which will replace the rest of the zeros, and explains why
        // sumTable has dimensions one greater than both height and width.
        sumTable.resize(height + 1, std::vector<T>(width + 1, 0));
        for (int r = 1; r < height + 1; r++) {
            // Read in an entire row at once
            sumTableFile.read(reinterpret_cast<char *>(&sumTable[r][1]), sizeof(T) * width);
        }
    }

    T sumWithinRectangle(const PixelRect &rect) {
        return sumTable[rect.N][rect.E] - sumTable[rect.S - 1][rect.E] -
               sumTable[rect.W - 1][rect.N] + sumTable[rect.W - 1][rect.S - 1];
    }

    T sumWithinCircle()

        int width;
    int height;

  private:
    std::vector<std::vector<T>> sumTable;
};
