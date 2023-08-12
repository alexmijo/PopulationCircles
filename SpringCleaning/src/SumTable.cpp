#include "Utility.cpp"

#include <fstream>
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

    // TODO: Something with the padding?
    // TODO: Move?
    SumTable(const std::vector<std::vector<T>> &sumTableOrValues, const bool isSumTable) {
        if (isSumTable) {
            sumTable = sumTableOrValues;
            height = sumTableOrValues.size() - 1;
            width = sumTableOrValues[0].size();
        } else {
            height = sumTableOrValues.size();
            width = height > 0 ? sumTableOrValues[0].size() : 0;

            // Initialize the sumTable with an extra row and column for padding
            sumTable.resize(height + 1, std::vector<T>(width + 1, 0));

            // Calculate the sumTable
            for (int r = 1; r <= height; ++r) {
                for (int c = 1; c <= width; ++c) {
                    sumTable[r][c] = sumTableOrValues[r - 1][c - 1] + sumTable[r - 1][c] +
                                     sumTable[r][c - 1] - sumTable[r - 1][c - 1];
                    std::cout << "[" << r << "][" << c << "]: " << sumTable[r][c] << std::endl;
                }
            }
        }
    }

    T sumWithinRectangle(const PixelRect &rect) {
        return sumTable[rect.N][rect.E] - sumTable[rect.S - 1][rect.E] -
               sumTable[rect.W - 1][rect.N] + sumTable[rect.W - 1][rect.S - 1];
    }

    int width;
    int height;

    //   private:
    std::vector<std::vector<T>> sumTable;
};
