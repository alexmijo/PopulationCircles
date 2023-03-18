#include "Utility.cpp"
#include <fstream>

constexpr int k30ArcSecsPerDegree = 2 * 60;
constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
constexpr int kNumRows = 180 * k30ArcSecsPerDegree;

struct Pixel {
    // x from left to right, y from bottom to top
    int x, y;
};

// Defines a rectangular region of the raster data. Ranges are inclusive.
class PixelRectangle {
  public:
    int west, east;
    int south, north;

    bool contains(Pixel p) const {
        return west <= p.x && p.x <= east && north <= p.y && p.y <= south;
    }
};

// Equirectangular raster data. Implemented using a summation table so it can get the sum of
// rectangular regions of the data in O(1) time.
class RasterData {
  public:
    // Returns the sum of the data in the specified rectangle.
    double sumWithinRectangle(PixelRectangle rect) {
        return sumTable[rect.south + 1][rect.east + 1] - sumTable[rect.north][rect.east + 1] -
               sumTable[rect.south + 1][rect.west] + sumTable[rect.north][rect.west];
    }

    // Get longitude of the center of the <x>th (0 indexed) column.
    double lon(int x) const { return ((x + 0.5) / width) * 360.0 - 180.0; }

    // Get lattitude of the center of the <y>th (0 indexed) row
    double lat(int y) const { return ((y + 0.5) / height) * 180.0 - 90.0; }

    // Get the latitude and longitude of the center of a pixel.
    util::Location pixelToLocation(Pixel p) const { return {lat(p.y), lon(p.x)}; }

    // Get the column containing a longitude.
    int lonToX(double lon) const { return ((lon + 180.0) / 360.0) * width; }

    // Get the row containing a lattitude.
    int latToY(double lat) const { return ((lat + 90.0) / 180.0) * height; }

    // Get the pixel containing a location.
    Pixel locationToPixel(util::Location loc) const { return {lonToX(loc.lon), latToY(loc.lat)}; }

    RasterData() {}
    RasterData(const std::string &sumTableFilename) { loadData(sumTableFilename); }

    // File specified by <sumTableFilename> must consist of an int for height, and int for
    //  width, and then height * width doubles representing all of the data in the summation
    //  table.
    // Projection must be equirectangular.
    // TODO: Complete spec
    void loadData(const std::string &sumTableFilename) {
        std::fstream sumTableFile;
        sumTableFile.open(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&height), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&width), sizeof(int));
        // TODO: Either use a unique_pointer to array of unique_pointers to arrays of ints, a vector
        //  of vectors, or the 1D versions of either of those (with an indexing function) for speed.
        //  No need for new and delete.
        sumTable = new double *[height + 1];
        for (int r = 0; r < height + 1; r++) {
            sumTable[r] = new double[width + 1];
            if (r == 0) {
                // Pad with a row of 0s to the north
                for (int c = 0; c < width + 1; c++) {
                    sumTable[0][c] = 0;
                }
                continue;
            }
            // Pad with a column of 0s to the west
            sumTable[r][0] = 0;
            // Reads in an entire row at once
            sumTableFile.read(reinterpret_cast<char *>(&sumTable[r][1]), sizeof(double) * width);
        }
        mustFree = true;
    }

    ~RasterData() {
        if (mustFree) {
            for (int r = 0; r < height + 1; r++) {
                delete[] sumTable[r];
            }
            delete[] sumTable;
        }
    }

    int getNumRows() const { return height; }

    int getNumCols() const { return width; }

    bool hasData() const { return mustFree; }

  private:
    // Dimensions of the original raster data.
    int height = kNumRows;
    int width = kNumCols;
    // First index is y (lat), second is x (lon). Matrix style indexing (top to bottom, left to
    //  right). Padded with 0s on the north and west, so it is one row and one column larger than
    //  the actual raster data. However, all methods of this class take indices into the original
    //  raster data (so one off in each dimension from the corresponding sumTable indices).
    // To be extra clear, the dimensions of sumTable are height + 1 by width + 1.
    // TODO: See if I can or should make this an array of arrays (library arrays). Can't be a 2D C
    //  array cause it's too large.
    // TODO: See if I should actually make this a 1D array (both for loading speed and using speed)
    double **sumTable;
    // True iff destructor must actually free resources
    bool mustFree{false};
};