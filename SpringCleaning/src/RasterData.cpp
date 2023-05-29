#pragma once
#include "EarthCircle.cpp"
#include "Utility.cpp"
#include <fstream>
#include <vector>

// Equirectangular raster data. Implemented using a summation table so it can get the sum of
// rectangular regions of the data in O(1) time.
class RasterData {
  public:
    // Returns the sum of the data in the specified rectangle.
    // TODO: Is PixelRectangle small enough that pass by value is faster?
    double sumWithinRectangle(const PixelRect &rect) {
        return sumTable[rect.south + 1][rect.east + 1] - sumTable[rect.north][rect.east + 1] -
               sumTable[rect.south + 1][rect.west] + sumTable[rect.north][rect.west];
    }

    // Returns the northernmost or southernmost row containing any pixel within <radius> km of the
    //  given row (specified by <y>).
    int boundingBoxEdge(const int y, const double radius, const bool isNorth) const {
        // Turn x and y (indices in the pop data) into geographic coordinates.
        const double cenLat = lat(y);

        // Start at center
        int edge = y;
        int edgeOfMap = kNumRows - 1;
        // Added to edge to move one pixel in desired direction
        int incOrDec = isNorth ? -1 : 1;

        while (edge >= 0 && edge <= edgeOfMap) {
            double currLat;
            currLat = lat(edge);
            // Since this function only works for north and south, longitude never changes
            const double distanceFromCen = distance({currLat, 0}, {cenLat, 0});
            if (distanceFromCen > radius) {
                edge -= incOrDec; // Went too far, walk it back
                break;
            } else {
                edge += incOrDec;
            }
        }
        if (edge < 0) {
            edge = 0;
        } else if (edge > edgeOfMap) {
            edge = edgeOfMap;
        }
        return edge;
    }

    double sumWithinCircle(const EarthCircle &circle) {}

    // Get longitude of the center of the <x>th (0 indexed) column.
    double lon(int x) const { return ((x + 0.5) / m_Width) * 360.0 - 180.0; }

    // Get latitude of the center of the <y>th (0 indexed) row
    double lat(int y) const { return ((y + 0.5) / m_Height) * 180.0 - 90.0; }

    // Get the latitude and longitude of the center of a pixel.
    // TODO: Is Pixel small enough that pass by value is faster?
    Location pixelToLocation(const Pixel &p) const { return {lat(p.y), lon(p.x)}; }

    // Get the column containing a longitude.
    int lonToX(double lon) const { return ((lon + 180.0) / 360.0) * m_Width; }

    // Get the row containing a latitude.
    int latToY(double lat) const { return ((lat + 90.0) / 180.0) * m_Height; }

    // Get the pixel containing a location.
    // TODO: Is util::Location small enough that pass by value is faster?
    Pixel locationToPixel(const Location &loc) const { return {lonToX(loc.lon), latToY(loc.lat)}; }

    RasterData() = default;
    explicit RasterData(const std::string &sumTableFilename) { loadData(sumTableFilename); }

    // File specified by <sumTableFilename> must consist of an int for height, an int for width, and
    // then height * width doubles representing all of the data in the summation table.
    // Projection must be equirectangular.
    void loadData(const std::string &sumTableFilename) {
        std::ifstream sumTableFile(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&m_Height), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&m_Width), sizeof(int));
        // The 0.0s will remain only as a column of 0.0s to the west and a row of 0.0s to the north.
        // This is padding for the data which will replace the rest of the zeros, and explains why
        // sumTable has dimensions one greater than both height and width.
        sumTable.resize(m_Height + 1, std::vector<double>(m_Width + 1, 0.0));
        for (int r = 1; r < m_Height + 1; r++) {
            // Read in an entire row at once
            sumTableFile.read(reinterpret_cast<char *>(&sumTable[r][1]), sizeof(double) * m_Width);
        }
        m_HasData = true;
    }

    int getNumRows() const { return m_Height; }

    int getNumCols() const { return m_Width; }

    bool hasData() const { return m_HasData; }

  private:
    // Dimensions of the original raster data.
    int m_Height = kNumRows;
    int m_Width = kNumCols;
    // First index is y (lat), second is x (lon). Matrix style indexing (top to bottom, left to
    //  right). Padded with 0s on the north and west, so it is one row and one column larger than
    //  the actual raster data. However, all methods of this class take indices into the original
    //  raster data (so one off in each dimension from the corresponding sumTable indices).
    // To be extra clear, the dimensions of sumTable are height + 1 by width + 1.
    std::vector<std::vector<double>> sumTable;
    // True iff the data has been loaded.
    bool m_HasData{false};
};