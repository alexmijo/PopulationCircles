#include <cmath>
#include <fstream>

constexpr int k30ArcSecsPerDegree = 2 * 60;
constexpr int kNumCols = 360 * k30ArcSecsPerDegree;
constexpr int kNumRows = 180 * k30ArcSecsPerDegree;

struct Pixel {
    // x from left to right, y from bottom to top
    int x, y;
};

struct Location {
    double lat, lon;
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
    Location pixelToLocation(Pixel p) const { return {lat(p.y), lon(p.x)}; }

    // Get the column containing a longitude.
    int lonToX(double lon) const { return ((lon + 180.0) / 360.0) * width; }

    // Get the row containing a lattitude.
    int latToY(double lat) const { return ((lat + 90.0) / 180.0) * height; }

    // Get the pixel containing a location.
    Pixel locationToPixel(Location loc) const { return {lonToX(loc.lon), latToY(loc.lat)}; }

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

// TODO: Put this in the circle stuff file or in its own file.
// Returns distance in kilometers between two points on Earth's surface.
// TODO: Find source and give credit
static double distance(Location loc1, Location loc2) {
    double req = 6378137.0;

    // alexmijo added code for dealing with equatorial points or identical points
    if (loc1.lat == loc2.lat && loc1.lon == loc2.lon) {
        return 0.0;
    } else if (loc1.lat == 0.0 && loc2.lat == 0.0) {
        double lonDiff = abs(loc1.lon - loc2.lon);
        double fractionOfEquator = lonDiff / 360.0;
        double equatorLength = 2.0 * 3.14159265358979323 * req;
        return fractionOfEquator * equatorLength;
    }

    const double latitude_01 = loc1.lat * M_PI / 180.0;
    const double longitude_01 = loc1.lon * M_PI / 180.0;

    const double latitude_02 = loc2.lat * M_PI / 180.0;
    const double longitude_02 = loc2.lon * M_PI / 180.0;

    const double a = 6378137.0;
    const double b = 6356752.314245;

    // Flattening
    const double f = (a - b) / a;

    // tan U1
    const double tan_U1 = (1 - f) * std::tan(latitude_01);
    const double tan_U2 = (1 - f) * std::tan(latitude_02);

    // Longitudinal Distance
    const double cos_U1 = 1 / std::sqrt(1 + tan_U1 * tan_U1);
    const double cos_U2 = 1 / std::sqrt(1 + tan_U2 * tan_U2);
    const double sin_U1 = tan_U1 * cos_U1;
    const double sin_U2 = tan_U2 * cos_U2;

    // Iterate until complete
    const double L = longitude_02 - longitude_01;
    double lambda = L;
    double diff, sigma;
    double cos_alpha_sq, cos_2sigma_m;
    double cos_sigma, sin_sigma;

    while (true) {

        double sin_lambda = std::sin(lambda);
        double cos_lambda = std::cos(lambda);

        double c1 = (cos_U2 * sin_lambda) * (cos_U2 * sin_lambda);
        double c2 = (cos_U1 * sin_U2);
        double c3 = (sin_U1 * cos_U2 * cos_lambda);

        //  sin sigma
        sin_sigma = std::sqrt(c1 + (c2 - c3) * (c2 - c3));

        // cos sigma
        cos_sigma = sin_U1 * sin_U2 + cos_U1 * cos_U2 * cos_lambda;

        // sigma
        sigma = std::atan2(sin_sigma, cos_sigma);

        // sin alpha
        double sin_alpha = (cos_U1 * cos_U2 * sin_lambda) / (sin_sigma);

        // cos^2 alpha
        cos_alpha_sq = 1 - (sin_alpha * sin_alpha);

        // cos^2 2sigmam
        cos_2sigma_m = cos_sigma - (2 * sin_U1 * sin_U2) / (cos_alpha_sq);

        // C
        double C = (f / 16.0) * cos_alpha_sq * (4 + f * (4 - 3 * cos_alpha_sq));

        // Update Lambda
        diff = lambda;
        lambda = L + (1 - C) * f * sin_alpha *
                         (sigma + C * sin_sigma *
                                      (cos_2sigma_m +
                                       C * cos_sigma * (-1 + 2 * cos_2sigma_m * cos_2sigma_m)));
        diff = lambda - diff;
        if (std::fabs(diff) < 0.00001) {
            break;
        }
    }

    // U2
    double u_sq = cos_alpha_sq * (a * a - b * b) / (b * b);

    // Compute A, B
    double A = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));

    double B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))));

    // Sigma
    double cos_2sigma_m_sq = cos_2sigma_m * cos_2sigma_m;
    double delta_sigma =
        B * sin_sigma *
        (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq) -
                                     (B / 6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma * sin_sigma) *
                                         (-3 + 4 * cos_2sigma_m_sq)));

    // Distance
    double s = b * A * (sigma - delta_sigma);

    // Convert from meters to kilometers before returning
    return s / 1000.;
}