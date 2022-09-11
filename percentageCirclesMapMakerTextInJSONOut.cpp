// Needs GDAL library
// This file is very spaghettiish right now, currently making it more like rasterCircleFinder.cpp,
//  and also making it so that their shared code isn't duplicated.

#include <cmath>
#include <cpl_conv.h>
#include <fstream>
#include <gdal.h>
#include <gdal_priv.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>

const double WORLD_POP = 7346242908.863955;
const int NUM_ROWS = 2500;
const int NUM_COLS = 5000;
constexpr double FLOAT_THRESHOLD = 0.000001;

class EquirectRasterData {

  private:
    const static int KERNEL_WIDTH = 4; // Num cols in each kernel (4 corners of a box)
    const static int SMALLEST_CIRCLES_VALUE_LENGTH = 3; // lon, lat, sum
    // key is radius, value is array holding lon, lat, sum
    std::map<int, double *> smallestCircleResults;
    std::string smallestCircleResultsFilename; // TODO: Make this a human readable text file instead
                                               //  of binary
    int numRows, numCols;
    // First index is y (lat), second is x (lon). Matrix style indexing (top to bottom, left to
    //  right)
    double **sumTable;

    // Used for boundingBoxEdge
    enum direction { north, south };

    // Returns the northmost or southmost row containing any pixel within <radius> km of the given
    //  pixel (specified by <x> and <y>).
    int boundingBoxEdge(const int x, const int y, const double radius, direction northOrSouth) {
        // Turn x and y (indices in the pop data) into geographic coordinates.
        const double cenLon = lon(x);
        const double cenLat = lat(y);

        // Start at center
        int edge = y;
        int edgeOfMap = numRows - 1;
        // Added to edge to move one pixel in desired direction. 0 = up, 1 = down
        int incOrDec;
        if (northOrSouth == north) {
            incOrDec = -1;
        } else if (northOrSouth == south) {
            incOrDec = 1;
        } else {
            // TODO: Throw an exception instead?
            std::cout << "boundingBoxEdge called with illegal direction" << std::endl;
        }

        while (edge >= 0 && edge <= edgeOfMap) {
            double currLat, currLon;
            currLat = lat(edge);
            currLon = cenLon; // Since this function only works for up (0) and down (1), lon never
                              //  changes
            const double distanceFromCen = distance(currLat, currLon, cenLat, cenLon);
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

    // For indexing into a kernel array
    inline int kernelIndex(const int row, const int col) { return KERNEL_WIDTH * row + col; }

    // Makes the kernel for a specific lattitude. Assigns the kernel's length to the reference
    //  parameter kernelLength.
    // TODO: Make the kernel an actual 2d array
    // TODO: See if this double counts part of 1 column when wrapping around all longitudes
    int *makeKernel(const int cenX, const int cenY, const double radius, int &kernelLength) {
        const double cenLon = lon(cenX);
        const double cenLat = lat(cenY);
        const int northEdge = boundingBoxEdge(cenX, cenY, radius, north);
        const int southEdge = boundingBoxEdge(cenX, cenY, radius, south);
        const int maxPossibleLength = southEdge - northEdge + 1;
        // A faster way of doing a 2D array of dimensions maxPossibleSize x KERNEL_WIDTH
        // Each row consists of: {westX, eastX, northY, southY} describing a summation table
        //  rectangle (so KERNEL_WIDTH must be 4) relative to cenX and cenY
        // Temp just cause we don't know length of real kernel yet
        int *tempKernel = new int[maxPossibleLength * KERNEL_WIDTH];

        int kernelRow = 0; // First index into kernel
        int y = northEdge;
        int horizontalOffset = 0; // From the verticle center line of the kernel
        while (y <= southEdge) {
            double currLon = lon(cenX + horizontalOffset);
            double currLat = lat(y);
            if (distance(currLat, currLon, cenLat, cenLon) > radius) {
                if (horizontalOffset == 0) {
                    // TODO: Probably some better way of logging/displaying this error
                    std::cout << "Something went wrong!1" << std::endl;
                }
                horizontalOffset--;
            } else {
                // See how much farther out we can go at this y level and still be in the circle
                while (distance(currLat, currLon, cenLat, cenLon) <= radius) {
                    horizontalOffset++;
                    if (horizontalOffset > numCols / 2) { // This rectangle wraps around the world
                        break;
                    }
                    currLon = lon(cenX + horizontalOffset);
                    currLat = lat(y);
                }
                horizontalOffset--; // horizontalOffset is now maximally far (after this decrement)

                // Start a new rectangle at y
                tempKernel[kernelIndex(kernelRow, 0)] = -horizontalOffset - 1;
                tempKernel[kernelIndex(kernelRow, 1)] = horizontalOffset;
                tempKernel[kernelIndex(kernelRow, 2)] = y - cenY - 1;
                if (y == southEdge) { // If we just started a new rectangle at southEdge, it must
                                      //  end there too
                    tempKernel[kernelIndex(kernelRow, 3)] = y - cenY;
                    break;
                } else { // Find where this new rectangle ends
                    while (true) {
                        y++;
                        if (y > southEdge) {
                            // Rectangle can't extend below south edge
                            tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1;
                            break;
                        }

                        // Check if the circle has widened
                        if (horizontalOffset < numCols / 2) { // Rectangles that wrap around the
                                                              //  whole world can't widen any more
                            currLon = lon(cenX + horizontalOffset + 1);
                            currLat = lat(y);
                            if (distance(currLat, currLon, cenLat, cenLon) <= radius) {
                                // The circle has widened; the rectangle is done
                                tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1;
                                kernelRow++;
                                break;
                            }
                        }

                        currLon = lon(cenX + horizontalOffset);
                        currLat = lat(y);
                        if (distance(currLat, currLon, cenLat, cenLon) > radius) {
                            // The y value can no longer be in the rectangle
                            tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1;
                            kernelRow++;
                            break;
                        }
                    }
                }
            }
        }

        kernelLength = kernelRow + 1; // Visible to the caller
        // Now that we know the length of kernel, we can construct it using the data in tempKernel
        //  and then return it
        int *kernel = new int[kernelLength * KERNEL_WIDTH];
        for (kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
            for (int kernelCol = 0; kernelCol < KERNEL_WIDTH; kernelCol++) {
                kernel[kernelIndex(kernelRow, kernelCol)] =
                    tempKernel[kernelIndex(kernelRow, kernelCol)];
            }
        }
        delete[] tempKernel;
        return kernel;
    }

    // Returns the total population in the specified rectangle. West and north are 1 pixel outside
    //  of the rectangle, east and south are 1 pixel inside of the rectangle.
    // TODO: Make this not use conditional logic by just padding the north and west sides of
    //  popSumTable with 0s
    // TODO: If this makes the program slow, see if making it inline helps
    double popWithinRectangle(const int west, const int east, const int north, const int south) {
        if (west == -1) {
            if (north == -1) {
                return sumTable[south][east]; // West side of rectangle on the antimeridian and
                                              //  north side of rectangle on the north pole
            }
            return sumTable[south][east] - sumTable[north][east]; // West side of rectangle on the
                                                                  //  antimeridian
        } else if (north == -1) { // North side of rectangle on the north pole
            return sumTable[south][east] - sumTable[south][west];
        }
        return sumTable[south][east] - sumTable[north][east] - sumTable[south][west] +
               sumTable[north][west];
    }

    // Returns the population of a circle with a radius of <radius> kilometers centered at the given
    //  latittude <cenLat> and longitude <cenLon>. TODO: Make this comment make more sense
    double popWithinKernel(const int cenX, const int cenY, int *kernel, const int kernelLength) {
        double totalPop = 0;

        for (int kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
            // Sides of the rectangle
            const int west = kernel[kernelIndex(kernelRow, 0)] + cenX;
            const int east = kernel[kernelIndex(kernelRow, 1)] + cenX;
            const int north = kernel[kernelIndex(kernelRow, 2)] + cenY;
            const int south = kernel[kernelIndex(kernelRow, 3)] + cenY;

            if (kernel[kernelIndex(kernelRow, 1)] == numCols / 2) { // This rectangle encircles the
                                                                    //  entire latitude
                totalPop += popWithinRectangle(-1, numCols - 1, north, south);
            } else if (west < -1) { // Need to wrap around the antimeridian, creating two rectangles
                totalPop += popWithinRectangle(-1, east, north, south);
                totalPop += popWithinRectangle(numCols + west, numCols - 1, north, south);
            } else if (east >= numCols) { // Need to wrap around the antimeridian, creating two
                                          //  rectangles
                totalPop += popWithinRectangle(west, numCols - 1, north, south);
                totalPop += popWithinRectangle(-1, east - numCols, north, south);
            } else {
                totalPop += popWithinRectangle(west, east, north, south);
            }
        }

        return totalPop;
    }

    // Get longitude of the center of the <x>th (0 indexed) column
    double lon(int x) { return (((double)x + 0.5) / (double)numCols) * 360.0 - 180.0; }

    // Get lattitude of the center of the <y>th (0 indexed) row
    double lat(int y) { return -((((double)y + 0.5) / (double)numRows) * 180.0 - 90.0); }

    // The inverse of the lon function
    int lonToX(double lon) { return ((lon + 180.0) / 360.0) * numCols - 0.5; }

    // The inverse of the lat function
    int latToY(double lat) { return ((-lat + 90.0) / 180.0) * numRows - 0.5; }

  public:
    // File specified by <sumTableFilename> must consist of an int for numRows, and int for
    //  numCols, and then numRows * numCols doubles representing all of the data in the summation
    //  table.
    // Projection must be equirectangular.
    // File specified by <smallestCircleResultsFilename> must consist of an int for
    //  numSmallestCircleResults, and then numSmallestCircleResults "rows" each consisting of 1 int
    //  followed by SMALLEST_CIRCLES_VALUE_LENGTH doubles.
    // TODO: Include some information in smallestCircleResultsFile to make sure it corresponds to
    //  this sumTableFile.
    // TODO: see if smallestCircleResultsFilename shouldn't be pass by reference
    EquirectRasterData(const std::string &sumTableFilename,
                       const std::string &smallestCircleResultsFilename)
        : smallestCircleResultsFilename(smallestCircleResultsFilename) {
        std::fstream sumTableFile;
        sumTableFile.open(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&numRows), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&numCols), sizeof(int));
        sumTable = new double *[numRows];
        for (int r = 0; r < numRows; r++) {
            sumTable[r] = new double[numCols];
            for (int c = 0; c < numCols; c++) {
                // TODO: See if it's faster to read in entire rows at a time
                sumTableFile.read(reinterpret_cast<char *>(&sumTable[r][c]), sizeof(double));
            }
        }
        sumTableFile.close();

        std::fstream smallestCircleResultsFile;
        smallestCircleResultsFile.open(smallestCircleResultsFilename,
                                       std::ios::in | std::ios::binary);
        // See if file is empty. If so, make it
        std::streampos begin, end;
        begin = smallestCircleResultsFile.tellg();
        smallestCircleResultsFile.seekg(0, std::ios::end);
        end = smallestCircleResultsFile.tellg();
        smallestCircleResultsFile.seekg(0, std::ios::beg);
        if (begin - end == 0) {
            std::cout << smallestCircleResultsFilename
                      << " is empty or doesn't exist. Making it "
                         "non-empty and existing."
                      << std::endl;
            smallestCircleResultsFile.close();
            // Make the file
            smallestCircleResultsFile.open(smallestCircleResultsFilename,
                                           std::ios::out | std::ios::binary);
            int numSmallestCircleResults = 0;
            smallestCircleResultsFile.write(reinterpret_cast<char *>(&numSmallestCircleResults),
                                            sizeof(int));
        } else {
            int numSmallestCircleResults;
            smallestCircleResultsFile.read(reinterpret_cast<char *>(&numSmallestCircleResults),
                                           sizeof(int));
            for (int i = 0; i < numSmallestCircleResults; i++) {
                int radius;
                smallestCircleResultsFile.read(reinterpret_cast<char *>(&radius), sizeof(int));
                // TODO: Don't really need this as a variable
                double *smallestCirclesValue = new double[SMALLEST_CIRCLES_VALUE_LENGTH];
                // TODO: See if I can just read in entire arrays at once
                // TODO: See if std::vector would make more sense here
                for (int j = 0; j < SMALLEST_CIRCLES_VALUE_LENGTH; j++) {
                    smallestCircleResultsFile.read(
                        reinterpret_cast<char *>(&smallestCirclesValue[j]), sizeof(double));
                }
                smallestCircleResults[radius] = smallestCirclesValue;
            }
        }
        smallestCircleResultsFile.close();
    }

    ~EquirectRasterData() {
        delete[] sumTable;
        for (std::map<int, double *>::iterator it = smallestCircleResults.begin();
             it != smallestCircleResults.end(); it++) {
            delete[](it->second);
        }
    }

    // Returns distance in kilometers between two points on Earth's surface
    static double distance(const double &lat1, const double &lon1, const double &lat2,
                           const double &lon2) {
        double req = 6378137.0;

        // alexmijo added code for dealing with equatorial points or identical points
        if (lat1 == lat2 && lon1 == lon2) {
            return 0.0;
        } else if (lat1 == 0.0 && lat2 == 0.0) {
            double lonDiff = abs(lon1 - lon2);
            double fractionOfEquator = lonDiff / 360.0;
            double equatorLength = 2.0 * 3.14159265358979323 * req;
            return fractionOfEquator * equatorLength;
        }
        // The rest of this function (except for the conversion to kilometers) was copy pasted from
        //  somewhere (I forgot where) online.

        const double latitude_01 = lat1 * M_PI / 180.0;
        const double longitude_01 = lon1 * M_PI / 180.0;

        const double latitude_02 = lat2 * M_PI / 180.0;
        const double longitude_02 = lon2 * M_PI / 180.0;

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

        const int maxIterations = 20000;
        for (int i = 0; i < maxIterations; i++) {

            //
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

            if (i == maxIterations - 1) {
                return 100000; // Failed to converge, nearly antipodal points
            }
        }

        // U2
        double u_sq = cos_alpha_sq * (a * a - b * b) / (b * b);

        // Compute A, B
        double A = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));

        double B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))));

        // Sigma
        double cos_2sigma_m_sq = cos_2sigma_m * cos_2sigma_m;
        double delta_sigma = B * sin_sigma *
                             (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq) -
                                                          (B / 6.0) * cos_2sigma_m *
                                                              (-3 + 4 * sin_sigma * sin_sigma) *
                                                              (-3 + 4 * cos_2sigma_m_sq)));

        // Distance
        double s = b * A * (sigma - delta_sigma);

        // Convert from meters to kilometers before returnings
        return s / 1000.;
    }

    int getNumRows() { return numRows; }

    int getNumCols() { return numCols; }

    // TODO: Check if the following two functions actually give the correct answers
    //  (by brute forcing some circles with the original pop data)

    // Finds the circle of the given radius that maximizes the sum of the data inside it. Returns
    //  a pointer to an array containing the longitude [0] and lattitude [1] of the center of the
    //  circle and the sum of the data inside the circle [2].
    // Can optionally constrain the center of the circle to be within given ranges of latitudes and
    //  longitudes with leftLon, rightLon, upLat and downLat.
    // The parameter initialLargestSum can speed up computation if it is passed in. Must be for sure
    //  known to be at least as small as what the largestSum will end up being.
    // Radius given in kilometers.
    double *largestSumCircleOfGivenRadius(const double radius, const double leftLon = -180,
                                          const double rightLon = 180, const double upLat = 90,
                                          const double downLat = -90,
                                          const double initialLargestSum = 0) {
        const int smallStep = 1;                                 // 1
        const int mediumStep = std::max((int)(radius / 128), 1); // 4
        const int largeStep = std::max((int)(radius / 32), 1);   // 16
        const int xLStep = std::max((int)(radius / 8), 1);       // 64
        const int xXLStep = std::max((int)(radius / 2), 1);      // 256
        const int printMod = 163; // Print lattitude when the Y index mod this is 0

        int step = smallStep;
        // Turn lat and lon into indices in the pop data.
        const int leftX = lonToX(leftLon);
        const int rightX = lonToX(rightLon);
        const int upY = latToY(upLat);
        const int downY = latToY(downLat);

        double largestSumCenLon;
        double largestSumCenLat;
        double largestSum = initialLargestSum;

        for (int cenY = upY; cenY <= downY; cenY += step) {
            if (cenY % printMod == 0) {
                std::cout << "Current lattitude: " << lat(cenY) << std::endl;
            }

            int kernelLength;
            int *kernel = makeKernel(1000, cenY, radius, kernelLength); // Initializes kernelLength

            for (int cenX = leftX; cenX <= rightX; cenX += step) {
                double popWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength);
                if (popWithinNKilometers >= largestSum) {
                    std::cout << "Sum within " << radius << " kilometers of (" << lat(cenY) << ", "
                              << lon(cenX) << "): " << ((long long)popWithinNKilometers)
                              << std::endl;
                    largestSumCenLon = lon(cenX);
                    largestSumCenLat = lat(cenY);
                    largestSum = popWithinNKilometers;
                }

                if (popWithinNKilometers > largestSum * 0.82) {
                    if (step > smallStep) {
                        cenX -= step;
                    }
                    step = smallStep;
                } else if (popWithinNKilometers > largestSum * 0.65) {
                    if (step > mediumStep) {
                        cenX -= step;
                    }
                    step = mediumStep;
                } else if (popWithinNKilometers > largestSum * 0.47) {
                    if (step > largeStep) {
                        cenX -= step;
                    }
                    step = largeStep;
                } else if (popWithinNKilometers > largestSum * 0.32) {
                    if (step > xLStep) {
                        cenX -= step;
                    }
                    step = xLStep;
                } else {
                    step = xXLStep;
                }
            }
            delete[] kernel;
            step = smallStep;
        }

        double *returnValues = new double[3];
        returnValues[0] = largestSumCenLon;
        returnValues[1] = largestSumCenLat;
        returnValues[2] = largestSum;
        return returnValues;
    }

    // Finds the circle of the given radius that minimizes the sum of the data inside it. Returns
    //  a pointer to an array containing the longitude [0] and lattitude [1] of the center of the
    //  circle and the sum of the data inside the circle [2].
    // Can optionally constrain the center of the circle to be within given ranges of latitudes and
    //  longitudes with leftLon, rightLon, upLat and downLat.
    // The parameter initialSmallestSum can speed up computation if it is passed in. Must be for
    // sure known to
    //  be at least as large as what the smallestSum will end up being.
    // Radius given in kilometers.
    double *smallestSumCircleOfGivenRadius(const double radius, const double leftLon = -180,
                                           const double rightLon = 180, const double upLat = 90,
                                           const double downLat = -90,
                                           const double initialSmallestSum = 1000000000000000) {
        const int smallStep = 1;                                 // 1
        const int mediumStep = std::max((int)(radius / 128), 1); // 4
        const int largeStep = std::max((int)(radius / 32), 1);   // 16
        const int xLStep = std::max((int)(radius / 8), 1);       // 64
        const int xXLStep = std::max((int)(radius / 2), 1);      // 256
        const int printMod = 163; // Print lattitude when the Y index mod this is 0

        int step = smallStep;
        // Turn lat and lon into indices in the pop data.
        const int leftX = lonToX(leftLon);
        const int rightX = lonToX(rightLon);
        const int upY = latToY(upLat);
        const int downY = latToY(downLat);

        double smallestSumCenLon;
        double smallestSumCenLat;
        double smallestSum = initialSmallestSum;

        for (int cenY = upY; cenY <= downY; cenY += step) {
            if (cenY % printMod == 0) {
                std::cout << "Current lattitude: " << lat(cenY) << std::endl;
            }

            int kernelLength;
            int *kernel = makeKernel(1000, cenY, radius, kernelLength); // Initializes kernelLength

            for (int cenX = leftX; cenX <= rightX; cenX += step) {
                double sumWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength);
                if (sumWithinNKilometers <= smallestSum) {
                    std::cout << "Sum within " << radius << " kilometers of (" << lat(cenY) << ", "
                              << lon(cenX) << "): " << ((long long)sumWithinNKilometers)
                              << std::endl;
                    smallestSumCenLon = lon(cenX);
                    smallestSumCenLat = lat(cenY);
                    smallestSum = sumWithinNKilometers;
                }

                if (sumWithinNKilometers < smallestSum * 1.2) {
                    if (step > smallStep) {
                        cenX -= step;
                    }
                    step = smallStep;
                } else if (sumWithinNKilometers < smallestSum * 1.4) {
                    if (step > mediumStep) {
                        cenX -= step;
                    }
                    step = mediumStep;
                } else if (sumWithinNKilometers < smallestSum * 1.6) {
                    if (step > largeStep) {
                        cenX -= step;
                    }
                    step = largeStep;
                } else if (sumWithinNKilometers < smallestSum * 1.8) {
                    if (step > xLStep) {
                        cenX -= step;
                    }
                    step = xLStep;
                } else {
                    step = xXLStep;
                }
            }
            delete[] kernel;
            step = smallStep;
        }

        double *returnValues = new double[3];
        returnValues[0] = smallestSumCenLon;
        returnValues[1] = smallestSumCenLat;
        returnValues[2] = smallestSum;
        return returnValues;
    }

    // TODO: Make a spec comment for this
    // Kind of dangerous to pass in the last 4 args to this function, since the
    //  smallestCircleResultsFile might get spurious data added to it then
    double *smallestCircleWithGivenSum(const double sum, const double leftLon = -180,
                                       const double rightLon = 180, const double upLat = 90,
                                       const double downLat = -90) {
        // Radii will be ints cause only interested in getting to the nearest kilometer
        const int EQUATOR_LEN = 40075;
        int upperBound = EQUATOR_LEN / 2;
        int lowerBound = 0;
        double *returnValues;
        bool returnValuesAssigned = false; // Pointer not yet pointing to an array
        double initialLargestSum = 0;      // To be passed into largestSumCircleOfGivenRadius()

        // Tighten bounds as much as possible using previous results
        for (std::map<int, double *>::iterator it = smallestCircleResults.begin();
             it != smallestCircleResults.end(); it++) {
            if ((it->second)[2] < sum) {
                lowerBound = it->first;
                initialLargestSum = (it->second)[2];
            } else if ((it->second)[2] >= sum) {
                upperBound = it->first;
                returnValues = new double[4];
                returnValues[0] = (it->second)[0];
                returnValues[1] = (it->second)[1];
                returnValues[2] = (it->second)[2];
                returnValues[3] = it->first;
                returnValuesAssigned = true;
                break;
            }
        }

        int radius = lowerBound + (upperBound - lowerBound) / 2; // Start of binary search
        while (upperBound - lowerBound > 1) {
            double narrowLeftLon = leftLon;
            double narrowRightLon = rightLon;
            double narrowUpLat = upLat;
            double narrowDownLat = downLat;
            // May be able to use knowledge from past map results (the ones on reddit) to narrow the
            //  search area
            if (radius < 117) {
                // Uncharted territory, don't narrow search area
            } else if (radius <= 8000) {
                narrowLeftLon = 42;
                narrowRightLon = 135;
                narrowUpLat = 53;
                narrowDownLat = 11;
            } else if (radius < 15000) {
                narrowDownLat = -53;
            }

            double *largestSumCircle =
                largestSumCircleOfGivenRadius(radius, narrowLeftLon, narrowRightLon, narrowUpLat,
                                              narrowDownLat, initialLargestSum);
            if (largestSumCircle[2] >= sum) {
                upperBound = radius;
                if (returnValuesAssigned) {
                    delete[] returnValues;
                }
                returnValues = new double[4];
                returnValues[0] = largestSumCircle[0];
                returnValues[1] = largestSumCircle[1];
                returnValues[2] = largestSumCircle[2];
                returnValues[3] = radius;
                returnValuesAssigned = true;
            } else {
                lowerBound = radius;
                initialLargestSum = largestSumCircle[2];
            }

            // Add result to smallestCircleResults and its file
            // TODO: Combine this line and the next one
            std::map<int, double *>::iterator it = smallestCircleResults.find(radius);
            if (it != smallestCircleResults.end()) {
                // TODO: Figure out a better way to do these sort of things (probably throw an
                //  exception)
                // TODO: See if there's a way to send a string to both streams
                std::cout << "smallestCircleWithGivenSum had an error involving "
                             "smallestCircleResults"
                          << std::endl;
                std::cerr << "smallestCircleWithGivenSum had an error involving "
                             "smallestCircleResults"
                          << std::endl; // Above line might get drowned out
            } else {
                std::fstream smallestCircleResultsFile;
                smallestCircleResultsFile.open(smallestCircleResultsFilename,
                                               std::ios::in | std::ios::out | std::ios::binary);
                int numSmallestCircleResults;
                smallestCircleResultsFile.read(reinterpret_cast<char *>(&numSmallestCircleResults),
                                               sizeof(int));
                numSmallestCircleResults++;
                smallestCircleResultsFile.seekg(0, std::ios::beg);
                smallestCircleResultsFile.write(reinterpret_cast<char *>(&numSmallestCircleResults),
                                                sizeof(int));
                smallestCircleResultsFile.seekg(0, std::ios::end);
                smallestCircleResultsFile.write(reinterpret_cast<char *>(&radius), sizeof(int));
                // TODO: See if I can just write the entire array at once
                for (int j = 0; j < SMALLEST_CIRCLES_VALUE_LENGTH; j++) {
                    smallestCircleResultsFile.write(reinterpret_cast<char *>(&largestSumCircle[j]),
                                                    sizeof(double));
                }
                smallestCircleResultsFile.close();
                // TODO: Instead of this workaround, just don't delete[] largestSumCircle in this
                //  case
                smallestCircleResults[radius] = new double[SMALLEST_CIRCLES_VALUE_LENGTH];
                smallestCircleResults[radius][0] = largestSumCircle[0];
                smallestCircleResults[radius][1] = largestSumCircle[1];
                smallestCircleResults[radius][2] = largestSumCircle[2];
            }

            delete[] largestSumCircle;
            radius = lowerBound + (upperBound - lowerBound) / 2; // Binary search
        }
        if (!returnValuesAssigned) {
            // TODO: Figure out a better way to do these sort of things (probably throw an
            //  exception)
            std::cout << "smallestCircleWithGivenSum wasn't able to find a circle with a large "
                         "enough sum. Either the sum ("
                      << sum << ") was too large, or a bug occurred." << std::endl;
        }
        return returnValues;
    }
};

// Get longitude of the center of the <x>th (0 indexed) column
double lon(int x) { return (((double)x + 0.5) / (double)NUM_COLS) * 360.0 - 180.0; }

// Get lattitude of the center of the <y>th (0 indexed) row
double lat(int y) { return -((((double)y + 0.5) / (double)NUM_ROWS) * 180.0 - 90.0); }

int main() {
    // TODO: Make this a function taking the percent as a parameter.
    double percent = 99.8; // Between 0.1 and 100 (inclusive). Must already have been found.
    double cenLat;
    double cenLon;
    double radius;
    bool centerCoordsAndRadiusInitialized = false;

    std::string foundPercentageCirclesFilename = "foundPercentageCircles2020.txt";
    std::ifstream foundPercentageCircles;
    foundPercentageCircles.open(foundPercentageCirclesFilename);
    std::string percentageCircleString;
    while (getline(foundPercentageCircles, percentageCircleString)) {
        std::stringstream percentageCircleSS(percentageCircleString);
        std::string percentString;
        getline(percentageCircleSS, percentString, ' ');
        if (std::abs(std::stod(percentString) - percent) < FLOAT_THRESHOLD) {
            std::string radiusString;
            // TODO: See how to get bool from these getlines to see if file is formatted correctly
            getline(percentageCircleSS, radiusString, ' ');
            radius = std::stoi(radiusString);
            std::string lonString;
            getline(percentageCircleSS, lonString, ' ');
            cenLon = std::stod(lonString);
            std::string latString;
            getline(percentageCircleSS, latString, ' ');
            cenLat = std::stod(latString);
            centerCoordsAndRadiusInitialized = true;
        }
    }
    foundPercentageCircles.close();
    if (!centerCoordsAndRadiusInitialized) {
        std::cout << "Desired percentage circle wasn't in " + foundPercentageCirclesFilename
                  << std::endl;
        return 1;
    }

    // TODO: Make a JSON-ize function for arrays, or find one
    // int output[NUM_ROWS][NUM_COLS];

    std::string landNBordersFileName = "landNBorders.txt";
    std::ifstream landNBorders;
    landNBorders.open(landNBordersFileName);
    std::ofstream colorsJSON;
    std::string colorsJSONFilename =
        "ColorsJSONFiles2020/colorsJSON" +
        std::to_string(percent).substr(0, std::to_string(percent).size() - 5) + ".txt";
    colorsJSON.open(colorsJSONFilename);
    colorsJSON << "[";
    for (int r = 0; r < NUM_ROWS; r++) {
        if (r % 10 == 0) {
            std::cout << r << std::endl;
        }
        std::string rString;
        getline(landNBorders, rString);
        colorsJSON << "[";
        for (int c = 0; c < NUM_COLS; c++) {
            double distance = EquirectRasterData::distance(cenLat, cenLon, lat(r), lon(c));
            char landOrBorder = rString.at(c);
            if (distance <= 53 && distance <= radius) {
                // TODO: Make color int constants
                // Dark red
                colorsJSON << '3';
            } else if (landOrBorder == '0') {
                // Black
                colorsJSON << '0';
            } else if (landOrBorder == '1') {
                if (distance <= radius) {
                    // Red
                    colorsJSON << '4';
                } else {
                    // Grey
                    colorsJSON << '1';
                }
            } else if (landOrBorder == '2') {
                if (distance <= radius && distance > radius - 80) {
                    // Light red
                    colorsJSON << '5';
                } else {
                    // White
                    colorsJSON << '2';
                }
            } else {
                std::cout << "Illegal LandOrBorder file!" << std::endl;
            }

            if (c != NUM_COLS - 1) {
                colorsJSON << ", ";
            }
        }
        colorsJSON << "]";
        if (r != NUM_ROWS - 1) {
            colorsJSON << ", ";
        }
    }
    colorsJSON << "]";
    colorsJSON.close();
    landNBorders.close();
}