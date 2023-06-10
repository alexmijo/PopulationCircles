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
            delete[] (it->second);
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

// TODO: Spec
void oldMain(int whichOne) {
    std::string colorsJSONFilename = std::to_string(whichOne) + ".txt";
    double radius = 1000;
    std::vector<double> centerLats = {
        27.7172,    28.613895,  23.763889,  27.472222,  39.906667,  33.693056,  39.019444,
        37.56,      25.066667,  19.7475,    21.028333,  34.525278,  49.611667,  47.141,
        46.948056,  48.2,       46.051389,  17.966667,  48.143889,  52.52,      50.0875,
        50.846667,  52.372778,  43.731111,  48.856613,  45.816667,  47.4925,    9.066667,
        43.9346,    6.497222,   44.817778,  51.507222,  43.856389,  38.536667,  6.131944,
        33.886944,  55.676111,  42.5,       13.515,     33.513056,  52.23,      4.85,
        35.1725,    -1.286389,  3.752064,   41.9025,    41.893333,  0.313611,   42.441286,
        5.55,       31.949722,  -3.428333,  3.866667,   31.778889,  31.778889,  13.7525,
        47.022778,  42.663333,  -1.943889,  39.93,      44.4325,    -6.173056,  33.315278,
        42.7,       41.996111,  12.368611,  -6.175,     41.328889,  6.934444,   9.03,
        54.687222,  53.9,       50.45,      15.322778,  40.181389,  11.569444,  30.044444,
        11.588333,  38.904722,  1.283333,   56.948889,  0.336111,   41.7225,    35.689167,
        41.311111,  6.816111,   53.35,      0.390278,   12.11,      40.416667,  29.369722,
        36.753889,  37.984167,  40.395278,  36.806389,  35.689722,  12.639167,  19.433333,
        3.147778,   55.755833,  -17.829167, 15.348333,  -13.983333, 45.424722,  -15.793889,
        15.500556,  37.9375,    14.5958,    42.874722,  4.175278,   59.437222,  59.329444,
        59.913333,  -15.416667, 34.020882,  60.170833,  26.225,     38.725267,  24.633333,
        35.898333,  4.711111,   6.313333,   -25.746111, 8.484444,   9.509167,   -24.658056,
        -25.966667, -26.416667, 25.286667,  14.613333,  -4.269444,  17.251389,  8.983333,
        -4.325,     4.373333,   2.039167,   -11.699,    -29.31,     -25.3,      13.698889,
        24.466667,  10.480556,  -0.22,      32.887222,  11.85,      14.1,       18.466667,
        12.136389,  23.588889,  13.453056,  17.971389,  25.078056,  -34.603333, -8.838333,
        -34.883611, 9.9325,     18.533333,  14.692778,  23.136667,  18.08581,   17.3,
        -18.933333, 4.890278,   -12.06,     -33.45,     12.05,      10.666667,  51.166667,
        15.301389,  13.157778,  17.116667,  14.016667,  -8.553611,  13.0975,    -19.0475,
        47.920278,  14.918,     -35.293056, 6.805833,   -22.57,     -9.478889,  5.852222,
        7.500556,   -20.164444, -41.288889, -21.133333, -9.431944,  -18.1416,   -17.733333,
        -13.833333, 64.146667,  -8.516667,  7.1167,     1.433333,   -0.5477,    -4.6167,
        6.917222};
    std::vector<double> centerLons = {
        85.324,     77.209006,  90.388889,  89.636111,  116.3975,   73.063889,  125.738056,
        126.99,     121.516667, 96.115,     105.854167, 69.178333,  6.131944,   9.521,
        7.4475,     16.366667,  14.506111,  102.6,      17.109722,  13.405,     14.421389,
        4.3525,     4.893611,   7.42,       2.352222,   15.983333,  19.051389,  7.483333,
        12.4473,    2.605,      20.456944,  -0.1275,    18.413056,  68.78,      1.222778,
        35.513056,  12.568333,  1.5,        2.1175,     36.291944,  21.011111,  31.6,
        33.365,     36.817222,  8.7737,     12.4525,    12.482778,  32.581111,  19.262892,
        -0.2,       35.932778,  29.925,     11.516667,  35.225556,  35.225556,  100.494167,
        28.835278,  21.162222,  30.059444,  32.85,      26.103889,  35.741944,  44.366111,
        23.33,      21.431667,  -1.5275,    106.8275,   19.817778,  79.842778,  38.74,
        25.28,      27.566667,  30.523333,  38.925,     44.514444,  104.921111, 31.235833,
        43.145,     -77.016389, 103.833333, 24.106389,  6.730556,   44.7925,    51.388889,
        69.279722,  -5.274167,  -6.260278,  9.454167,   15.05,      -3.7025,    47.978333,
        3.058889,   23.728056,  49.882222,  10.181667,  139.692222, -8.002778,  -99.133333,
        101.695278, 37.617222,  31.052222,  44.206389,  33.783333,  -75.695,    -47.882778,
        32.56,      58.38,      120.9772,   74.612222,  73.508889,  24.745278,  18.068611,
        10.738889,  28.283333,  -6.84165,   24.9375,    50.5775,    -9.150019,  46.716667,
        14.5125,    -74.072222, -10.801389, 28.188056,  -13.234444, -13.712222, 25.912222,
        32.583333,  31.166667,  51.533333,  -90.535278, 15.271389,  -88.766944, -79.516667,
        15.322222,  18.562778,  45.341944,  43.256,     27.48,      -57.633333, -89.191389,
        54.366667,  -66.903611, -78.5125,   13.191389,  -15.566667, -87.216667, -69.95,
        -86.251389, 58.408333,  -16.5775,   -76.793056, -77.338611, -58.381667, 13.234444,
        -56.181944, -84.08,     -72.333333, -17.446667, -82.358889, -15.9785,   -62.733333,
        47.516667,  114.942222, -77.0375,   -70.666667, -61.75,     -61.516667, 71.433333,
        -61.388333, -61.225,    -61.85,     -60.983333, 125.578333, -59.616667, -65.26,
        106.917222, -23.509,    149.126944, -58.150833, 17.083611,  147.149444, -55.203889,
        134.624167, 57.504167,  174.777222, -175.2,     159.955556, 178.4419,   168.316667,
        -171.75,    -21.94,     179.2,      171.3667,   173,        166.920867, 55.45,
        158.158889};
    if (centerLats.size() != centerLons.size()) {
        std::cout << "BADDDD" << std::endl;
        return;
    }
    whichOne = centerLats.size() - whichOne;
    // bool centerCoordsAndRadiusInitialized = false;

    // std::string foundPercentageCirclesFilename = "foundPercentageCircles2020.txt";
    // std::ifstream foundPercentageCircles;
    // foundPercentageCircles.open(foundPercentageCirclesFilename);
    // std::string percentageCircleString;
    // while (getline(foundPercentageCircles, percentageCircleString)) {
    //     std::stringstream percentageCircleSS(percentageCircleString);
    //     std::string percentString;
    //     getline(percentageCircleSS, percentString, ' ');
    //     if (std::abs(std::stod(percentString) - percent) < FLOAT_THRESHOLD) {
    //         std::string radiusString;
    //         // TODO: See how to get bool from these getlines to see if file is formatted
    //         correctly getline(percentageCircleSS, radiusString, ' '); radius =
    //         std::stoi(radiusString); std::string lonString; getline(percentageCircleSS,
    //         lonString, ' '); cenLon = std::stod(lonString); std::string latString;
    //         getline(percentageCircleSS, latString, ' ');
    //         cenLat = std::stod(latString);
    //         centerCoordsAndRadiusInitialized = true;
    //     }
    // }
    // foundPercentageCircles.close();
    // if (!centerCoordsAndRadiusInitialized) {
    //     std::cout << "Desired percentage circle wasn't in " + foundPercentageCirclesFilename
    //               << std::endl;
    //     return;
    // }

    // TODO: Make a JSON-ize function for arrays, or find one
    // int output[NUM_ROWS][NUM_COLS];

    std::string landNBordersFileName = "landNBorders.txt";
    std::ifstream landNBorders;
    landNBorders.open(landNBordersFileName);
    std::ofstream colorsJSON;
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
            double distance = EquirectRasterData::distance(centerLats[whichOne],
                                                           centerLons[whichOne], lat(r), lon(c));
            char landOrBorder = rString.at(c);
            if (distance <= 60 && distance <= radius) {
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
                    char toSend = '1'; // Grey
                    for (int i = whichOne + 1; i < centerLats.size(); ++i) {
                        double d = EquirectRasterData::distance(centerLats[i], centerLons[i],
                                                                lat(r), lon(c));
                        if (d <= 60 && d <= radius) {
                            toSend = '6'; // Blue
                        } else if (d <= radius) {
                            if (toSend != '6') {
                                toSend = '7'; // Light blue
                            }
                        }
                    }

                    colorsJSON << toSend;
                }
            } else if (landOrBorder == '2') {
                if (distance <= radius && distance > radius - 60) {
                    // Light red
                    colorsJSON << '5';
                } else {
                    char toSend = '2'; // White
                    for (int i = whichOne + 1; i < centerLats.size(); ++i) {
                        double d = EquirectRasterData::distance(centerLats[i], centerLons[i],
                                                                lat(r), lon(c));
                        if (d <= radius && d > radius - 60) {
                            toSend = '8'; // Very light blue
                        }
                    }
                    colorsJSON << toSend;
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

int main() {
    // Federated States of Micronesia
    // oldMain(6.917222, 158.158889);
    // Seychelles
    for (int i = 177; i <= 197; ++i) {
        oldMain(i);
    }
    // Russia
    // oldMain(55.755833, 37.617222);
    // Nepal
    // oldMain(27.7172, 85.324);
}