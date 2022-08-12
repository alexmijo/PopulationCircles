// Needs GDAL library
// This is an early version of the program. It's still kinda spaghetti code and has code I copied
//  from random places on the internet unattributed (GeoTiff and most of distance()). If you're a
//  prospective employer, please look at BlackjackSim on my github instead for C++ stuff I've
//  written with better code style.

#include <iostream>
#include <gdal.h>
#include <string>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <map>
#include <iomanip>
#include <vector>
#include <array>

const double WORLD_POP_2015 = 7346242908.863955;
// See https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout#comment99267684_554134
const int DOUBLE_ROUND_TRIP_PRECISION =
    (std::numeric_limits<double>::digits10 == 15) ? 17 : std::numeric_limits<double>::digits10 + 3;
const int EQUATOR_LEN = 40075;

class EquirectRasterData {

private:

    const static int KERNEL_WIDTH = 4; // Num cols in each kernel (4 corners of a box)
    const static int SMALLEST_CIRCLES_VALUE_LENGTH = 3; // lon, lat, sum
    // key is radius, value is array holding lon, lat, sum, isGreaterThanOrEqualToResult
    std::map<int, double*> smallestCircleResults;
    std::string smallestCircleResultsFilename;
    int numRows, numCols;
    // First index is y (lat), second is x (lon). Matrix style indexing (top to bottom, left to
    //  right)
    double **sumTable;

    // Used for boundingBoxEdge
    enum direction {
        north,
        south
    };

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
            }
            else {
                edge += incOrDec;
            }
        }
        if (edge < 0) {
            edge = 0;
        }
        else if (edge > edgeOfMap) {
            edge = edgeOfMap;
        }
        return edge;
    }

    // For indexing into a kernel array
    inline int kernelIndex(const int row, const int col) {
        return KERNEL_WIDTH * row + col;
    }

    // Makes the kernel for a specific lattitude. Assigns the kernel's length to the reference
    //  parameter kernelLength.
    // TODO: Make the kernel an actual 2d array
    // TODO: See if this double counts part of 1 column when wrapping around all longitudes
    int* makeKernel(const int cenX, const int cenY, const double radius, int& kernelLength) {
        const double cenLon = lon(cenX);
        const double cenLat = lat(cenY);
        const int northEdge = boundingBoxEdge(cenX, cenY, radius, north);
        const int southEdge = boundingBoxEdge(cenX, cenY, radius, south);
        const int maxPossibleLength = southEdge - northEdge + 1;
        // A faster way of doing a 2D array of dimensions maxPossibleSize x KERNEL_WIDTH
        // Each row consists of: {westX, eastX, northY, southY} describing a summation table
        //  rectangle (so KERNEL_WIDTH must be 4) relative to cenX and cenY
        // Temp just cause we don't know length of real kernel yet
        int* tempKernel = new int[maxPossibleLength * KERNEL_WIDTH];

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
            }
            else {
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
                }
                else { // Find where this new rectangle ends
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
        int* kernel = new int[kernelLength * KERNEL_WIDTH];
        for (kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
            for (int kernelCol = 0; kernelCol < KERNEL_WIDTH; kernelCol++) {
                kernel[kernelIndex(kernelRow, kernelCol)] 
                    = tempKernel[kernelIndex(kernelRow, kernelCol)];
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
        return sumTable[south][east] - sumTable[north][east] - sumTable[south][west] 
               + sumTable[north][west];
    }

    // Returns the population of a circle with a radius of <radius> kilometers centered at the given
    //  latittude <cenLat> and longitude <cenLon>. TODO: Make this comment make more sense
    double popWithinKernel(const int cenX, const int cenY, int* kernel, const int kernelLength) {
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

    // Does the same thing as popWithinKernel, except faster for circles that cover more than half
    //  of the earth. It does this by summing up the area outside of the circle and then subtracting
    //  that from the total world population instead of summing up the area inside the circle.
    // TODO: Wait why did I think this would be any faster? It still has to go through the whole
    //  kernel, the only difference is that it might have less splitting along the antimeridian on
    //  average (although even that I'm not sure about). So I'm pretty sure this was sort of
    //  useless. I need to test to find out.
    // TODO: Instead of caring about if the circle has more than half the world's area, I should be
    //  doing something different if the circle has more than half the world's population.
    double popWithinKernelLargeCircle(const int cenX, const int cenY, int* kernel,
                                      const int kernelLength) {
        double popOutsideCircle = 0;
        // TODO: Fewer magic numbers please
        const int firstNorth = kernel[kernelIndex(0, 2)] + cenY;
        if (firstNorth > -1) {
            // Add all population north of the circle
            popOutsideCircle += popWithinRectangle(-1, numCols - 1, -1, firstNorth);
        }

        for (int kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
            // Sides of the original rectangle
            const int west = kernel[kernelIndex(kernelRow, 0)] + cenX;
            const int east = kernel[kernelIndex(kernelRow, 1)] + cenX;
            const int north = kernel[kernelIndex(kernelRow, 2)] + cenY;
            const int south = kernel[kernelIndex(kernelRow, 3)] + cenY;

            if (kernel[kernelIndex(kernelRow, 1)] == numCols / 2) {
                // Original rectangle encircles the entire latitude
                continue;
            } else if (west < -1) {
                // Originally needed to wrap around the antimeridian, creating two rectangles.
                //  Therefore we don't do that.
                popOutsideCircle += popWithinRectangle(east, numCols + west, north, south);
            } else if (east >= numCols) {
                // Originally needed to wrap around the antimeridian, creating two rectangles.
                //  Therefore we don't do that.
                popOutsideCircle += popWithinRectangle(east - numCols, west, north, south);
            } else {
                // Originally didn't need to wrap around the antimeridian, so now we do need to.
                popOutsideCircle += popWithinRectangle(-1, west, north, south);
                popOutsideCircle += popWithinRectangle(east, numCols - 1, north, south);
            }
        }

        const int lastSouth = kernel[kernelIndex(kernelLength - 1, 3)] + cenY;
        if (lastSouth < numRows - 1) {
            // Add all population south of the circle
            popOutsideCircle += popWithinRectangle(-1, numCols - 1, lastSouth, numRows - 1);
        }

        return WORLD_POP_2015 - popOutsideCircle;
    }

    // Get longitude of the center of the <x>th (0 indexed) column
    double lon(int x) {
        return (((double)x + 0.5) / (double)numCols) * 360.0 - 180.0;
    }

    // Get lattitude of the center of the <y>th (0 indexed) row
    double lat(int y) {
        return -((((double)y + 0.5) / (double)numRows) * 180.0 - 90.0);
    }
    
    // The inverse of the lon function
    int lonToX(double lon) {
        return ((lon + 180.0) / 360.0) * numCols - 0.5;
    }

    // The inverse of the lat function
    int latToY(double lat) {
        return ((-lat + 90.0) / 180.0) * numRows - 0.5;
    }

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
    // TODO: Make it just use the text file instead of sometimes creating it from the binary one
    EquirectRasterData(const std::string& sumTableFilename,
                       const std::string& smallestCircleResultsFilename) : 
                       smallestCircleResultsFilename(smallestCircleResultsFilename) {
        std::fstream sumTableFile;
        sumTableFile.open(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&numRows), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&numCols), sizeof(int));
        sumTable = new double*[numRows];
        for(int r = 0; r < numRows; r++) {
            sumTable[r] = new double[numCols];
            for (int c = 0; c < numCols; c++) {
                // TODO: See if it's faster to read in entire rows at a time
                sumTableFile.read(reinterpret_cast<char *>(&sumTable[r][c]), sizeof(double));
            }
        }
        sumTableFile.close();

        std::fstream smallestCircleResultsFile;
        smallestCircleResultsFile.open(smallestCircleResultsFilename);
        // See if file is empty. If so, make it
        std::streampos begin, end;
        begin = smallestCircleResultsFile.tellg();
        smallestCircleResultsFile.seekg(0, std::ios::end);
        end = smallestCircleResultsFile.tellg();
        smallestCircleResultsFile.seekg(0, std::ios::beg);
        if (begin - end == 0) {
            // TODO: See if any of this is even necessary
            std::cout << smallestCircleResultsFilename << " is empty or doesn't exist. Making it "
                "non-empty and existing." << std::endl;
            smallestCircleResultsFile.close();
            // Make the file
            // TODO: See if this is still needed now that it's text not binary
            smallestCircleResultsFile.open(smallestCircleResultsFilename, std::ios::out);
        } else {
            std::string resultString;
            while (getline(smallestCircleResultsFile, resultString)) {
                std::stringstream resultSS(resultString);
                std::string radiusString;
                getline(resultSS, radiusString, ' ');
                const int radius = std::stoi(radiusString);
                // TODO: Don't really need this as a variable
                // The + 1 is to hold whether or not it's a >= result
                // TODO: Hold the boolean some better way, probably by making this whole thing a
                //  class
                double *smallestCirclesValue = new double[SMALLEST_CIRCLES_VALUE_LENGTH + 1];
                // TODO: See if std::vector would make more sense here
                for (int j = 0; j < SMALLEST_CIRCLES_VALUE_LENGTH; j++) {
                    std::string doubleString;
                    getline(resultSS, doubleString, ' ');
                    smallestCirclesValue[j] = std::stod(doubleString);
                }
                std::string possibleGreaterThanOrEqualString;
                if (getline(resultSS, possibleGreaterThanOrEqualString)) {
                    if (possibleGreaterThanOrEqualString != ">=") {
                        // TODO: Print the actual name of the file
                        std::cout << "Incorrectly formatted smallestCircleResultsFile" << std::endl;
                        std::cerr << "Incorrectly formatted smallestCircleResultsFile" << std::endl;
                        throw 1776; // TODO: Make an actually descriptive exception here
                    }
                    smallestCirclesValue[SMALLEST_CIRCLES_VALUE_LENGTH] = true;
                } else {
                    smallestCirclesValue[SMALLEST_CIRCLES_VALUE_LENGTH] = false;
                }
                smallestCircleResults[radius] = smallestCirclesValue;
            }
        }
        smallestCircleResultsFile.close();
    }

    ~EquirectRasterData() {
        for (int r = 0; r < numRows; r++) {
            delete[] sumTable[r];
        }
        delete[] sumTable;
        for (std::map<int, double*>::iterator it = smallestCircleResults.begin();
             it != smallestCircleResults.end(); it++) {
            delete[] (it->second);
        }
    }

    // Returns distance in kilometers between two points on Earth's surface
    static double distance(const double& lat1, const double& lon1, const double& lat2, 
                           const double& lon2) {
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

        while (true) {

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
            lambda = L + (1 - C) * f * sin_alpha * (sigma + C * sin_sigma * (cos_2sigma_m + C 
                     * cos_sigma * (-1 + 2 * cos_2sigma_m * cos_2sigma_m)));
            diff = lambda - diff;
            if (std::fabs(diff) < 0.00001) { break; }
        }

        // U2
        double u_sq = cos_alpha_sq * (a * a - b * b) / (b * b);

        // Compute A, B
        double A = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)));

        double B = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (-128 + u_sq * (74 - 47 * u_sq))));

        // Sigma
        double cos_2sigma_m_sq = cos_2sigma_m * cos_2sigma_m;
        double delta_sigma = B * sin_sigma * (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 
                             * cos_2sigma_m_sq) - (B / 6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma 
                             * sin_sigma) * (-3 + 4 * cos_2sigma_m_sq)));

        // Distance
        double s = b * A * (sigma - delta_sigma);

        // Convert from meters to kilometers before returnings
        return s / 1000.;
    }

    int getNumRows() {
        return numRows;
    }

    int getNumCols() {
        return numCols;
    }

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
    // TODO: Update or remove this.
    double* largestSumCircleOfGivenRadius(const double radius, const double leftLon=-180, 
                                          const double rightLon=180, const double upLat=90, 
                                          const double downLat=-90, 
                                          const double initialLargestSum=0) {
        const int smallStep = 1; //1
        const int mediumStep = std::max((int)(radius / 128), 1); //4
        const int largeStep = std::max((int)(radius / 32), 1); //16
        const int xLStep = std::max((int)(radius / 8), 1); //64
        const int xXLStep = std::max((int)(radius / 2), 1); //256

        const double smallCutoff = 0.82; // 0.82
        const double mediumCutoff = 0.65; // 0.65
        const double largeCutoff = 0.47; // 0.47
        const double xLCutoff = 0.32; // 0.32
        
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
            int* kernel = makeKernel(1000, cenY, radius, kernelLength); // Initializes kernelLength

            for (int cenX = leftX; cenX <= rightX; cenX += step) {
                double popWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength);
                if (popWithinNKilometers >= largestSum) {
                    std::cout << "Sum within " << radius << " kilometers of ("
                        << lat(cenY) << ", " << lon(cenX) << "): "
                        << ((long long)popWithinNKilometers) << std::endl;
                    largestSumCenLon = lon(cenX);
                    largestSumCenLat = lat(cenY);
                    largestSum = popWithinNKilometers;
                }

                if (popWithinNKilometers > largestSum * smallCutoff) {
                    if (step > smallStep) {
                        cenX -= step;
                    }
                    step = smallStep;
                }
                else if (popWithinNKilometers > largestSum * mediumCutoff) {
                    if (step > mediumStep) {
                        cenX -= step;
                    }
                    step = mediumStep;
                }
                else if (popWithinNKilometers > largestSum * largeCutoff) {
                    if (step > largeStep) {
                        cenX -= step;
                    }
                    step = largeStep;
                }
                else if (popWithinNKilometers > largestSum * xLCutoff) {
                    if (step > xLStep) {
                        cenX -= step;
                    }
                    step = xLStep;
                }
                else {
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

    // First new approach to optimizing the function
    // TODO: Make this less spagettiish
    // Doesn't really restrict strictly to range yet (TODO)
    // TODO: Doesnt check easternmost or southernmost part of the world
    double* largestSumCircleOfGivenRadiusOpt1(const double radius, const double desiredSum,
                                              const double leftLon=-180, const double rightLon=180,
                                              const double upLat=90, const double downLat=-90) {
        std::cout << "Radius:" << radius << std::endl;
        // TODO: Don't look at centers outside these ranges
        // Turn lat and lon into indices in the pop data.
        const int leftX = lonToX(leftLon);
        const int rightX = lonToX(rightLon);
        const int upY = latToY(upLat);
        const int downY = latToY(downLat);

        int initialStep;
        double cutoff256;
        double cutoff64;
        double cutoff16;
        double cutoff4;
        if (radius >= 10400) {
            initialStep = 256;
            cutoff256 = 0.55;
            cutoff64 = 0.85;
            cutoff16 = 0.996;
            cutoff4 = 0.9985;
        } else if (radius >= 6700) {
            initialStep = 256;
            cutoff256 = 0.52;
            cutoff64 = 0.845;
            cutoff16 = 0.982;
            cutoff4 = 0.9921;
        } else if (radius >= 1150) {
            initialStep = 256;
            cutoff256 = 0.51;
            cutoff64 = 0.84;
            cutoff16 = 0.98;
            cutoff4 = 0.99;
        } else if (radius >= 300) {
            initialStep = 64;
            cutoff64 = 0.93;
            cutoff16 = 0.97;
            cutoff4 = 0.98;
        } else if (radius >= 100) {
            initialStep = 16;
            cutoff16 = 0.95;
            cutoff4 = 0.97;
        } else if (radius >= 20) {
            initialStep = 4;
            cutoff4 = 0.7;
        } else {
            initialStep = 1;
            cutoff4 = 0.999;
        }

        int step = initialStep; // Must be a power of 4, I think
        double cutoff = cutoff256;
        if (step <= 4) {
            cutoff = cutoff4;
        } else if (step <= 16) {
            cutoff = cutoff16;
        } else if (step <= 64) {
            cutoff = cutoff64;
        }
        // TODO: Make this one vector instead of three
        std::vector<int> topCenXs;
        std::vector<int> topCenYs;
        std::vector<double> topSums;
        double largestSum = 0;

        std::cout << "step: " << step << std::endl;
        largestSumCirclesOfGivenRadiusOpt1PixelBoundaries(radius, leftX, rightX, upY, downY, step,
                                                          topCenXs, topCenYs, topSums, largestSum,
                                                          cutoff, -100, -100);
        std::cout << "topCenXs.size(): " << topCenXs.size() << std::endl;
        cutUnderperformingCircles(topCenXs, topCenYs, topSums, largestSum, cutoff);

        while (step > 1) {
            step /= 4;
            std::cout << "step: " << step << std::endl;
            std::cout << "topCenXs.size(): " << topCenXs.size() << std::endl;
            std::cout << "largestSum: " << ((long)largestSum) << std::endl;
            
            if (largestSum > desiredSum) {
                for (int i = 0; i < topCenXs.size(); i++) {
                    if (topSums[i] == largestSum) {
                        double *returnValues = new double[4]; // TODO: No magic 4
                        returnValues[0] = lon(topCenXs[i]);
                        returnValues[1] = lat(topCenYs[i]);
                        returnValues[2] = largestSum;
                        returnValues[3] = true;
                        std::cout << "(SS)Sum within " << radius << " kilometers of ("
                            << returnValues[1] << ", " << returnValues[0] << "): "
                            << ((long)largestSum) << std::endl;
                        std::cout << ">= " << ((long)desiredSum) << ". Short circuiting."
                            << std::endl;
                        return returnValues;
                    }
                }
            }

            if (step <= 4) {
                cutoff = cutoff4;
            } else if (step <= 16) {
                cutoff = cutoff16;
            } else if (step <= 64) {
                cutoff = cutoff64;
            }
            const int numCirclesToSum = topCenXs.size();
            for (int i = 0; i < numCirclesToSum; i++) {
                int sectionLeftX = topCenXs[i] - step * 2;
                int sectionRightX = topCenXs[i] + step * 2 - 1;
                int sectionUpY = topCenYs[i] - step * 2;
                int sectionDownY = topCenYs[i] + step * 2 - 1;
                largestSumCirclesOfGivenRadiusOpt1PixelBoundaries(radius, sectionLeftX,
                                                                  sectionRightX, sectionUpY,
                                                                  sectionDownY, step, topCenXs,
                                                                  topCenYs, topSums, largestSum,
                                                                  cutoff, topCenXs[i], topCenYs[i]);
            }
            cutUnderperformingCircles(topCenXs, topCenYs, topSums, largestSum, cutoff);
        }
        std::cout << "topCenXs.size(): " << topCenXs.size() << std::endl;

        for (int i = 0; i < topCenXs.size(); i++) {
            if (topSums[i] == largestSum) {
                double *returnValues = new double[4]; // TODO: No magic 4
                returnValues[0] = lon(topCenXs[i]);
                returnValues[1] = lat(topCenYs[i]);
                returnValues[2] = largestSum;
                returnValues[3] = false;
                std::cout << "Sum within " << radius << " kilometers of (" << returnValues[1]
                    << ", " << returnValues[0] << "): " << ((long)largestSum) << std::endl;
                return returnValues;
            }
        }

        // TODO: Throw an exception here. Maybe nullptr return.
        std::cerr << "Never get here!!! No return value for Opt1" << std::endl;
        std::cout << "Never get here!!! No return value for Opt1" << std::endl;
        double *placeholder;
        return placeholder;
    }

    // TODO: Make this private, probably
    // Passed in ranges are inclusive.
    // Adds all centers that are at least 90% the sum of the largest center to topCenXs and topCenYs
    // Reassigns largestSum if necessary
    void largestSumCirclesOfGivenRadiusOpt1PixelBoundaries(const double radius, const int leftX,
                                                              const int rightX, const int upY,
                                                              const int downY, const int step,
                                                              std::vector<int>& topCenXs,
                                                              std::vector<int>& topCenYs,
                                                              std::vector<double>& topSums,
                                                              double& largestSum,
                                                              const double cutoff, const int skipX,
                                                              const int skipY) {
        for (int cenY = upY; cenY <= downY; cenY += step) {
            if (cenY < 0 || cenY >= numRows) {
                continue;
            }
            int kernelLength;
            int* kernel = makeKernel(1000, cenY, radius, kernelLength); // Initializes kernelLength

            for (int cenX = leftX; cenX <= rightX; cenX += step) {
                if (cenX < 0 || cenX >= numCols || (cenX == skipX && cenY == skipY)) {
                    continue;
                }
                double popWithinNKilometers;
                if (radius <= EQUATOR_LEN / 4) {
                    popWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength);
                } else {
                    popWithinNKilometers =
                        popWithinKernelLargeCircle(cenX, cenY, kernel, kernelLength);
                }
                // TODO: Make this logic simpler lol
                if (step != 1 && popWithinNKilometers > largestSum * cutoff) {
                    // std::cout << "Sum within " << radius << " kilometers of ("
                    //     << lat(cenY) << ", " << lon(cenX) << "): "
                    //     << ((long long)popWithinNKilometers) << std::endl;
                    topCenXs.push_back(cenX);
                    topCenYs.push_back(cenY);
                    topSums.push_back(popWithinNKilometers);
                    if (popWithinNKilometers > largestSum) {
                        largestSum = popWithinNKilometers;
                    }
                } else if (step == 1 && popWithinNKilometers >= largestSum) {
                    topCenXs.push_back(cenX);
                    topCenYs.push_back(cenY);
                    topSums.push_back(popWithinNKilometers);
                    largestSum = popWithinNKilometers;
                }
            }
            delete[] kernel;
        }
    }

    // Removes from topCenXs, topCenYs, and topSums any circles whose sum is less than largestSum *
    //  cutoff.
    // TODO: remove duplicates (actually, probably pass center into above function to prevent that)
    void cutUnderperformingCircles(std::vector<int>& topCenXs, std::vector<int>& topCenYs,
                                   std::vector<double>& topSums, const double largestSum,
                                   const double cutoff) {
        std::vector<int> indicesToErase;
        for (int i = topSums.size() - 1; i >= 0; i--) {
            if (topSums[i] < largestSum * cutoff) {
                indicesToErase.push_back(i);
            }
        }
        for (int indexToErase : indicesToErase) {
            topCenXs.erase(topCenXs.begin() + indexToErase);
            topCenYs.erase(topCenYs.begin() + indexToErase);
            topSums.erase(topSums.begin() + indexToErase);
        }
    }

    // Finds the circle of the given radius that minimizes the sum of the data inside it. Returns
    //  a pointer to an array containing the longitude [0] and lattitude [1] of the center of the
    //  circle and the sum of the data inside the circle [2].
    // Can optionally constrain the center of the circle to be within given ranges of latitudes and
    //  longitudes with leftLon, rightLon, upLat and downLat.
    // The parameter initialSmallestSum can speed up computation if it is passed in. Must be for sure known to
    //  be at least as large as what the smallestSum will end up being.
    // Radius given in kilometers.
    double* smallestSumCircleOfGivenRadius(const double radius, const double leftLon=-180, 
                                           const double rightLon=180, const double upLat=90, 
                                           const double downLat=-90, 
                                           const double initialSmallestSum=1000000000000000) {
        const int smallStep = 1; //1
        const int mediumStep = std::max((int)(radius / 128), 1); //4
        const int largeStep = std::max((int)(radius / 32), 1); //16
        const int xLStep = std::max((int)(radius / 8), 1); //64
        const int xXLStep = std::max((int)(radius / 2), 1); //256
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
            int* kernel = makeKernel(1000, cenY, radius, kernelLength); // Initializes kernelLength

            for (int cenX = leftX; cenX <= rightX; cenX += step) {
                double sumWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength);
                if (sumWithinNKilometers <= smallestSum) {
                    std::cout << "Sum within " << radius << " kilometers of ("
                        << lat(cenY) << ", " << lon(cenX) << "): "
                        << ((long long)sumWithinNKilometers) << std::endl;
                    smallestSumCenLon = lon(cenX);
                    smallestSumCenLat = lat(cenY);
                    smallestSum = sumWithinNKilometers;
                }

                if (sumWithinNKilometers < smallestSum * 1.2) {
                    if (step > smallStep) {
                        cenX -= step;
                    }
                    step = smallStep;
                }
                else if (sumWithinNKilometers < smallestSum * 1.4) {
                    if (step > mediumStep) {
                        cenX -= step;
                    }
                    step = mediumStep;
                }
                else if (sumWithinNKilometers < smallestSum * 1.6) {
                    if (step > largeStep) {
                        cenX -= step;
                    }
                    step = largeStep;
                }
                else if (sumWithinNKilometers < smallestSum * 1.8) {
                    if (step > xLStep) {
                        cenX -= step;
                    }
                    step = xLStep;
                }
                else {
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
    // TODO: Figure out why largestSum seems to stay the same after searching at a step that's 1/4
    //  the size something like 1/4 of the time rather than 1/16 of the time like expected.
    double* smallestCircleWithGivenSum(const double sum, const double leftLon=-180, 
                                       const double rightLon=180, const double upLat=90, 
                                       const double downLat=-90) {
        // Radii will be ints cause only interested in getting to the nearest kilometer
        int upperBound = EQUATOR_LEN / 2;
        int lowerBound = 0;
        double *returnValues;
        bool returnValuesAssigned = false; // Pointer not yet pointing to an array
        double initialLargestSum = 0; // To be passed into largestSumCircleOfGivenRadius()
        bool softUpperBound = false;

        // Tighten bounds as much as possible using previous results
        for (std::map<int, double*>::iterator it = smallestCircleResults.begin();
             it != smallestCircleResults.end(); it++) {
            const bool isAGreaterThanOrEqualToResult = it->second[3];
            if ((it->second)[2] < sum && it->first >= lowerBound
                && !isAGreaterThanOrEqualToResult) {
                lowerBound = it->first;
                initialLargestSum = (it->second)[2];
            } else if ((it->second)[2] >= sum && it->first <= upperBound) {
                if (it->first == upperBound && isAGreaterThanOrEqualToResult) {
                    // Don't replace hard upper bound with a soft one
                    continue;
                }
                upperBound = it->first;
                returnValues = new double[4];
                returnValues[0] = (it->second)[0];
                returnValues[1] = (it->second)[1];
                returnValues[2] = (it->second)[2];
                returnValues[3] = it->first;
                returnValuesAssigned = true;
                softUpperBound = isAGreaterThanOrEqualToResult;
                break;
            }
        }

        int radius = lowerBound + (upperBound - lowerBound) / 2; // Start of binary search
        while (upperBound - lowerBound > 1 || softUpperBound) {
            double narrowLeftLon = leftLon;
            double narrowRightLon = rightLon;
            double narrowUpLat = upLat;
            double narrowDownLat = downLat;
            // May be able to use knowledge from past map results (the ones on reddit) to narrow the
            //  search area
            if (radius < 117) {
                // Uncharted territory, don't narrow search area
            } else if (radius <= 3500) {
                narrowLeftLon = 15;
                narrowRightLon = 156;
                narrowUpLat = 67;
                narrowDownLat = -10;
            } else if (radius < 16000) {
                narrowDownLat = -50;
            }

            double *largestSumCircle;
            if (upperBound - lowerBound <= 3) {
                if (upperBound - lowerBound == 1) {
                    // Need to make upperBound not soft
                    radius = upperBound;
                }
                // Huge desired sum since we don't want it to short circuit
                largestSumCircle = largestSumCircleOfGivenRadiusOpt1(radius, 10000000000000,
                                                                     narrowLeftLon, narrowRightLon,
                                                                     narrowUpLat, narrowDownLat);
            } else {
                largestSumCircle = largestSumCircleOfGivenRadiusOpt1(radius, sum, narrowLeftLon,
                                                                     narrowRightLon, narrowUpLat,
                                                                     narrowDownLat);
            }
            if (largestSumCircle[2] >= sum) {
                if (largestSumCircle[3]) {
                    softUpperBound = true;
                } else {
                    softUpperBound = false;
                }
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

            const bool isNewGreaterThanOrEqualResult = largestSumCircle[3];
            // Add result to smallestCircleResults and its file
            std::map<int, double*>::iterator it = smallestCircleResults.find(radius);
            if (it != smallestCircleResults.end() && !it->second[3]) {
                // TODO: Figure out a better way to do these sort of things (probably throw an
                //  exception)
                // TODO: See if there's a way to send a string to both streams
                std::cout << "smallestCircleWithGivenSum had an error involving "
                    "smallestCircleResults" << std::endl;
                std::cerr << "smallestCircleWithGivenSum had an error involving "
                    "smallestCircleResults" << std::endl; // Above line might get drowned out
            } else {
                std::fstream smallestCircleResultsFile;
                smallestCircleResultsFile.open(smallestCircleResultsFilename);
                smallestCircleResultsFile.seekg(0, std::ios::end);
                smallestCircleResultsFile << radius
                                          << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION);
                // TODO: See if I can just write the entire array at once
                // Note the lack of a + 1 here
                for (int j = 0; j < SMALLEST_CIRCLES_VALUE_LENGTH; j++) {
                    smallestCircleResultsFile << " " << largestSumCircle[j];
                }
                if (isNewGreaterThanOrEqualResult) {
                    smallestCircleResultsFile << " >=";
                }
                smallestCircleResultsFile << std::setprecision(6) << "\n";
                smallestCircleResultsFile.close();
                // TODO: Instead of this workaround, just don't delete[] largestSumCircle in
                //  this case
                // TODO: Get rid of the magic + 1, here to hold >= boolean
                smallestCircleResults[radius] = new double[SMALLEST_CIRCLES_VALUE_LENGTH + 1];
                smallestCircleResults[radius][0] = largestSumCircle[0];
                smallestCircleResults[radius][1] = largestSumCircle[1];
                smallestCircleResults[radius][2] = largestSumCircle[2];
                smallestCircleResults[radius][3] = isNewGreaterThanOrEqualResult;
            }

            delete[] largestSumCircle;
// TODO: This was done incorrectly. ------------------------------------------------------------------------------------------------------------
            const double popAbove = smallestCircleResults[upperBound][2] - smallestCircleResults[radius][2];
            const double popBelow = smallestCircleResults[radius][2] - smallestCircleResults[lowerBound][2];
            const double sumMinusPop = sum - smallestCircleResults[radius][2];
            const double popMinusSum = smallestCircleResults[radius][2] - sum;
            const int range = upperBound - lowerBound;
            const std::array<int, 20> nonBinarySearches = {3, 4, 5, 8, 12, 16, 20, 24, 28, 32, 38, 44, 54, 64, 73, 100, 128, 256, 512, 1024};
            const std::array<int, 20> nonBinaryCutoffs = {1/3, 1/4, 1/5, 1/8, 1/12, 1/16, 1/20, 1/24, 1/28, 1/32, 1/38, 1/44, 1/54, 1/64, 1/73, 1/100, 1/128, 1/256, 1/512, 1/1024};
            //--------------------------------------------------------------------------------------
            // 0
            /*/ TODO
            if (sumMinusPop < nonBinaryCutoffs[0] * popMinusSum && range > nonBinarySearches[0]) {
                radius = upperBound - range / nonBinarySearches[0];
            } else if (popMinusSum < nonBinaryCutoffs[0] * sumMinusPop && range > nonBinarySearches[0]) {
                radius = lowerBound + range / nonBinarySearches[0];
            //--------------------------------------------------------------------------------------
            // 1
            } else if (sumMinusPop < nonBinaryCutoffs[1] * popMinusSum && range > nonBinarySearches[1]) {
                radius = upperBound - range / nonBinarySearches[1];
            } else if (popMinusSum < nonBinaryCutoffs[1] * sumMinusPop && range > nonBinarySearches[1]) {
                radius = lowerBound + range / nonBinarySearches[1];
            //--------------------------------------------------------------------------------------
            // 2
            } else if (sumMinusPop < nonBinaryCutoffs[2] * popMinusSum && range > nonBinarySearches[2]) {
                radius = upperBound - range / nonBinarySearches[2];
            } else if (popMinusSum < nonBinaryCutoffs[2] * sumMinusPop && range > nonBinarySearches[2]) {
                radius = lowerBound + range / nonBinarySearches[2];
            //--------------------------------------------------------------------------------------
            // 3
            // TODO
            } else {
                radius = lowerBound + (upperBound - lowerBound) / 2; // Binary search
            }*/
            radius = lowerBound + (upperBound - lowerBound) / 2; // Binary search
        }
        if (!returnValuesAssigned) {
            // TODO: Figure out a better way to do these sort of things (probably throw an
            //  exception)
            std::cout << "smallestCircleWithGivenSum wasn't able to find a circle with a large "
                "enough sum. Either the sum (" << sum << ") was too large, or a bug occurred." 
                << std::endl;
        }
        return returnValues;
    }
};

// See what I can get away with
void testCircleSkipping() {
    double radius = 11029;

    std::cout << "Loading population summation table." << std::endl;
    std::string sumTableFilename = "popSumTable.bin";
    std::string smallestCircleResultsFilename = "popSmallestCircleResults.txt";
    EquirectRasterData data(sumTableFilename, smallestCircleResultsFilename);
    std::cout << "Loaded population summation table." << std::endl;

    // Huge desiredSum so that it doesn't short circuit
    double *smallestCircle = data.largestSumCircleOfGivenRadiusOpt1(radius, 10000000000000);
    std::cout << "Population within " << radius << " km of (" << smallestCircle[1] 
        << ", " << smallestCircle[0] << "): " << ((long long)(smallestCircle[2])) << std::endl;
    delete[] smallestCircle;
}

// What was in the main function before.
void normalMain() {
    //-------------------------------Parameters---------------------------------------------
    // TODO: Make it so that the program works even if populationMode is false
    const bool populationMode = true; // TODO: make this an enum. false means GDP PPP mode
    // TODO: Actually use smallestPopMode
    // const bool smallestPopMode = false;
    //-------------------------------Parameters-end-----------------------------------------

    std::cout << "Loading population summation table." << std::endl;
    std::string sumTableFilename = "popSumTable.bin";
    std::string smallestCircleResultsFilename = "popSmallestCircleResults.txt";
    EquirectRasterData data(sumTableFilename, smallestCircleResultsFilename);

    std::cout << "Loaded " << (populationMode ? "population" : "GDP PPP") << " summation table." 
        << std::endl;

    // Used by imageManipulation.java so that it knows what text to add, as well as by
    //  percentageCirclesMapMakerTextInJSONOut.cpp so it knows where and how large to draw the
    //  circles
    std::string percentCirclesFilename = "foundPercentageCircles.txt";
    std::ofstream percentCirclesFile;
    // TODO: See if this should be opened and reclosed in each iteration to make sure recent results
    //  are written before any keyboard interrupt.
    percentCirclesFile.open(percentCirclesFilename);
    for (int percent = 1; percent <= 100; percent++) {
        const double desiredPopulation = (WORLD_POP_2015 / 100.0) * percent;
        std::cout << std::endl << "Now finding smallest possible circle with " << percent
                  << "\% of the world's population (" << ((long)desiredPopulation) << " people)"
                  << std::endl;
        double *smallestCircle = data.smallestCircleWithGivenSum(desiredPopulation);
        // TODO: Fewer magic numbers for indexing into smallestCircle (probably make a circle class)
        const double longitude = smallestCircle[0];
        const double lattitude = smallestCircle[1];
        const double population = smallestCircle[2];
        const int radius = smallestCircle[3];
        delete[] smallestCircle;

        std::cout << "Smallest possible circle with " << percent << "\% of the world's population ("
                  << ((long)desiredPopulation) << " people):" << std::endl;
        std::cout << "Population within " << radius << " km of (" << lattitude << ", " << longitude
                  << "): " << ((long)population) << std::endl; 
        // Used to be for entering into the python map making code, now just to match the format in
        //  the google doc
        std::cout << percent << ": (" << radius << ", (" << std::setprecision(8) << lattitude
                  << ", " << longitude << "))" << std::setprecision(6) << std::endl;

        // TODO: Remove the magic number 6
        percentCirclesFile << ((int)percent) << " " << radius << " "
                           << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION) << longitude << " "
                           << lattitude << " " << population << std::setprecision(6) << std::endl;
    }
    percentCirclesFile.close();
}

int main() {
    normalMain();
}