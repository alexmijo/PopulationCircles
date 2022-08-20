// Uses GDAL library
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
#include <stdexcept>

const double WORLD_POP_2015 = 7346242908.863955;
// See https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout#comment99267684_554134
const int DOUBLE_ROUND_TRIP_PRECISION =
    (std::numeric_limits<double>::digits10 == 15) ? 17 : std::numeric_limits<double>::digits10 + 3;
const int EQUATOR_LEN = 40075;

struct SmallestCircleResult {
    // The lattitude and longitude of the center of the circle.
    double lat, lon;
    // The radius of the circle in kilometers.
    double radius;
    // The population within the circle.
    double pop;
};

// Represents a previously found population circle, except instead of the radius it has whether or
//  the result was a >= result cause the algorithm short circuited.
struct SmallestCircleResultMaybeShortCircuit {
    // The lattitude and longitude of the center of the circle.
    // TODO: Now that the order is reversed, reverse it in the results file too.
    double lat, lon;
    // The population within the circle.
    double pop;
    // True iff this circle was found by the algorithm running to completion, false if the
    //  algorithm "short circuited" and was a >= result (i.e., there may be other circles with a
    //  bigger population but the same radius).
    bool shortCircuited;
};

class EquirectRasterData {
private:
    const static int KERNEL_WIDTH = 4; // Num cols in each kernel (4 corners of a box)
    // key is radius
    std::map<int, SmallestCircleResultMaybeShortCircuit> smallestCircleResults;
    std::string smallestCircleResultsFilename;
    int numRows, numCols;
    // First index is y (lat), second is x (lon). Matrix style indexing (top to bottom, left to
    //  right)
    // TODO: See if I can or should make this an array of arrays (library arrays). Can't be a 2D C
    //  array cause it's too large.
    // TODO: See if I should actually make this a 1D array (both for loading speed and using speed)
    double **sumTable;

    // Used for boundingBoxEdge
    enum direction {
        north,
        south
    };

    // Returns the northernmost or southernmost row containing any pixel within <radius> km of the
    //  given pixel (specified by <x> and <y>).
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
            throw std::invalid_argument("BoundingBoxEdge called with bad direction.");
        }

        while (edge >= 0 && edge <= edgeOfMap) {
            double currLat, currLon;
            currLat = lat(edge);
            // Since this function only works for north and south, longitude never changes
            currLon = cenLon;
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
    // TODO: Make the kernel a vector instead of a dynamically allocated c array
    // TODO: Make the kernel an actual 2d array (maybe)
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
                    throw std::logic_error("Kernel making failed. Should never get here.");
                }
                horizontalOffset--;
            }
            else {
                // See how much farther out we can go at this y level and still be in the circle
                while (distance(currLat, currLon, cenLat, cenLon) <= radius) {
                    horizontalOffset++;
                    if (horizontalOffset > numCols / 2) {
                        // This rectangle wraps around the world
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
                if (y == southEdge) {
                    // If we just started a new rectangle at southEdge, it must end there too
                    tempKernel[kernelIndex(kernelRow, 3)] = y - cenY;
                    break;
                }
                else {
                    // Find where this new rectangle ends
                    while (true) {
                        y++;
                        if (y > southEdge) {
                            // Rectangle can't extend below south edge
                            tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1;
                            break;
                        }

                        // Check if the circle has widened
                        if (horizontalOffset < numCols / 2) {
                            // Rectangles that wrap around the whole world can't widen any more
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
                // West side of rectangle on the antimeridian and north side of rectangle on the
                //  north pole
                return sumTable[south][east];
            }
            // West side of rectangle on the antimeridian
            return sumTable[south][east] - sumTable[north][east];
        } else if (north == -1) {
            // North side of rectangle on the north pole
            return sumTable[south][east] - sumTable[south][west];
        }
        return sumTable[south][east] - sumTable[north][east] - sumTable[south][west] +
               sumTable[north][west];
    }

    // Returns the population contained within a circle who's boundary is specified by <kernel> and
    //  which is centered at the (<cenX>, <cenY>) pixel.
    double popWithinKernel(const int cenX, const int cenY, int* kernel, const int kernelLength) {
        double totalPop = 0;
        for (int kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
            // Sides of the rectangle
            const int west = kernel[kernelIndex(kernelRow, 0)] + cenX;
            const int east = kernel[kernelIndex(kernelRow, 1)] + cenX;
            const int north = kernel[kernelIndex(kernelRow, 2)] + cenY;
            const int south = kernel[kernelIndex(kernelRow, 3)] + cenY;
            if (kernel[kernelIndex(kernelRow, 1)] == numCols / 2) {
                // This rectangle encircles the entire latitude
                totalPop += popWithinRectangle(-1, numCols - 1, north, south);
            } else if (west < -1) {
                // Need to wrap around the antimeridian, creating two rectangles
                totalPop += popWithinRectangle(-1, east, north, south);
                totalPop += popWithinRectangle(numCols + west, numCols - 1, north, south);
            } else if (east >= numCols) {
                // Need to wrap around the antimeridian, creating two rectangles
                totalPop += popWithinRectangle(west, numCols - 1, north, south);
                totalPop += popWithinRectangle(-1, east - numCols, north, south);
            } else {
                totalPop += popWithinRectangle(west, east, north, south);
            }
        }
        return totalPop;
    }

    // Get longitude of the center of the <x>th (0 indexed) column
    double lon(int x) {
        return ((x + 0.5) / numCols) * 360.0 - 180.0;
    }

    // Get lattitude of the center of the <y>th (0 indexed) row
    double lat(int y) {
        return -(((y + 0.5) / numRows) * 180.0 - 90.0);
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
    // TODO: Complete spec
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
            // Reads in an entire row at once
            sumTableFile.read(reinterpret_cast<char *>(sumTable[r]), sizeof(double) * numCols);
        }
        sumTableFile.close();
        std::fstream smallestCircleResultsFile;
        smallestCircleResultsFile.open(smallestCircleResultsFilename);
        // See if file is empty. If so, make it
        // TODO: Have this be a function
        std::streampos begin = smallestCircleResultsFile.tellg();
        smallestCircleResultsFile.seekg(0, std::ios::end);
        std::streampos end = smallestCircleResultsFile.tellg();
        smallestCircleResultsFile.seekg(0, std::ios::beg);
        if (begin - end == 0) {
            std::cout << smallestCircleResultsFilename << " is empty or doesn't exist. Making it "
                "non-empty and existing." << std::endl;
            smallestCircleResultsFile.close();
            // Make the file
            // TODO: See if this is still needed now that it's text not binary
            smallestCircleResultsFile.open(smallestCircleResultsFilename, std::ios::out);
        }
        std::string resultString;
        while (getline(smallestCircleResultsFile, resultString)) {
            std::stringstream resultSS(resultString);
            std::string radiusString;
            getline(resultSS, radiusString, ' ');
            const int radius = std::stoi(radiusString);
            // TODO: Don't really need this as a variable
            // The + 1 is to hold whether or not it's a >= result
            // TODO: Hold the boolean some better way, probably by making this whole thing a
            //  class;
            SmallestCircleResultMaybeShortCircuit result;
            std::string dummyString;
            getline(resultSS, dummyString, ' ');
            result.lon = std::stod(dummyString);
            getline(resultSS, dummyString, ' ');
            result.lat = std::stod(dummyString);
            getline(resultSS, dummyString, ' ');
            result.pop = std::stod(dummyString);
            result.shortCircuited = false;
            if (getline(resultSS, dummyString)) {
                if (dummyString != ">=") {
                    throw std::runtime_error("Incorrectly formatted smallest circle results file: "
                                             + smallestCircleResultsFilename);
                }
                result.shortCircuited = true;
            }
            smallestCircleResults[radius] = result;
        }
        smallestCircleResultsFile.close();
    }

    ~EquirectRasterData() {
        for (int r = 0; r < numRows; r++) {
            delete[] sumTable[r];
        }
        delete[] sumTable;
    }

    // Returns distance in kilometers between two points on Earth's surface
    // TODO: Find source and give credit
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
        double delta_sigma = B * sin_sigma * (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 
                             * cos_2sigma_m_sq) - (B / 6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma 
                             * sin_sigma) * (-3 + 4 * cos_2sigma_m_sq)));

        // Distance
        double s = b * A * (sigma - delta_sigma);

        // Convert from meters to kilometers before returning
        return s / 1000.;
    }

    int getNumRows() {
        return numRows;
    }

    int getNumCols() {
        return numCols;
    }

    // TODO: Fix the spec below, which is outdated 
    // Finds the circle of the given radius containing maximal population. Returns
    //  a pointer to an array containing the longitude [0] and lattitude [1] of the center of the
    //  circle and the population inside the circle [2].
    // Can optionally constrain the center of the circle to be within given ranges of latitudes and
    //  longitudes with leftLon, rightLon, upLat and downLat.
    // The parameter initialLargestSum can speed up computation if it is passed in. Must be for sure
    //  known to be at least as small as what the largestPop will end up being.
    // Radius given in kilometers.
    // TODO: Make this less spagettiish
    // Doesn't really restrict strictly to range yet (TODO)
    // TODO: Doesnt check easternmost or southernmost part of the world
    // TODO: Keep this public but have a more normal return value, and then make a private version
    //  that short circuits. Also then make SmallestCircleResultMaybeShortCircuit private. This will
    //  also take care of the spec and desiredPop not making sense here.
    SmallestCircleResultMaybeShortCircuit mostPopulousCircleOfGivenRadius(
        const double radius, const double desiredPop, const double leftLon=-180,
        const double rightLon=180, const double upLat=90, const double downLat=-90) {
        std::cout << "Radius: " << radius << std::endl;
        // TODO: Don't look at centers outside these ranges.
        // Turn lat and lon into indices in the pop data.
        const int leftX = lonToX(leftLon);
        const int rightX = lonToX(rightLon);
        const int upY = latToY(upLat);
        const int downY = latToY(downLat);

        // TODO: Make this nicer
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
        std::vector<double> topPops;
        double largestPop = 0;

        std::cout << "step: " << step << std::endl;
        // -100 for skipY and skipX since we don't want to skip any pixels
        mostPopulousCirclesOfGivenRadiusPixelBoundaries(
            radius, leftX, rightX, upY, downY, step, topCenXs, topCenYs, topPops, largestPop,
            cutoff, -100, -100);
        std::cout << "topCenXs.size(): " << topCenXs.size() << std::endl;
        cutUnderperformingCircles(topCenXs, topCenYs, topPops, largestPop, cutoff);
        SmallestCircleResultMaybeShortCircuit result;
        while (step > 1) {
            step /= 4;
            std::cout << "step: " << step << std::endl;
            std::cout << "topCenXs.size(): " << topCenXs.size() << std::endl;
            std::cout << "largestPop: " << ((long)largestPop) << std::endl;
            
            if (largestPop > desiredPop) {
                for (int i = 0; i < topCenXs.size(); i++) {
                    if (topPops[i] == largestPop) {
                        result.lat = lat(topCenYs[i]);
                        result.lon = lon(topCenXs[i]);
                        result.pop = largestPop;
                        result.shortCircuited = true;
                        std::cout << "(SS)Population within " << radius << " kilometers of ("
                            << result.lat << ", " << result.lon << "): "
                            << ((long)result.pop) << std::endl;
                        std::cout << ">= " << ((long)desiredPop) << ". Short circuiting."
                            << std::endl;
                        return result;
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
                mostPopulousCirclesOfGivenRadiusPixelBoundaries(
                    radius, sectionLeftX, sectionRightX, sectionUpY, sectionDownY, step, topCenXs,
                    topCenYs, topPops, largestPop, cutoff, topCenXs[i], topCenYs[i]);
            }
            cutUnderperformingCircles(topCenXs, topCenYs, topPops, largestPop, cutoff);
        }
        std::cout << "topCenXs.size(): " << topCenXs.size() << std::endl;
        // Now find and return the most populous circle in the list (all three vectors).
        for (int i = 0; i < topCenXs.size(); i++) {
            if (topPops[i] == largestPop) {
                result.lat = lat(topCenYs[i]);
                result.lon = lon(topCenXs[i]);
                result.pop = largestPop;
                result.shortCircuited = false;
                std::cout << "Population within " << radius << " kilometers of (" << result.lat
                    << ", " << result.lon << "): " << ((long)result.pop) << std::endl;
                return result;
            }
        }
        // largestPop should be in topPops
        throw std::logic_error("Should never get here");
    }

    // TODO: Make this private, probably
    // Passed in ranges are inclusive.
    // Adds all centers that are at least 90% the pop of the largest center to topCenXs and topCenYs
    // Reassigns largestPop if necessary
    void mostPopulousCirclesOfGivenRadiusPixelBoundaries(
        const double radius, const int leftX, const int rightX, const int upY, const int downY, 
        const int step, std::vector<int>& topCenXs, std::vector<int>& topCenYs,
        std::vector<double>& topPops, double& largestPop, const double cutoff, const int skipX,
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
                double popWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength);
                // TODO: Make this logic simpler lol
                if (step != 1 && popWithinNKilometers > largestPop * cutoff) {
                    topCenXs.push_back(cenX);
                    topCenYs.push_back(cenY);
                    topPops.push_back(popWithinNKilometers);
                    if (popWithinNKilometers > largestPop) {
                        largestPop = popWithinNKilometers;
                    }
                } else if (step == 1 && popWithinNKilometers >= largestPop) {
                    topCenXs.push_back(cenX);
                    topCenYs.push_back(cenY);
                    topPops.push_back(popWithinNKilometers);
                    largestPop = popWithinNKilometers;
                }
            }
            delete[] kernel;
        }
    }

    // Removes from topCenXs, topCenYs, and topPops any circles whose pop is less than largestPop *
    //  cutoff.
    void cutUnderperformingCircles(std::vector<int>& topCenXs, std::vector<int>& topCenYs,
                                   std::vector<double>& topPops, const double largestPop,
                                   const double cutoff) {
        std::vector<int> indicesToErase;
        for (int i = topPops.size() - 1; i >= 0; i--) {
            if (topPops[i] < largestPop * cutoff) {
                indicesToErase.push_back(i);
            }
        }
        for (int indexToErase : indicesToErase) {
            topCenXs.erase(topCenXs.begin() + indexToErase);
            topCenYs.erase(topCenYs.begin() + indexToErase);
            topPops.erase(topPops.begin() + indexToErase);
        }
    }

    // TODO: Make a spec comment for this
    // Kind of dangerous to pass in the last 4 args to this function, since the
    //  smallestCircleResultsFile might get spurious data added to it then
    // TODO: Figure out why largestPop seems to stay the same after searching at a step that's 1/4
    //  the size something like 1/4 of the time rather than 1/16 of the time like expected.
    SmallestCircleResult smallestCircleWithGivenPopulation(
        const double pop, const double leftLon=-180, const double rightLon=180,
        const double upLat=90, const double downLat=-90) {
        // Radii will be ints cause only interested in getting to the nearest kilometer
        int upperBound = EQUATOR_LEN / 2;
        int lowerBound = 0;
        SmallestCircleResult result;
        bool foundSuitableCircle = false; // Haven't yet found a circle with sufficiently large pop
        // TODO: Is this still even used?
        double initialLargestPop = 0; // To be passed into mostPopulousCircleOfGivenRadius()
        // A "soft" upper bound is one which comes from a short circuited previous result
        bool softUpperBound = false;
        // Tighten bounds as much as possible using previous results
        for (auto it = smallestCircleResults.begin(); it != smallestCircleResults.end(); it++) {
            if (it->second.pop < pop && it->first >= lowerBound && !it->second.shortCircuited) {
                lowerBound = it->first;
                initialLargestPop = it->second.pop;
            } else if (it->second.pop >= pop && it->first <= upperBound) {
                if (it->first == upperBound && it->second.shortCircuited) {
                    // Don't replace hard upper bound with an equal (due to smallestCircleResults
                    //  being sorted by radius) soft one
                    continue;
                }
                upperBound = it->first;
                result.lon = it->second.lon;
                result.lat = it->second.lat;
                result.radius = it->first;
                result.pop = it->second.pop;
                softUpperBound = it->second.shortCircuited;
                // Should only return this if it's not a short circuited result
                foundSuitableCircle = !softUpperBound;
                // Can break here since smallestCircleResults is sorted by radius
                break;
            }
        }

        int radius = lowerBound + (upperBound - lowerBound) / 2; // Start of binary search
        while (upperBound - lowerBound > 1 || softUpperBound) {
            // TODO: Could make a boundaries class and a function taking radius and returning
            //  boundaries
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
            SmallestCircleResultMaybeShortCircuit largestSumCircle;
            if (upperBound - lowerBound <= 3) {
                if (upperBound - lowerBound == 1) {
                    // Need to make upperBound not soft
                    radius = upperBound;
                }
                // Huge desired pop since we don't want it to short circuit
                largestSumCircle = mostPopulousCircleOfGivenRadius(radius, 10000000000000,
                                                                   narrowLeftLon, narrowRightLon,
                                                                   narrowUpLat, narrowDownLat);
            } else {
                largestSumCircle = mostPopulousCircleOfGivenRadius(radius, pop, narrowLeftLon,
                                                                   narrowRightLon, narrowUpLat,
                                                                   narrowDownLat);
            }
            if (largestSumCircle.pop >= pop) {
                softUpperBound = largestSumCircle.shortCircuited;
                upperBound = radius;
                result.lon = largestSumCircle.lon;
                result.lat = largestSumCircle.lat;
                result.radius = radius;
                result.pop = largestSumCircle.pop;
                foundSuitableCircle = true;
            } else {
                lowerBound = radius;
                initialLargestPop = largestSumCircle.pop;
            }
            // Add result to smallestCircleResults and its file
            // TODO: Put this in one or more functions
            auto it = smallestCircleResults.find(radius);
            if (it != smallestCircleResults.end() && !it->second.shortCircuited) {
                throw std::logic_error("Redid a previous result, which should never happen");
            } else {
                // TODO: Put this in a function
                std::fstream smallestCircleResultsFile;
                smallestCircleResultsFile.open(smallestCircleResultsFilename);
                smallestCircleResultsFile.seekg(0, std::ios::end);
                smallestCircleResultsFile << radius
                                          << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION);
                smallestCircleResultsFile << " " << largestSumCircle.lon << " "
                                          << largestSumCircle.lat << " " << largestSumCircle.pop
                                          << (largestSumCircle.shortCircuited ? " >=\n" : "\n");
                smallestCircleResultsFile.close();
                smallestCircleResults[radius] = largestSumCircle;
            }
            radius = lowerBound + (upperBound - lowerBound) / 2; // Binary search
        }
        if (!foundSuitableCircle) {
            // TODO: Fail fast
            if (pop > WORLD_POP_2015) {
                throw std::invalid_argument("Desired pop is larger than the world population, and a"
                    " suitable circle wasn't found");
            }
            throw std::logic_error("Should never get here");
        }
        return result;
    }
};

// See what I can get away with
void testCircleSkipping() {
    double radius = 11029;

    std::cout << "Loading population summation table." << std::endl;
    std::string sumTableFilename = "popSumTable.bin";
    std::string smallestCircleResultsFilename = "popSmallestCircleResults.txt";
    EquirectRasterData popData(sumTableFilename, smallestCircleResultsFilename);
    std::cout << "Loaded population summation table." << std::endl;

    // Huge desiredPop so that it doesn't short circuit
    SmallestCircleResultMaybeShortCircuit smallestCircle =
        popData.mostPopulousCircleOfGivenRadius(radius, 10000000000000);
    std::cout << "Population within " << radius << " km of (" << smallestCircle.lat << ", "
              << smallestCircle.lon << "): " << (long long)smallestCircle.pop << std::endl;
}

// Used for finding the smallest circles containing 1-100% of the world's population.
void findPercentCircles() {
    std::cout << "Loading population summation table." << std::endl;
    std::string sumTableFilename = "popSumTable.bin";
    std::string smallestCircleResultsFilename = "popSmallestCircleResults.txt";
    EquirectRasterData popData(sumTableFilename, smallestCircleResultsFilename);
    std::cout << "Loaded population summation table." << std::endl;

    // Used by imageManipulation.java so that it knows what text to add, as well as by
    //  percentageCirclesMapMakerTextInJSONOut.cpp so it knows where and how large to draw the
    //  circles
    std::string percentCirclesFilename = "foundPercentageCircles.txt";
    std::ofstream percentCirclesFile;
    percentCirclesFile.open(percentCirclesFilename);
    for (int percent = 1; percent <= 100; percent++) {
        const double desiredPopulation = (WORLD_POP_2015 / 100.0) * percent;
        std::cout << std::endl << "Now finding smallest possible circle with " << percent
                  << "\% of the world's population (" << ((long)desiredPopulation) << " people)"
                  << std::endl;
        SmallestCircleResult smallestCircle = 
            popData.smallestCircleWithGivenPopulation(desiredPopulation);
        std::cout << "Smallest possible circle with " << percent << "\% of the world's population ("
                  << ((long)desiredPopulation) << " people):" << std::endl;
        std::cout << "Population within " << smallestCircle.radius << " km of ("
                  << smallestCircle.lat << ", " << smallestCircle.lon << "): "
                  << ((long)smallestCircle.pop) << std::endl;
        // Used to be for entering into the python map making code, now just to match the format in
        //  the google doc
        std::cout << percent << ": (" << smallestCircle.radius << ", (" << std::setprecision(8)
                  << smallestCircle.lat << ", " << smallestCircle.lon << "))"
                  << std::setprecision(6) << std::endl;
        // TODO: Remove the magic number 6
        // TODO: Probably only wanna write if it's not already in there
        percentCirclesFile << ((int)percent) << " " << smallestCircle.radius << " "
                           << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION) << smallestCircle.lon
                           << " " << smallestCircle.lat << " " << smallestCircle.pop
                           << std::setprecision(6) << std::endl;
    }
    percentCirclesFile.close();
}

int main() {
    findPercentCircles();
}