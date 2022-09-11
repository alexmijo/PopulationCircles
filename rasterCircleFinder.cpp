// Uses GDAL library
// This is an early version of the program. It's still kinda spaghetti code and has code I copied
//  from random places on the internet unattributed (GeoTiff and most of distance()). If you're a
//  prospective employer, please look at BlackjackSim on my github instead for C++ stuff I've
//  written with better code style.

#include <algorithm>
#include <array>
#include <cmath>
#include <cpl_conv.h>
#include <fstream>
#include <gdal.h>
#include <gdal_priv.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

// Uses 2015 population data if false.
constexpr bool USE_2020_DATA = true;
constexpr double WORLD_POP = USE_2020_DATA ? 7757982599.3135586 : 7346242908.863955;
// See
//  https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout#comment99267684_554134
constexpr int DOUBLE_ROUND_TRIP_PRECISION =
    (std::numeric_limits<double>::digits10 == 15) ? 17 : std::numeric_limits<double>::digits10 + 3;
constexpr int EQUATOR_LEN = 40075;

struct CircleResult {
    // The lattitude and longitude of the center of the circle.
    double lat, lon;
    // The radius of the circle in kilometers.
    double radius;
    // The population within the circle.
    double pop;
};

// Defines a rectangular (on an equirectangular projection) region of the world, by the range of
//  latitudes and range of longitudes that the rectangle covers. Ranges are inclusive.
struct LatLonBoundaries {
    double leftLon, rightLon, upLat, downLat;
};

// TODO: Have this actually just be the data, and then a separate class has the actual circle
//  finding code. Could also then have another class with band finding code.
class RasterDataCircleFinder {
  private:
    constexpr static int KERNEL_WIDTH = 4; // Num cols in each kernel (4 corners of a box)
    int numRows, numCols;
    // First index is y (lat), second is x (lon). Matrix style indexing (top to bottom, left to
    //  right)
    // TODO: See if I can or should make this an array of arrays (library arrays). Can't be a 2D C
    //  array cause it's too large.
    // TODO: See if I should actually make this a 1D array (both for loading speed and using speed)
    double **sumTable;

    // Represents a previously found population circle, except instead of the radius it has whether
    //  not the result was a >= result cause the algorithm short circuited.
    struct CircleResultMaybeShortCircuit {
        // The lattitude and longitude of the center of the circle.
        // TODO: Now that the order is reversed, reverse it in the results file too.
        double lat, lon;
        // The population within the circle.
        double pop;
        // False iff this circle was found by the algorithm running to completion, true if the
        //  algorithm "short circuited" and was a >= result (i.e., there may be other circles with a
        //  bigger population but the same radius).
        bool shortCircuited;

        CircleResultMaybeShortCircuit(CircleResult result) {
            lat = result.lat;
            lon = result.lon;
            pop = result.pop;
            shortCircuited = false;
        }

        CircleResultMaybeShortCircuit() {}
    };

    // key is radius
    std::map<int, CircleResultMaybeShortCircuit> smallestCircleResults;
    std::string smallestCircleResultsFilename;

    // Defines a rectangular region of the raster data, by the rows and columns of the data that the
    //  rectangle covers. Ranges are inclusive.
    struct PixelBoundaries {
        int leftX, rightX, upY, downY;
    };

    // Used for boundingBoxEdge
    enum direction { north, south };

    PixelBoundaries pixelBoundariesFromLatLonBoundaries(const LatLonBoundaries &boundaries) {
        return PixelBoundaries{lonToX(boundaries.leftLon), lonToX(boundaries.rightLon),
                               latToY(boundaries.upLat), latToY(boundaries.downLat)};
    }

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

    // Makes the kernel for a specific lattitude.
    // TODO: Make the kernel an actual 2d array (maybe)
    // TODO: See if this double counts part of 1 column when wrapping around all
    // TODO: Get rid of the cenX parameter.
    std::vector<int> makeKernel(const int cenX, const int cenY, const double radius) {
        const double cenLon = lon(cenX);
        const double cenLat = lat(cenY);
        const int northEdge = boundingBoxEdge(cenX, cenY, radius, north);
        const int southEdge = boundingBoxEdge(cenX, cenY, radius, south);
        const int maxPossibleLength = southEdge - northEdge + 1;
        // A faster way of doing a 2D array of dimensions maxPossibleSize x KERNEL_WIDTH
        // Each row consists of: {westX, eastX, northY, southY} describing a summation table
        //  rectangle (so KERNEL_WIDTH must be 4) relative to cenX and cenY
        std::vector<int> kernel;

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
            } else {
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
                kernel.emplace_back(-horizontalOffset - 1);
                kernel.emplace_back(horizontalOffset);
                kernel.emplace_back(y - cenY - 1);
                if (y == southEdge) {
                    // If we just started a new rectangle at southEdge, it must end there too
                    kernel.emplace_back(y - cenY);
                    break;
                } else {
                    // Find where this new rectangle ends
                    while (true) {
                        y++;
                        if (y > southEdge) {
                            // Rectangle can't extend below south edge
                            kernel.emplace_back(y - cenY - 1);
                            break;
                        }

                        // Check if the circle has widened
                        if (horizontalOffset < numCols / 2) {
                            // Rectangles that wrap around the whole world can't widen any more
                            currLon = lon(cenX + horizontalOffset + 1);
                            currLat = lat(y);
                            if (distance(currLat, currLon, cenLat, cenLon) <= radius) {
                                // The circle has widened; the rectangle is done
                                kernel.emplace_back(y - cenY - 1);
                                break;
                            }
                        }

                        currLon = lon(cenX + horizontalOffset);
                        currLat = lat(y);
                        if (distance(currLat, currLon, cenLat, cenLon) > radius) {
                            // The y value can no longer be in the rectangle
                            kernel.emplace_back(y - cenY - 1);
                            break;
                        }
                    }
                }
            }
        }
        // This is ok because it will either be moved or copy-elisioned.
        return kernel;
    }

    // Returns the total population in the specified rectangle. West and north are 1 pixel outside
    //  of the rectangle, east and south are 1 pixel inside of the rectangle.
    // TODO: Make this not use conditional logic by just padding the north and west sides of
    //  popSumTable with 0s
    // TODO: If this makes the program slow, see if making it inline helps
    // TODO: See if it's actually even any slower to use PixelBoundaries here
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
    double popWithinKernel(const int cenX, const int cenY, const std::vector<int> &kernel) {
        double totalPop = 0;
        for (int kernelRow = 0; kernelRow < kernel.size() / KERNEL_WIDTH; kernelRow++) {
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
    double lon(int x) { return ((x + 0.5) / numCols) * 360.0 - 180.0; }

    // Get lattitude of the center of the <y>th (0 indexed) row
    double lat(int y) { return -(((y + 0.5) / numRows) * 180.0 - 90.0); }

    // The inverse of the lon function
    int lonToX(double lon) { return ((lon + 180.0) / 360.0) * numCols - 0.5; }

    // The inverse of the lat function
    int latToY(double lat) { return ((-lat + 90.0) / 180.0) * numRows - 0.5; }

    struct PixelCenterAndPop {
        // The center pixel of the circle
        int x, y;
        // The population of the circle
        double pop;

        PixelCenterAndPop(int x, int y, double pop) : x(x), y(y), pop(pop) {}

        // Makes it easy to sort PixelCenterAndPops by latitude.
        bool operator<(const PixelCenterAndPop &other) const { return y < other.y; }
    };

    // Passed in ranges are inclusive.
    // Adds all circles who's pop is at least cutoff * pop of the most populous circle to topCircles
    // Reassigns largestPop if necessary
    // Can short circuit and return early as soon as it finds a circle with at least desiredPop.
    void mostPopulousCirclesOfGivenRadiusPixelBoundaries(
        const double radius, const PixelBoundaries &boundaries, const int step,
        std::vector<PixelCenterAndPop> &topCircles, double &largestPop, const double cutoff,
        const int skipX, const int skipY, std::map<int, std::vector<int>> &kernels,
        const double desiredPop) {
        for (int cenY = boundaries.upY + step / 2; cenY <= boundaries.downY; cenY += step) {
            if (cenY < 0 || cenY >= numRows) {
                continue;
            }
            auto it = kernels.find(cenY);
            if (it == kernels.end()) {
                kernels[cenY] = makeKernel(1000, cenY, radius);
            }
            for (int cenX = boundaries.leftX + step / 2; cenX <= boundaries.rightX; cenX += step) {
                if (cenX < 0 || cenX >= numCols || (cenX == skipX && cenY == skipY)) {
                    continue;
                }
                double popWithinNKilometers = popWithinKernel(cenX, cenY, kernels[cenY]);
                if (popWithinNKilometers > largestPop * cutoff) {
                    topCircles.emplace_back(cenX, cenY, popWithinNKilometers);
                    if (popWithinNKilometers > largestPop) {
                        largestPop = popWithinNKilometers;
                        if (largestPop >= desiredPop) {
                            return;
                        }
                    }
                }
            }
        }
    }

    // Removes from topCircles any circles whose pop is less than largestPop * cutoff.
    // TODO: If this is slowing things down, use a linked list.
    void cutUnderperformingCircles(std::vector<PixelCenterAndPop> &topCircles,
                                   const double largestPop, const double cutoff) {
        std::vector<int> indicesToErase;
        for (int i = topCircles.size() - 1; i >= 0; i--) {
            if (topCircles[i].pop < largestPop * cutoff) {
                indicesToErase.push_back(i);
            }
        }
        for (int indexToErase : indicesToErase) {
            topCircles.erase(topCircles.begin() + indexToErase);
        }
    }

    // Finds the circle of the given radius containing maximal population, unless that would be a
    //  larger population that the passed in desiredPop, in which case this "short circuits" and the
    //  result will be the first circle (of the given radius) it found with a population greater
    //  than or equal to desiredPop.
    // Can optionally constrain the center of the circle to be within the given boundaries.
    // Radius given in kilometers.
    // The return value's shortCircuited field will be true iff this short circuited.
    // TODO: Make this less spagettiish
    // TODO: Doesnt check easternmost and southernmost column/row of the world
    CircleResultMaybeShortCircuit shortCircuitingMostPopulousCircleOfGivenRadius(
        const double radius, const double desiredPop,
        const LatLonBoundaries &boundaries = LatLonBoundaries{-180, 180, 90, -90}) {
        std::cout << "Radius: " << radius << std::endl;
        // TODO: Don't look at centers outside these ranges.
        // Turn lat and lon into indices in the pop data.
        PixelBoundaries pixelBoundaries = pixelBoundariesFromLatLonBoundaries(boundaries);

        // TODO: Make this nicer
        int initialStep;
        double cutoff256;
        double cutoff64;
        double cutoff16;
        double cutoff4;
        if (radius >= 14100) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.9) * 1;   // 0.913
            cutoff64 = 1 - (1 - 0.99) * 1;   // 0.9953
            cutoff16 = 1 - (1 - 0.9993) * 1; // 0.9999
            cutoff4 = 1 - (1 - 0.99993) * 1; // 0.999995
        } else if (radius >= 13910) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.8867) * 1;
            cutoff64 = 1 - (1 - 0.9633) * 1;
            cutoff16 = 1 - (1 - 0.9986) * 1;
            cutoff4 = 1 - (1 - 0.99985) * 1;
        } else if (radius >= 13750) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.88) * 1;
            cutoff64 = 1 - (1 - 0.95) * 1;
            cutoff16 = 1 - (1 - 0.9983) * 1;
            cutoff4 = 1 - (1 - 0.99984) * 1;
        } else if (radius >= 13625) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.8733) * 1;
            cutoff64 = 1 - (1 - 0.9433) * 1;
            cutoff16 = 1 - (1 - 0.9976) * 1;
            cutoff4 = 1 - (1 - 0.99983) * 1;
        } else if (radius >= 13500) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.87) * 1;
            cutoff64 = 1 - (1 - 0.94) * 1;
            cutoff16 = 1 - (1 - 0.9973) * 1;
            cutoff4 = 1 - (1 - 0.9998) * 1;
        } else if (radius >= 13250) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.8567) * 1;
            cutoff64 = 1 - (1 - 0.9333) * 1;
            cutoff16 = 1 - (1 - 0.9971) * 1;
            cutoff4 = 1 - (1 - 0.9996) * 1;
        } else if (radius >= 13000) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.85) * 1;
            cutoff64 = 1 - (1 - 0.93) * 1;
            cutoff16 = 1 - (1 - 0.997) * 1;
            cutoff4 = 1 - (1 - 0.9995) * 1;
        } else if (radius >= 12500) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.7833) * 1;
            cutoff64 = 1 - (1 - 0.9167) * 1;
            cutoff16 = 1 - (1 - 0.9967) * 1;
            cutoff4 = 1 - (1 - 0.9993) * 1;
        } else if (radius >= 12000) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.75) * 1;
            cutoff64 = 1 - (1 - 0.91) * 1;
            cutoff16 = 1 - (1 - 0.9965) * 1;
            cutoff4 = 1 - (1 - 0.9992) * 1;
        } else if (radius >= 11500) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.7167) * 1;
            cutoff64 = 1 - (1 - 0.9033) * 1;
            cutoff16 = 1 - (1 - 0.9962) * 1;
            cutoff4 = 1 - (1 - 0.9991) * 1;
        } else if (radius >= 11000) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.7) * 1;
            cutoff64 = 1 - (1 - 0.9) * 1;
            cutoff16 = 1 - (1 - 0.9961) * 1;
            cutoff4 = 1 - (1 - 0.9987) * 1;
        } else if (radius >= 9925) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.6667) * 1;
            cutoff64 = 1 - (1 - 0.8833) * 1;
            cutoff16 = 1 - (1 - 0.9914) * 1;
            cutoff4 = 1 - (1 - 0.9969) * 1;
        } else if (radius >= 8850) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.65) * 1;
            cutoff64 = 1 - (1 - 0.875) * 1;
            cutoff16 = 1 - (1 - 0.98905) * 1;
            cutoff4 = 1 - (1 - 0.9958) * 1;
        } else if (radius >= 7775) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.6167) * 1;
            cutoff64 = 1 - (1 - 0.865) * 1;
            cutoff16 = 1 - (1 - 0.9843) * 1;
            cutoff4 = 1 - (1 - 0.9937) * 1;
        } else if (radius >= 6700) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.6) * 1;
            cutoff64 = 1 - (1 - 0.86) * 1;
            cutoff16 = 1 - (1 - 0.982) * 1;
            cutoff4 = 1 - (1 - 0.9926) * 1;
        } else if (radius >= 5313) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.57) * 1;
            cutoff64 = 1 - (1 - 0.8511) * 1;
            cutoff16 = 1 - (1 - 0.9811) * 1;
            cutoff4 = 1 - (1 - 0.991) * 1;
        } else if (radius >= 3925) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.555) * 1;
            cutoff64 = 1 - (1 - 0.8467) * 1;
            cutoff16 = 1 - (1 - 0.9807) * 1;
            cutoff4 = 1 - (1 - 0.9902) * 1;
        } else if (radius >= 2538) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.525) * 1;
            cutoff64 = 1 - (1 - 0.8422) * 1;
            cutoff16 = 1 - (1 - 0.9802) * 1;
            cutoff4 = 1 - (1 - 0.9894) * 1;
        } else if (radius >= 1150) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.51) * 1;
            cutoff64 = 1 - (1 - 0.84) * 1;
            cutoff16 = 1 - (1 - 0.98) * 1;
            cutoff4 = 1 - (1 - 0.989) * 1;
        } else if (radius >= 725) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.17) * 1;
            cutoff64 = 1 - (1 - 0.8333) * 1;
            cutoff16 = 1 - (1 - 0.97) * 1;
            cutoff4 = 1 - (1 - 0.98) * 1;
        } else if (radius >= 300) {
            initialStep = 64;
            cutoff64 = 1 - (1 - 0.83) * 1;
            cutoff16 = 1 - (1 - 0.965) * 1;
            cutoff4 = 1 - (1 - 0.973) * 1;
        } else if (radius >= 100) {
            initialStep = 16;
            cutoff16 = 1 - (1 - 0.8) * 1;
            cutoff4 = 1 - (1 - 0.965) * 1;
        } else if (radius >= 20) {
            initialStep = 4;
            cutoff4 = 1 - (1 - 0.7) * 1;
        } else {
            initialStep = 1;
        }

        int step = initialStep; // Must be a power of 4, I think
        double cutoff = cutoff256;
        if (step == 1) {
            // Only want to keep most populous circle when step is 1
            cutoff = 1;
        } else if (step <= 4) {
            cutoff = cutoff4;
        } else if (step <= 16) {
            cutoff = cutoff16;
        } else if (step <= 64) {
            cutoff = cutoff64;
        }
        std::vector<PixelCenterAndPop> topCircles;
        double largestPop = 0;
        std::map<int, std::vector<int>> kernels;

        std::cout << "step: " << step << std::endl;
        // -100 for skipY and skipX since we don't want to skip any pixels
        mostPopulousCirclesOfGivenRadiusPixelBoundaries(radius, pixelBoundaries, step, topCircles,
                                                        largestPop, cutoff, -100, -100, kernels,
                                                        desiredPop);
        cutUnderperformingCircles(topCircles, largestPop, cutoff);
        CircleResultMaybeShortCircuit result;
        if (largestPop >= desiredPop) {
            for (const auto &circle : topCircles) {
                if (circle.pop == largestPop) {
                    result.lat = lat(circle.y);
                    result.lon = lon(circle.x);
                    result.pop = largestPop;
                    result.shortCircuited = true;
                    std::cout << "(SS)Population within " << radius << " kilometers of ("
                              << result.lat << ", " << result.lon << "): " << ((long)result.pop)
                              << std::endl;
                    std::cout << ">= " << ((long)desiredPop) << ". Short circuiting." << std::endl;
                    return result;
                }
            }
        }
        while (step > 1) {
            step /= 4;
            std::cout << "topCircles.size(): " << topCircles.size() << std::endl;
            std::cout << "largestPop: " << ((long)largestPop) << std::endl;
            std::cout << "step: " << step << std::endl;
            if (step == 1) {
                // Only want to keep most populous circle when step is 1
                cutoff = 1;
            } else if (step <= 4) {
                cutoff = cutoff4;
            } else if (step <= 16) {
                cutoff = cutoff16;
            } else if (step <= 64) {
                cutoff = cutoff64;
            }
            // Sorting by lattitude allows kernels to be reused between calls to
            //  mostPopulousCirclesOfGivenRadiusPixelBoundaries().
            std::sort(topCircles.begin(), topCircles.end());
            const int numCirclesAtStartOfLoop = topCircles.size();
            int currY = -1;
            for (int i = 0; i < numCirclesAtStartOfLoop; i++) {
                if (topCircles[i].y != currY) {
                    kernels.clear();
                    currY = topCircles[i].y;
                }
                PixelBoundaries sectionBoundaries{
                    topCircles[i].x - step * 2, topCircles[i].x + step * 2 - 1,
                    topCircles[i].y - step * 2, topCircles[i].y + step * 2 - 1};
                mostPopulousCirclesOfGivenRadiusPixelBoundaries(
                    radius, sectionBoundaries, step, topCircles, largestPop, cutoff,
                    topCircles[i].x, topCircles[i].y, kernels, desiredPop);
                if (largestPop >= desiredPop) {
                    for (const auto &circle : topCircles) {
                        if (circle.pop == largestPop) {
                            result.lat = lat(circle.y);
                            result.lon = lon(circle.x);
                            result.pop = largestPop;
                            result.shortCircuited = true;
                            std::cout << "(SS)Population within " << radius << " kilometers of ("
                                      << result.lat << ", " << result.lon
                                      << "): " << ((long)result.pop) << std::endl;
                            std::cout << ">= " << ((long)desiredPop) << ". Short circuiting."
                                      << std::endl;
                            return result;
                        }
                    }
                }
            }
            cutUnderperformingCircles(topCircles, largestPop, cutoff);
        }
        // Only the most populous circle remains in the list, since the cutoff for step size 1 is 1.
        result.lat = lat(topCircles[0].y);
        result.lon = lon(topCircles[0].x);
        if (topCircles[0].pop != largestPop) {
            std::cout << "topCircles[0].pop: " << topCircles[0].pop << std::endl;
            std::cout << "largestPop: " << largestPop << std::endl;
            throw std::logic_error(
                "Should never get here. largestPop should be first element of topCircles.");
        }
        result.pop = largestPop;
        result.shortCircuited = false;
        std::cout << "Population within " << radius << " kilometers of (" << result.lat << ", "
                  << result.lon << "): " << ((long)result.pop) << std::endl;
        return result;
    }

  public:
    // File specified by <sumTableFilename> must consist of an int for numRows, and int for
    //  numCols, and then numRows * numCols doubles representing all of the data in the summation
    //  table.
    // Projection must be equirectangular.
    // TODO: Complete spec
    RasterDataCircleFinder(const std::string &sumTableFilename,
                           const std::string &smallestCircleResultsFilename)
        : smallestCircleResultsFilename(smallestCircleResultsFilename) {
        std::fstream sumTableFile;
        sumTableFile.open(sumTableFilename, std::ios::in | std::ios::binary);
        sumTableFile.read(reinterpret_cast<char *>(&numRows), sizeof(int));
        sumTableFile.read(reinterpret_cast<char *>(&numCols), sizeof(int));
        // TODO: Either use a unique_pointer to array of unique_pointers to arrays of ints, a vector
        //  of vectors, or the 1D versions of either of those (with an indexing function) for speed.
        //  No need for new and delete.
        sumTable = new double *[numRows];
        for (int r = 0; r < numRows; r++) {
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
            std::cout << smallestCircleResultsFilename
                      << " is empty or doesn't exist. Making it non-empty and existing."
                      << std::endl;
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
            CircleResultMaybeShortCircuit result;
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
                    throw std::runtime_error(
                        "Incorrectly formatted smallest circle results file: " +
                        smallestCircleResultsFilename);
                }
                result.shortCircuited = true;
            }
            smallestCircleResults[radius] = result;
        }
        smallestCircleResultsFile.close();
    }

    ~RasterDataCircleFinder() {
        for (int r = 0; r < numRows; r++) {
            delete[] sumTable[r];
        }
        delete[] sumTable;
    }

    // Returns distance in kilometers between two points on Earth's surface
    // TODO: Find source and give credit
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
        double delta_sigma = B * sin_sigma *
                             (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq) -
                                                          (B / 6.0) * cos_2sigma_m *
                                                              (-3 + 4 * sin_sigma * sin_sigma) *
                                                              (-3 + 4 * cos_2sigma_m_sq)));

        // Distance
        double s = b * A * (sigma - delta_sigma);

        // Convert from meters to kilometers before returning
        return s / 1000.;
    }

    int getNumRows() { return numRows; }

    int getNumCols() { return numCols; }

    // Finds the circle of the given radius containing maximal population.
    // Can optionally constrain the center of the circle to be within the given boundaries.
    // Radius given in kilometers.
    CircleResult mostPopulousCircleOfGivenRadius(
        const double radius,
        const LatLonBoundaries &boundaries = LatLonBoundaries{-180, 180, 90, -90}) {
        // Maximally large desiredPop parameter is passed so that it won't short circuit
        CircleResultMaybeShortCircuit shouldntHaveShortCircuited =
            shortCircuitingMostPopulousCircleOfGivenRadius(
                radius, std::numeric_limits<double>::max(), boundaries);
        return CircleResult{shouldntHaveShortCircuited.lat, shouldntHaveShortCircuited.lon, radius,
                            shouldntHaveShortCircuited.pop};
    }

    // TODO: Make a spec comment for this
    // Kind of dangerous to pass in the last 4 args to this function, since the
    //  smallestCircleResultsFile might get spurious data added to it then
    // TODO: Figure out why largestPop seems to stay the same after searching at a step that's 1/4
    //  the size something like 1/4 of the time rather than 1/16 of the time like expected.
    CircleResult smallestCircleWithGivenPopulation(
        const double pop,
        const LatLonBoundaries boundaries = LatLonBoundaries{-180, 180, 90, -90}) {
        // Radii will be ints cause only interested in getting to the nearest kilometer
        int upperBound = EQUATOR_LEN / 2;
        int lowerBound = 0;
        CircleResult result;
        bool foundSuitableCircle = false; // Haven't yet found a circle with sufficiently large pop
        // A "soft" upper bound is one which comes from a short circuited previous result
        bool softUpperBound = false;
        // Tighten bounds as much as possible using previous results
        for (auto it = smallestCircleResults.begin(); it != smallestCircleResults.end(); it++) {
            if (it->second.pop < pop && it->first >= lowerBound && !it->second.shortCircuited) {
                lowerBound = it->first;
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
        int radius = (upperBound + lowerBound) / 2;
        if (smallestCircleResults[upperBound].pop / 3.0 +
                    smallestCircleResults[lowerBound].pop * 2.0 / 3.0 >
                pop &&
            upperBound - lowerBound >= 3) {
            // + 0.00001 just in case floating point error
            radius = upperBound / 3.0 + lowerBound * 2.0 / 3.0 + 0.00001;
        } else if (smallestCircleResults[upperBound].pop * 2.0 / 3.0 +
                           smallestCircleResults[lowerBound].pop / 3.0 <
                       pop &&
                   upperBound - lowerBound >= 3) {
            // + 0.00001 just in case floating point error
            radius = upperBound * 2.0 / 3.0 + lowerBound / 3.0 + 0.00001;
        }
        if (smallestCircleResults[upperBound].pop / 4.0 +
                    smallestCircleResults[lowerBound].pop * 3.0 / 4.0 >
                pop &&
            upperBound - lowerBound >= 4) {
            // + 0.00001 just in case floating point error
            radius = upperBound / 4.0 + lowerBound * 3.0 / 4.0 + 0.00001;
        } else if (smallestCircleResults[upperBound].pop * 3.0 / 4.0 +
                           smallestCircleResults[lowerBound].pop / 4.0 <
                       pop &&
                   upperBound - lowerBound >= 4) {
            // + 0.00001 just in case floating point error
            radius = upperBound * 3.0 / 4.0 + lowerBound / 4.0 + 0.00001;
        }
        while (upperBound - lowerBound > 1 || softUpperBound) {
            CircleResultMaybeShortCircuit largestSumCircle;
            if (upperBound - lowerBound <= 3) {
                if (upperBound - lowerBound == 1) {
                    // Need to make upperBound not soft
                    radius = upperBound;
                }
                // Huge desired pop since we don't want it to short circuit
                largestSumCircle = mostPopulousCircleOfGivenRadius(radius, boundaries);
            } else {
                largestSumCircle =
                    shortCircuitingMostPopulousCircleOfGivenRadius(radius, pop, boundaries);
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
            if (!previousResultsConsistent()) {
                throw std::logic_error("Previous results inconsistent.");
            }
            radius = (upperBound + lowerBound) / 2;
            if (smallestCircleResults[upperBound].pop / 3.0 +
                        smallestCircleResults[lowerBound].pop * 2.0 / 3.0 >
                    pop &&
                upperBound - lowerBound >= 3) {
                // + 0.00001 just in case floating point error
                radius = upperBound / 3.0 + lowerBound * 2.0 / 3.0 + 0.00001;
            } else if (smallestCircleResults[upperBound].pop * 2.0 / 3.0 +
                               smallestCircleResults[lowerBound].pop / 3.0 <
                           pop &&
                       upperBound - lowerBound >= 3) {
                // + 0.00001 just in case floating point error
                radius = upperBound * 2.0 / 3.0 + lowerBound / 3.0 + 0.00001;
            }
            if (smallestCircleResults[upperBound].pop / 4.0 +
                        smallestCircleResults[lowerBound].pop * 3.0 / 4.0 >
                    pop &&
                upperBound - lowerBound >= 4) {
                // + 0.00001 just in case floating point error
                radius = upperBound / 4.0 + lowerBound * 3.0 / 4.0 + 0.00001;
            } else if (smallestCircleResults[upperBound].pop * 3.0 / 4.0 +
                               smallestCircleResults[lowerBound].pop / 4.0 <
                           pop &&
                       upperBound - lowerBound >= 4) {
                // + 0.00001 just in case floating point error
                radius = upperBound * 3.0 / 4.0 + lowerBound / 4.0 + 0.00001;
            }
        }
        if (!foundSuitableCircle) {
            // TODO: Fail fast
            if (pop > WORLD_POP) {
                throw std::invalid_argument("Desired pop is larger than the world population, and a"
                                            " suitable circle wasn't found");
            }
            throw std::logic_error("Should never get here");
        }
        return result;
    }

    bool previousResultsConsistent() {
        double maxPop = -1;
        int maxPopRadius = -1;
        for (const auto &[key, value] : smallestCircleResults) {
            if (value.pop < maxPop && !value.shortCircuited) {
                std::cout << "Previous results inconsistent:" << std::endl;
                std::cout << maxPopRadius << ": " << ((long)maxPop) << std::endl;
                std::cout << key << ": " << ((long)value.pop) << std::endl;
                return false;
            }
            if (value.pop > maxPop) {
                maxPop = value.pop;
                maxPopRadius = key;
            }
        }
        return true;
    }
};

// See what I can get away with
void testCircleSkipping() {
    double radius = 11029;

    std::cout << "Loading population summation table." << std::endl;
    std::string sumTableFilename = "popSumTable.bin";
    std::string smallestCircleResultsFilename = "popSmallestCircleResults.txt";
    RasterDataCircleFinder popData(sumTableFilename, smallestCircleResultsFilename);
    std::cout << "Loaded population summation table." << std::endl;

    CircleResult smallestCircle = popData.mostPopulousCircleOfGivenRadius(radius);
    std::cout << "Population within " << radius << " km of (" << smallestCircle.lat << ", "
              << smallestCircle.lon << "): " << (long long)smallestCircle.pop << std::endl;
}

// Used for finding the smallest circles containing 1-100% of the world's population.
// TODO: Allow for keyboard interrupt for manual entering of radius.
void findPercentCircles() {
    std::cout << "Loading population summation table." << std::endl;
    if (USE_2020_DATA) {
        std::cout << "Using 2020 data." << std::endl;
    } else {
        std::cout << "Using 2015 data." << std::endl;
    }
    std::string sumTableFilename = "popSumTable.bin";
    if (USE_2020_DATA) {
        sumTableFilename = "popSumTable2020.bin";
    }
    std::string smallestCircleResultsFilename = "popSmallestCircleResults.txt";
    if (USE_2020_DATA) {
        smallestCircleResultsFilename = "popSmallestCircleResults2020.txt";
    }
    RasterDataCircleFinder popData(sumTableFilename, smallestCircleResultsFilename);
    std::cout << "Loaded population summation table." << std::endl;
    if (!popData.previousResultsConsistent()) {
        return;
    }
    // Used by imageManipulation.java so that it knows what text to add, as well as by
    //  percentageCirclesMapMakerTextInJSONOut.cpp so it knows where and how large to draw the
    //  circles
    std::string percentCirclesFilename = "foundPercentageCircles.txt";
    if (USE_2020_DATA) {
        percentCirclesFilename = "foundPercentageCircles2020.txt";
    }
    std::ofstream percentCirclesFile;
    percentCirclesFile.open(percentCirclesFilename);
    for (double percent = 100; percent > 0.09; percent -= 0.1) {
        double desiredPopulation = (WORLD_POP / 100.0) * percent;
        if (percent > 99.99) {
            desiredPopulation = (long long)WORLD_POP;
        }
        std::cout << std::endl
                  << "Now finding smallest possible circle with " << percent
                  << "\% of the world's population (" << ((long)desiredPopulation) << " people)"
                  << std::endl;
        CircleResult smallestCircle;
        smallestCircle = popData.smallestCircleWithGivenPopulation(desiredPopulation);
        std::cout << "Smallest possible circle with " << percent << "\% of the world's population ("
                  << ((long)desiredPopulation) << " people):" << std::endl;
        std::cout << "Population within " << smallestCircle.radius << " km of ("
                  << smallestCircle.lat << ", " << smallestCircle.lon
                  << "): " << ((long)smallestCircle.pop) << std::endl;
        // Used to be for entering into the python map making code, now just to match the format in
        //  the google doc
        std::cout << percent << ": (" << smallestCircle.radius << ", (" << std::setprecision(8)
                  << smallestCircle.lat << ", " << smallestCircle.lon << "))"
                  << std::setprecision(6) << std::endl;
        // TODO: Remove the magic number 6
        // TODO: Probably only wanna write if it's not already in there
        percentCirclesFile << percent << " " << smallestCircle.radius << " "
                           << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION) << smallestCircle.lon
                           << " " << smallestCircle.lat << " " << smallestCircle.pop
                           << std::setprecision(6) << std::endl;
    }
    percentCirclesFile.close();
}

int main() { findPercentCircles(); }