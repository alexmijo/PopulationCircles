// Uses GDAL library
// This is an early version of the program. It's still kinda spaghetti code. If you're a
//  prospective employer, please look at BlackjackSim on my github instead for C++ stuff I've
//  written with better code style.

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <unordered_set>
#include <vector>
#include <chrono>

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

struct Pixel {
    int x, y;
};

// Defines a rectangular region of the raster data, by the rows and columns of the data that the
//  rectangle covers. Ranges are inclusive.
class PixelBoundaries {
  public:
    int leftX, rightX, upY, downY;

    bool contains(Pixel p) const {
        return leftX <= p.x && p.x <= rightX && upY <= p.y && p.y <= downY;
    }
};

// Defines a rectangular (on an equirectangular projection) region of the world, by the range of
//  latitudes and range of longitudes that the rectangle covers. Ranges are inclusive.
struct LatLonBoundaries {
    double leftLon, rightLon, upLat, downLat;
};

class RasterDataCircleFinder;

class CirclePixelSet {
  public:
    CirclePixelSet(const RasterDataCircleFinder &, CircleResult);
    CirclePixelSet(const RasterDataCircleFinder &, Pixel, double);
    bool containsY(const int) const;
    bool contains(Pixel) const;
    PixelBoundaries boundingBox;

  private:
    Pixel center;
    // Inclusive on both sides
    std::multimap<int, std::pair<int, int>> xRanges;
};

// TODO: Have this actually just be the data, and then a separate class has the actual circle
//  finding code. Could also then have another class with band finding code.
class RasterDataCircleFinder {
    friend class CirclePixelSet;

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

    // Used for boundingBoxEdge
    enum direction { north, south };

    PixelBoundaries pixelBoundariesFromLatLonBoundaries(const LatLonBoundaries &boundaries) {
        return PixelBoundaries{lonToX(boundaries.leftLon), lonToX(boundaries.rightLon),
                               latToY(boundaries.upLat), latToY(boundaries.downLat)};
    }

    // Returns the northernmost or southernmost row containing any pixel within <radius> km of the
    //  given pixel (specified by <x> and <y>).
    int boundingBoxEdge(const int x, const int y, const double radius,
                        direction northOrSouth) const {
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
    inline int kernelIndex(const int row, const int col) const { return KERNEL_WIDTH * row + col; }

    // Makes the kernel for a specific lattitude.
    // TODO: Make the kernel an actual 2d array (maybe)
    // TODO: See if this double counts part of 1 column when wrapping around all
    // TODO: Get rid of the cenX parameter.
    std::vector<int> makeKernel(const int cenX, const int cenY, const double radius) const {
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

    std::vector<int> makeKernel(const int cenY, const double radius) const {
        // TODO: Magic number.
        return makeKernel(1000, cenY, radius);
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
    double lon(int x) const { return ((x + 0.5) / numCols) * 360.0 - 180.0; }

    // Get lattitude of the center of the <y>th (0 indexed) row
    double lat(int y) const { return -(((y + 0.5) / numRows) * 180.0 - 90.0); }

    // The inverse of the lon function
    int lonToX(double lon) const { return ((lon + 180.0) / 360.0) * numCols; }

    // The inverse of the lat function
    int latToY(double lat) const { return ((-lat + 90.0) / 180.0) * numRows; }

    struct PixelCenterAndPop {
        // The center pixel of the circle
        int x, y;
        // The population of the circle
        double pop;

        PixelCenterAndPop(int x, int y, double pop) : x(x), y(y), pop(pop) {}

        // Makes it easy to sort PixelCenterAndPops by latitude.
        bool operator<(const PixelCenterAndPop &other) const {
            return y < other.y || (y == other.y && x < other.x);
        }
    };

    struct intPairHash {
        size_t operator()(const std::pair<int, int> &intPair) const {
            std::hash<int> intHash;
            return intHash(intPair.first) * 31 + intHash(intPair.second);
        }
    };

    // Passed in ranges are inclusive.
    // Adds all circles who's pop is at least cutoff * pop of the most populous circle to topCircles
    // Reassigns largestPop if necessary
    // Can short circuit and return early as soon as it finds a circle with at least desiredPop.
    void mostPopulousCirclesOfGivenRadiusPixelBoundaries(
        const double radius, const PixelBoundaries &boundaries, const int step,
        std::set<PixelCenterAndPop> &topCircles, double &largestPop, const double cutoff,
        std::map<int, std::vector<int>> &kernels, const double desiredPop,
        std::unordered_set<std::pair<int, int>, intPairHash> &alreadyChecked,
        bool useAlreadyChecked = true) {
        for (int cenY = boundaries.upY + step / 2; cenY <= boundaries.downY; cenY += step) {
            if (cenY < 0 || cenY >= numRows) {
                continue;
            }
            bool checkedForKernel = false;
            std::vector<int> kernel;
            for (int cenX = boundaries.leftX + step / 2; cenX <= boundaries.rightX; cenX += step) {
                if (cenX < 0 || cenX >= numCols ||
                    (useAlreadyChecked && alreadyChecked.find(std::pair<int, int>(cenX, cenY)) !=
                                              alreadyChecked.end())) {
                    continue;
                }
                if (!checkedForKernel) {
                    auto it = kernels.find(cenY);
                    if (it == kernels.end()) {
                        if (step <= 4) {
                            kernel = makeKernel(1000, cenY, radius);
                        } else {
                            kernels[cenY] = makeKernel(1000, cenY, radius);
                        }
                    }
                    checkedForKernel = true;
                }
                double popWithinNKilometers = (step <= 4)
                                                  ? popWithinKernel(cenX, cenY, kernel)
                                                  : popWithinKernel(cenX, cenY, kernels[cenY]);
                if (useAlreadyChecked) {
                    alreadyChecked.emplace(cenX, cenY);
                }
                if (popWithinNKilometers > largestPop * cutoff) {
                    topCircles.emplace(cenX, cenY, popWithinNKilometers);
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

    // TODO: Maybe there are rather easy ways to reduce work here for speed. Definitely should
    // clean.
    double getRectanglesOverlappingPop(const int west, const int east, const int north,
                                       const int south,
                                       const std::vector<CirclePixelSet> &largestCitiesPixels,
                                       const std::vector<int> &couldIntersectWithIndices) {
        double overlappingPop = 0;
        for (int y = north + 1; y <= south; y++) {
            bool skipThisY = true;
            for (int i : couldIntersectWithIndices) {
                const CirclePixelSet &cityPixels = largestCitiesPixels[i];
                if (largestCitiesPixels[i].containsY(y)) {
                    skipThisY = false;
                    break;
                }
            }
            if (skipThisY) {
                continue;
            }
            for (int x = west + 1; x <= east; x++) {
                // Could add larger rectangles at a time in order to try to speed things up.
                for (int i : couldIntersectWithIndices) {
                    if (largestCitiesPixels[i].contains({x, y})) {
                        overlappingPop += popWithinRectangle(x - 1, x, y - 1, y);
                        break;
                    }
                }
            }
        }
        return overlappingPop;
    }

    // TODO: Clean/spec
    double getOverlappingPop(int cenX, int cenY, const std::vector<int> &kernel,
                             PixelBoundaries boundingBox,
                             const std::vector<CirclePixelSet> &largestCitiesPixels) {
        std::vector<int> couldIntersectWithIndices;
        for (int i = 0; i < largestCitiesPixels.size(); i++) {
            if ((boundingBox.rightX >= largestCitiesPixels[i].boundingBox.leftX &&
                 boundingBox.leftX <= largestCitiesPixels[i].boundingBox.rightX) &&
                (boundingBox.downY >= largestCitiesPixels[i].boundingBox.upY &&
                 boundingBox.upY <= largestCitiesPixels[i].boundingBox.downY)) {
                couldIntersectWithIndices.emplace_back(i);
            }
        }
        if (couldIntersectWithIndices.empty()) {
            return 0;
        }
        double overlappingPop = 0;
        for (int kernelRow = 0; kernelRow < kernel.size() / KERNEL_WIDTH; kernelRow++) {
            // Sides of the rectangle
            const int west = kernel[kernelIndex(kernelRow, 0)] + cenX;
            const int east = kernel[kernelIndex(kernelRow, 1)] + cenX;
            const int north = kernel[kernelIndex(kernelRow, 2)] + cenY;
            const int south = kernel[kernelIndex(kernelRow, 3)] + cenY;
            if (kernel[kernelIndex(kernelRow, 1)] == numCols / 2) {
                // This rectangle encircles the entire latitude
                overlappingPop += getRectanglesOverlappingPop(
                    -1, numCols - 1, north, south, largestCitiesPixels, couldIntersectWithIndices);
            } else if (west < -1) {
                // Need to wrap around the antimeridian, creating two rectangles
                overlappingPop += getRectanglesOverlappingPop(
                    -1, east, north, south, largestCitiesPixels, couldIntersectWithIndices);
                overlappingPop +=
                    getRectanglesOverlappingPop(numCols + west, numCols - 1, north, south,
                                                largestCitiesPixels, couldIntersectWithIndices);
            } else if (east >= numCols) {
                // Need to wrap around the antimeridian, creating two rectangles
                overlappingPop +=
                    getRectanglesOverlappingPop(west, numCols - 1, north, south,
                                                largestCitiesPixels, couldIntersectWithIndices);
                overlappingPop +=
                    getRectanglesOverlappingPop(-1, east - numCols, north, south,
                                                largestCitiesPixels, couldIntersectWithIndices);
            } else {
                overlappingPop += getRectanglesOverlappingPop(
                    west, east, north, south, largestCitiesPixels, couldIntersectWithIndices);
            }
        }
        return overlappingPop;
    }

    // TODO: Spec
    int maxEastWestExtent(const std::vector<int> &kernel) {
        int maxExtent = 0;
        for (int kernelRow = 0; kernelRow < kernel.size() / KERNEL_WIDTH; kernelRow++) {
            maxExtent = std::max(kernel[kernelIndex(kernelRow, 1)], maxExtent); // East
        }
        return maxExtent;
    }

    void mostPopulousCirclesOfGivenRadiusPixelBoundariesForFindingLargestCities(
        const double radius, const PixelBoundaries &boundaries, const int step,
        std::set<PixelCenterAndPop> &topCircles, double &largestPop, const double cutoff,
        std::map<int, std::vector<int>> &kernels,
        const std::vector<CirclePixelSet> &largestCitiesPixels,
        std::unordered_set<std::pair<int, int>, intPairHash> &alreadyChecked,
        bool useAlreadyChecked = true) {
        for (int cenY = boundaries.upY + step / 2; cenY <= boundaries.downY; cenY += step) {
            if (cenY < 0 || cenY >= numRows) {
                continue;
            }
            bool checkedForKernel = false;
            std::vector<int> kernel;
            PixelBoundaries boundingBox;
            int maxEWExtent;
            for (int cenX = boundaries.leftX + step / 2; cenX <= boundaries.rightX; cenX += step) {
                if (cenX < 0 || cenX >= numCols ||
                    (useAlreadyChecked && alreadyChecked.find(std::pair<int, int>(cenX, cenY)) !=
                                              alreadyChecked.end())) {
                    continue;
                }
                if (!checkedForKernel) {
                    auto it = kernels.find(cenY);
                    if (it == kernels.end()) {
                        if (step <= 4) {
                            kernel = makeKernel(1000, cenY, radius);
                            maxEWExtent = maxEastWestExtent(kernel);
                            if (cenX - maxEWExtent >= 0 && cenX + maxEWExtent < numCols) {
                                // Doesn't cross the antimeridian
                                boundingBox = PixelBoundaries{
                                    cenX - maxEWExtent, cenX + maxEWExtent,
                                    cenY + kernel[kernelIndex(0, 2)] + 1,
                                    cenY + kernel[kernelIndex(kernel.size() / KERNEL_WIDTH - 1, 3)]};
                            } else {
                                boundingBox = PixelBoundaries{
                                    cenX - maxEWExtent, cenX + maxEWExtent,
                                    cenY + kernel[kernelIndex(0, 2)] + 1,
                                    cenY + kernel[kernelIndex(kernel.size() / KERNEL_WIDTH - 1, 3)]};
                            }
                        } else {
                            kernels[cenY] = makeKernel(1000, cenY, radius);
                            maxEWExtent = maxEastWestExtent(kernels[cenY]);
                            boundingBox = PixelBoundaries{
                                cenX - maxEWExtent, cenX + maxEWExtent,
                                cenY + kernels[cenY][kernelIndex(0, 2)] + 1,
                                cenY + kernels[cenY][kernelIndex(
                                           kernels[cenY].size() / KERNEL_WIDTH - 1, 3)]};
                        }
                    }
                    checkedForKernel = true;
                }
                boundingBox.leftX = cenX - maxEWExtent;
                boundingBox.rightX = cenX + maxEWExtent;
                double popWithinNKilometers = (step <= 4)
                                                  ? popWithinKernel(cenX, cenY, kernel)
                                                  : popWithinKernel(cenX, cenY, kernels[cenY]);
                popWithinNKilometers -=
                    (step <= 4)
                        ? getOverlappingPop(cenX, cenY, kernel, boundingBox, largestCitiesPixels)
                        : getOverlappingPop(cenX, cenY, kernels[cenY], boundingBox,
                                            largestCitiesPixels);
                if (useAlreadyChecked) {
                    alreadyChecked.emplace(cenX, cenY);
                }
                if (popWithinNKilometers > largestPop * cutoff) {
                    topCircles.emplace(cenX, cenY, popWithinNKilometers);
                    if (popWithinNKilometers > largestPop) {
                        largestPop = popWithinNKilometers;
                    }
                }
            }
        }
    }

    // Removes from topCircles any circles whose pop is less than largestPop * cutoff.
    // TODO: If this is slowing things down, use a linked list.
    void cutUnderperformingCircles(std::set<PixelCenterAndPop> &topCircles, const double largestPop,
                                   const double cutoff) {
        for (auto it = topCircles.begin(); it != topCircles.end();) {
            if (it->pop < largestPop * cutoff) {
                it = topCircles.erase(it);
            } else {
                ++it;
            }
        }
    }

    // TODO: Return these with a tuple and structured binding.
    void setCutoffsAndInitialStep(int radius, double &cutoff256, double &cutoff64, double &cutoff16,
                                  double &cutoff4, int &initialStep) {
        constexpr double diff256 = 0.207;
        constexpr double diff64 = 0.142;
        constexpr double diff16 = 0.142;
        constexpr double diff4 = 0.148;
        if (radius >= 14100) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.9) * 0.9;
            cutoff64 = 1 - (1 - 0.99) * 0.21;
            cutoff16 = 1 - (1 - 0.9993) * 0.19;
            cutoff4 = 1 - (1 - 0.99993) * 0.15;
        } else if (radius >= 13910) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.8867) * 0.9;
            cutoff64 = 1 - (1 - 0.9633) * 0.22;
            cutoff16 = 1 - (1 - 0.9986) * 0.2;
            cutoff4 = 1 - (1 - 0.99985) * 0.15;
        } else if (radius >= 13750) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.88) * 0.9;
            cutoff64 = 1 - (1 - 0.95) * 0.195;
            cutoff16 = 1 - (1 - 0.9983) * 0.18;
            cutoff4 = 1 - (1 - 0.99984) * 0.14;
        } else if (radius >= 13625) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.8733) * 0.9;
            cutoff64 = 1 - (1 - 0.9433) * 0.2;
            cutoff16 = 1 - (1 - 0.9976) * 0.185;
            cutoff4 = 1 - (1 - 0.99983) * 0.145;
        } else if (radius >= 13500) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.87) * 0.63;
            cutoff64 = 1 - (1 - 0.94) * 0.18;
            cutoff16 = 1 - (1 - 0.9973) * 0.19;
            cutoff4 = 1 - (1 - 0.9998) * 0.145;
        } else if (radius >= 13250) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.8567) * 0.64;
            cutoff64 = 1 - (1 - 0.9333) * 0.1825;
            cutoff16 = 1 - (1 - 0.9971) * 0.19;
            cutoff4 = 1 - (1 - 0.9996) * 0.145;
        } else if (radius >= 13000) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.85) * 0.65;
            cutoff64 = 1 - (1 - 0.93) * 0.185;
            cutoff16 = 1 - (1 - 0.997) * 0.19;
            cutoff4 = 1 - (1 - 0.9995) * 0.145;
        } else if (radius >= 12500) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.7833) * 0.66;
            cutoff64 = 1 - (1 - 0.9167) * 0.19;
            cutoff16 = 1 - (1 - 0.9967) * 0.19;
            cutoff4 = 1 - (1 - 0.9993) * 0.145;
        } else if (radius >= 12000) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.75) * 0.67;
            cutoff64 = 1 - (1 - 0.91) * 0.195;
            cutoff16 = 1 - (1 - 0.9965) * 0.19;
            cutoff4 = 1 - (1 - 0.9992) * 0.145;
        } else if (radius >= 11500) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.7167) * 0.677;
            cutoff64 = 1 - (1 - 0.9033) * 0.202;
            cutoff16 = 1 - (1 - 0.9962) * 0.187;
            cutoff4 = 1 - (1 - 0.9991) * 0.142;
        } else if (radius >= 11000) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.7) * 0.676;
            cutoff64 = 1 - (1 - 0.9) * 0.211;
            cutoff16 = 1 - (1 - 0.9961) * 0.186;
            cutoff4 = 1 - (1 - 0.9987) * 0.141;
        } else if (radius >= 9925) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.6667) * (0.676 - diff256);
            cutoff64 = 1 - (1 - 0.8833) * (0.231 - diff64);
            cutoff16 = 1 - (1 - 0.9914) * (0.2 - diff16);
            cutoff4 = 1 - (1 - 0.9969) * (0.169 - diff4);
        } else if (radius >= 8850) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.65) * (0.676 - diff256);
            cutoff64 = 1 - (1 - 0.875) * (0.251 - diff64);
            cutoff16 = 1 - (1 - 0.98905) * (0.21 - diff16);
            cutoff4 = 1 - (1 - 0.9958) * (0.17 - diff4);
        } else if (radius >= 7775) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.6167) * (0.686 - diff256);
            cutoff64 = 1 - (1 - 0.865) * (0.271 - diff64);
            cutoff16 = 1 - (1 - 0.9843) * (0.221 - diff16);
            cutoff4 = 1 - (1 - 0.9937) * (0.161 - diff4);
        } else if (radius >= 6700) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.6) * (0.696 - diff256);
            cutoff64 = 1 - (1 - 0.86) * (0.291 - diff64);
            cutoff16 = 1 - (1 - 0.982) * (0.241 - diff16);
            cutoff4 = 1 - (1 - 0.9926) * (0.171 - diff4);
        } else if (radius >= 5313) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.57) * (0.696 - diff256);
            cutoff64 = 1 - (1 - 0.8511) * (0.326 - diff64);
            cutoff16 = 1 - (1 - 0.9811) * (0.276 - diff16);
            cutoff4 = 1 - (1 - 0.991) * (0.181 - diff4);
        } else if (radius >= 3925) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.555) * (0.696 - diff256);
            cutoff64 = 1 - (1 - 0.8467) * (0.366 - diff64);
            cutoff16 = 1 - (1 - 0.9807) * (0.296 - diff16);
            cutoff4 = 1 - (1 - 0.9902) * (0.186 - diff4);
        } else if (radius >= 2538) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.525) * (0.696 - diff256);
            cutoff64 = 1 - (1 - 0.8422) * (0.396 - diff64);
            cutoff16 = 1 - (1 - 0.9802) * (0.296 - diff16);
            cutoff4 = 1 - (1 - 0.9894) * (0.186 - diff4);
        } else if (radius >= 1150) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.51) * (0.696 - diff256);
            cutoff64 = 1 - (1 - 0.84) * (0.396 - diff64);
            cutoff16 = 1 - (1 - 0.98) * (0.296 - diff16);
            cutoff4 = 1 - (1 - 0.989) * (0.186 - diff4);
        } else if (radius >= 725) {
            initialStep = 256;
            cutoff256 = 1 - (1 - 0.17) * (0.796 - diff256);
            cutoff64 = 1 - (1 - 0.8333) * (0.496 - diff64);
            cutoff16 = 1 - (1 - 0.97) * (0.396 - diff16);
            cutoff4 = 1 - (1 - 0.98) * (0.286 - diff4);
        } else if (radius >= 300) {
            initialStep = 64;
            cutoff64 = 1 - (1 - 0.83) * (0.996 - diff64);
            cutoff16 = 1 - (1 - 0.965) * (0.996 - diff16);
            cutoff4 = 1 - (1 - 0.973) * (0.986 - diff4);
        } else if (radius >= 100) {
            initialStep = 16;
            cutoff16 = 1 - (1 - 0.8) * (1 - diff16);
            cutoff4 = 1 - (1 - 0.965) * (0.99 - diff4);
        } else if (radius >= 95) {
            initialStep = 16;
            cutoff16 = 0.22;
            cutoff4 = 0.8;
        } else if (radius >= 40) {
            initialStep = 16;
            cutoff16 = 0.4;
            cutoff4 = 0.7;
        } else if (radius >= 20) {
            initialStep = 4;
            cutoff4 = 0.4;
        } else if (radius >= 10) {
            initialStep = 4;
            cutoff4 = 0.15;
        } else {
            initialStep = 1;
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

        // TODO: Make this nicer (combine into setCutoffsAndInitialStep, among other things)
        double cutoff256;
        double cutoff64;
        double cutoff16;
        double cutoff4;
        int initialStep;
        setCutoffsAndInitialStep(radius, cutoff256, cutoff64, cutoff16, cutoff4, initialStep);
        // TODO: Come up with a better solution here, this is at least doubling work.
        bool useAlreadyChecked = initialStep > 16;

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
        std::unordered_set<std::pair<int, int>, intPairHash> alreadyChecked;
        std::set<PixelCenterAndPop> topCircles;
        double largestPop = 0;
        std::map<int, std::vector<int>> kernels;

        std::cout << "step: " << step << std::endl;
        mostPopulousCirclesOfGivenRadiusPixelBoundaries(radius, pixelBoundaries, step, topCircles,
                                                        largestPop, cutoff, kernels, desiredPop,
                                                        alreadyChecked, useAlreadyChecked);
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
            const std::set<PixelCenterAndPop> topCirclesCopy{topCircles};
            int currY = -1;
            for (const auto &topCircle : topCirclesCopy) {
                if (topCircle.y != currY) {
                    currY = topCircle.y;
                    kernels.clear();
                }
                PixelBoundaries sectionBoundaries{
                    topCircle.x - step * 4, topCircle.x + step * 4 - 1, topCircle.y - step * 4,
                    topCircle.y + step * 4 - 1};
                mostPopulousCirclesOfGivenRadiusPixelBoundaries(
                    radius, sectionBoundaries, step, topCircles, largestPop, cutoff, kernels,
                    desiredPop, alreadyChecked);
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
        result.lat = lat(topCircles.begin()->y);
        result.lon = lon(topCircles.begin()->x);
        if (topCircles.begin()->pop != largestPop) {
            std::cout << "topCircles[0].pop: " << topCircles.begin()->pop << std::endl;
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
    double popWithinCircle(const double lat, const double lon, const double radius) {
        const int y = latToY(lat);
        const int x = lonToX(lon);
        std::set<PixelCenterAndPop> topCircles;
        double largestPop = 0;
        std::map<int, std::vector<int>> kernels;
        std::unordered_set<std::pair<int, int>, intPairHash> alreadyChecked;
        mostPopulousCirclesOfGivenRadiusPixelBoundaries(radius, PixelBoundaries{x, x, y, y}, 1,
                                                        topCircles, largestPop, 0, kernels,
                                                        100000000000, alreadyChecked, false);
        return largestPop;
    }

    CircleResult
    findNextLargestCity(const double radius, const std::vector<CircleResult> &largestCities,
                        const LatLonBoundaries &boundaries = LatLonBoundaries{-180, 180, 90, -90}) {
        std::cout << "Radius: " << radius << std::endl;
        std::vector<CirclePixelSet> largestCitiesPixels;
        for (const CircleResult &city : largestCities) {
            largestCitiesPixels.emplace_back(*this, city);
        }
        // TODO: Don't look at centers outside these ranges.
        // Turn lat and lon into indices in the pop data.
        PixelBoundaries pixelBoundaries = pixelBoundariesFromLatLonBoundaries(boundaries);
        // TODO: Make this nicer (combine into setCutoffsAndInitialStep, among other things)
        double cutoff256;
        double cutoff64;
        double cutoff16;
        double cutoff4;
        int initialStep;
        setCutoffsAndInitialStep(radius, cutoff256, cutoff64, cutoff16, cutoff4, initialStep);
        // TODO: Come up with a better solution here, this is at least doubling work.
        bool useAlreadyChecked = initialStep > 16;

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
        std::unordered_set<std::pair<int, int>, intPairHash> alreadyChecked;
        std::set<PixelCenterAndPop> topCircles;
        double largestPop = 0;
        std::map<int, std::vector<int>> kernels;

        std::cout << "step: " << step << std::endl;
        mostPopulousCirclesOfGivenRadiusPixelBoundariesForFindingLargestCities(
            radius, pixelBoundaries, step, topCircles, largestPop, cutoff, kernels,
            largestCitiesPixels, alreadyChecked, useAlreadyChecked);
        cutUnderperformingCircles(topCircles, largestPop, cutoff);
        CircleResult result;
        result.radius = radius;
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
            const std::set<PixelCenterAndPop> topCirclesCopy{topCircles};
            int currY = -1;
            for (const auto &topCircle : topCirclesCopy) {
                if (topCircle.y != currY) {
                    currY = topCircle.y;
                    kernels.clear();
                }
                PixelBoundaries sectionBoundaries{
                    topCircle.x - step * 4, topCircle.x + step * 4 - 1, topCircle.y - step * 4,
                    topCircle.y + step * 4 - 1};
                mostPopulousCirclesOfGivenRadiusPixelBoundariesForFindingLargestCities(
                    radius, sectionBoundaries, step, topCircles, largestPop, cutoff, kernels,
                    largestCitiesPixels, alreadyChecked);
            }
            cutUnderperformingCircles(topCircles, largestPop, cutoff);
        }
        // Only the most populous circle remains in the list, since the cutoff for step size 1 is 1.
        result.lat = lat(topCircles.begin()->y);
        result.lon = lon(topCircles.begin()->x);
        if (topCircles.begin()->pop != largestPop) {
            std::cout << "topCircles[0].pop: " << topCircles.begin()->pop << std::endl;
            std::cout << "largestPop: " << largestPop << std::endl;
            throw std::logic_error(
                "Should never get here. largestPop should be first element of topCircles.");
        }
        result.pop = largestPop;
        std::cout << "Population within " << radius << " kilometers of (" << result.lat << ", "
                  << result.lon << "): " << ((long)result.pop) << std::endl;
        return result;
    }

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

    int getNumRows() const { return numRows; }

    int getNumCols() const { return numCols; }

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
            if (radius - lowerBound <= 1) {
                if (upperBound - lowerBound == 1) {
                    // Need to make upperBound not soft
                    radius = upperBound;
                }
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

CirclePixelSet::CirclePixelSet(const RasterDataCircleFinder &circleFinder, CircleResult circle)
    : CirclePixelSet(circleFinder,
                     Pixel{circleFinder.lonToX(circle.lon), circleFinder.latToY(circle.lat)},
                     circle.radius) {}

CirclePixelSet::CirclePixelSet(const RasterDataCircleFinder &circleFinder, Pixel center,
                               double radius)
    : center(center) {
    boundingBox = {circleFinder.getNumCols(), -1, circleFinder.getNumRows(), -1};
    std::vector<int> kernel = circleFinder.makeKernel(center.y, radius);
    for (int kernelRow = 0; kernelRow < kernel.size() / circleFinder.KERNEL_WIDTH; kernelRow++) {
        // TODO: Make rectangle class.
        // Sides of the rectangle
        const int west = kernel[circleFinder.kernelIndex(kernelRow, 0)] + center.x;
        const int east = kernel[circleFinder.kernelIndex(kernelRow, 1)] + center.x;
        const int north = kernel[circleFinder.kernelIndex(kernelRow, 2)] + center.y;
        const int south = kernel[circleFinder.kernelIndex(kernelRow, 3)] + center.y;
        boundingBox.upY = std::min(north + 1, boundingBox.upY);
        boundingBox.downY = std::max(south, boundingBox.downY);
        if (kernel[circleFinder.kernelIndex(kernelRow, 1)] == circleFinder.getNumCols() / 2) {
            // This rectangle encircles the entire latitude
            for (int y = north + 1; y <= south; y++) {
                xRanges.emplace(y, std::pair{0, circleFinder.getNumCols() - 1});
            }
            boundingBox.leftX = std::min(0, boundingBox.leftX);
            boundingBox.rightX = std::max(circleFinder.getNumCols() - 1, boundingBox.rightX);
        } else if (west < -1) {
            // Need to wrap around the antimeridian, creating two rectangles
            for (int y = north + 1; y <= south; y++) {
                xRanges.emplace(y, std::pair{0, east});
                xRanges.emplace(y, std::pair{circleFinder.getNumCols() + west + 1,
                                             circleFinder.getNumCols() - 1});
            }
            boundingBox.leftX = std::min(0, boundingBox.leftX);
            boundingBox.rightX = std::max(circleFinder.getNumCols() - 1, boundingBox.rightX);
        } else if (east >= circleFinder.getNumCols()) {
            // Need to wrap around the antimeridian, creating two rectangles
            for (int y = north + 1; y <= south; y++) {
                xRanges.emplace(y, std::pair{west + 1, circleFinder.getNumCols() - 1});
                xRanges.emplace(y, std::pair{0, east - circleFinder.getNumCols()});
            }
            boundingBox.leftX = std::min(0, boundingBox.leftX);
            boundingBox.rightX = std::max(circleFinder.getNumCols() - 1, boundingBox.rightX);
        } else {
            for (int y = north + 1; y <= south; y++) {
                xRanges.emplace(y, std::pair{west + 1, east});
            }
            boundingBox.leftX = std::min(west + 1, boundingBox.leftX);
            boundingBox.rightX = std::max(east, boundingBox.rightX);
        }
    }
}

bool CirclePixelSet::contains(Pixel p) const {
    if (boundingBox.contains(p)) {
        // Could only find left bound to try and speed things up maybe.
        auto ranges = xRanges.equal_range(p.y);
        for (auto range = ranges.first; range != ranges.second && range != xRanges.end(); range++) {
            if (range->second.first <= p.x && p.x <= range->second.second) {
                return true;
            }
        }
    }
    return false;
}

bool CirclePixelSet::containsY(const int y) const {
    return boundingBox.upY <= y && y <= boundingBox.downY;
}

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
    for (double percent = 100.0; percent > 0.09; percent -= 0.1) {
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

void findPop(double lat, double lon, std::vector<std::string>& results, const std::string& cityName) {
    //std::cout << "Loading population summation table." << std::endl;
    if (USE_2020_DATA) {
        //std::cout << "Using 2020 data." << std::endl;
    } else {
        std::cout << "Using 2015 data." << std::endl;
    }
    std::string sumTableFilename = "popSumTable.bin";
    if (USE_2020_DATA) {
        sumTableFilename = "popSumTable2020.bin";
    }
    static RasterDataCircleFinder popData(sumTableFilename, "placeholder.txt");
    //std::cout << "Loaded population summation table." << std::endl;
    if (!popData.previousResultsConsistent()) {
        return;
    }
    // int radius = 2;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 5;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 10;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 20;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 50;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 100;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 200;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 500;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    int radius = 1000;
    std::stringstream result;
    result << lat << "~" << lon << "~" << ((long long)(popData.popWithinCircle(lat, lon, radius))) << "~" << cityName;
    results.emplace_back(result.str());
    // radius = 2000;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 5000;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    // radius = 10000;
    // std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
}

void findLargestCities(int radius) {
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
    RasterDataCircleFinder popData(sumTableFilename, "placeholder.txt");
    std::cout << "Loaded population summation table." << std::endl;
    if (!popData.previousResultsConsistent()) {
        return;
    }
    // Used by imageManipulation.java so that it knows what text to add, as well as by
    //  percentageCirclesMapMakerTextInJSONOut.cpp so it knows where and how large to draw the
    //  circles
    std::string largestCitiesFilename = "largestCities2015Radius" + std::to_string(radius) + ".txt";
    if (USE_2020_DATA) {
        largestCitiesFilename = "largestCities2020Radius" + std::to_string(radius) + ".txt";
    }
    std::ifstream largestCitiesFile;
    largestCitiesFile.open(largestCitiesFilename);
    // TODO: Make if not exists and also write full pop with overlap to file.
    std::vector<CircleResult> largestCities;
    std::string cityString;
    while (getline(largestCitiesFile, cityString)) {
        std::stringstream citySS(cityString);
        std::string dummyString;
        getline(citySS, dummyString, ' '); // Rank
        getline(citySS, dummyString, ' ');
        double lat = std::stod(dummyString);
        getline(citySS, dummyString, ' ');
        double lon = std::stod(dummyString);
        getline(citySS, dummyString, ' ');
        double pop = std::stod(dummyString);
        largestCities.emplace_back(CircleResult{lat, lon, (double)radius, pop});
    }
    largestCitiesFile.close();
    int rank = 1;
    double smallestPopSoFar = std::numeric_limits<double>::max();
    for (; rank <= largestCities.size(); rank++) {
        std::cout << "\n" << rank << "th largest city:" << std::endl;
        std::cout << "Population within " << largestCities[rank - 1].radius << " km of ("
                  << largestCities[rank - 1].lat << ", " << largestCities[rank - 1].lon
                  << "): " << ((long)largestCities[rank - 1].pop) << std::endl;
        if (largestCities[rank - 1].pop > smallestPopSoFar) {
            std::cout << "WEE WOO WEE WOO INCONSISTENCY" << std::endl;
            std::cout << "WEE WOO WEE WOO INCONSISTENCY" << std::endl;
            std::cout << "WEE WOO WEE WOO INCONSISTENCY" << std::endl;
            return;
        } else {
            smallestPopSoFar = largestCities[rank - 1].pop;
        }
    }
    auto t0 = std::chrono::system_clock::now();
    for (; rank <= 40000; rank++) {
        std::cout << std::endl << "Now finding the " << rank << "th largest city." << std::endl;
        // WARNING: Can't yet touch antimeridian or poles
        CircleResult city =
            popData.findNextLargestCity(radius, largestCities, {-180, 180, 80, -60});
        largestCities.emplace_back(city);
        std::cout << rank << "th largest city:" << std::endl;
        std::cout << "Population within " << city.radius << " km of (" << city.lat << ", "
                  << city.lon << "): " << ((long)city.pop) << std::endl;
        const auto t1 = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << t1 << "\n\n\n";
        const auto tn = std::chrono::system_clock::now();
        std::cout << (std::chrono::duration_cast<std::chrono::seconds>(tn - t0).count() / 60) << "\n\n\n";
        t0 = tn;
        std::fstream largestCitiesOFile;
        largestCitiesOFile.open(largestCitiesFilename);
        largestCitiesOFile.seekg(0, std::ios::end);
        // TODO: Remove the magic number 6
        largestCitiesOFile << rank << " " << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION)
                           << city.lat << " " << city.lon << " " << city.pop << std::setprecision(6)
                           << std::endl;
        largestCitiesOFile.close();
        if (largestCities[rank - 1].pop > smallestPopSoFar) {
            std::cout << "WEE WOO WEE WOO INCONSISTENCY" << std::endl;
            std::cout << "WEE WOO WEE WOO INCONSISTENCY" << std::endl;
            std::cout << "WEE WOO WEE WOO INCONSISTENCY" << std::endl;
            return;
        } else {
            smallestPopSoFar = largestCities[rank - 1].pop;
        }
    }
}

void makeGoogleMapsJavascript(int radius) {
    std::string largestCitiesFilename = "largestCities2015Radius" + std::to_string(radius) + ".txt";
    std::string largestCitiesJavascriptFilename =
        "largestCities2015Radius" + std::to_string(radius) + "Javascript.txt";
    if (USE_2020_DATA) {
        largestCitiesFilename = "largestCities2020Radius" + std::to_string(radius) + ".txt";
        largestCitiesJavascriptFilename =
            "largestCities2020Radius" + std::to_string(radius) + "Javascript.txt";
    }
    std::ifstream largestCitiesFile;
    largestCitiesFile.open(largestCitiesFilename);
    // TODO: Make if not exists and also write full pop with overlap to file.
    std::vector<CircleResult> largestCities;
    std::string cityString;
    while (getline(largestCitiesFile, cityString)) {
        std::stringstream citySS(cityString);
        std::string dummyString;
        getline(citySS, dummyString, ' '); // Rank
        getline(citySS, dummyString, ' ');
        double lat = std::stod(dummyString);
        getline(citySS, dummyString, ' ');
        double lon = std::stod(dummyString);
        getline(citySS, dummyString, ' ');
        double pop = std::stod(dummyString);
        largestCities.emplace_back(CircleResult{lat, lon, (double)radius, pop});
    }
    largestCitiesFile.close();
    int rank = 1;
    std::fstream largestCitiesJavascriptFile;
    largestCitiesJavascriptFile.open(largestCitiesJavascriptFilename);
    largestCitiesJavascriptFile << std::setprecision(DOUBLE_ROUND_TRIP_PRECISION);
    for (const CircleResult &city : largestCities) {
        largestCitiesJavascriptFile << rank++ << ": { center: { lat: " << city.lat
                                    << ", lng: " << city.lon << " } }, " << std::endl;
    }
    largestCitiesJavascriptFile.close();
}

void findLargestCitiesManyRadiuses() {
    for (int radius = 45; radius <= 100; radius += 10) {
        findLargestCities(radius);
    }
}

int main() {
    // makeGoogleMapsJavascript(45);
    // findLargestCitiesManyRadiuses();

    std::vector<std::string> results;
    std::string cityName;

    cityName = "Moscow, Russia:";
    double lat = 55.755833;
    double lon = 37.617222;
    findPop(lat, lon, results, cityName);

    cityName = "Ottawa, Canada";
    lat = 45.424722;
    lon = -75.695;
    findPop(lat, lon, results, cityName);

    cityName = "Beijing, China";
    lat = 39.906667;
    lon = 116.3975;
    findPop(lat, lon, results, cityName);

    cityName = "Washington, D.C., USA";
    lat = 38.904722;
    lon = -77.016389;
    findPop(lat, lon, results, cityName);

    cityName = "Brasilia, Brazil";
    lat = -15.793889;
    lon = -47.882778;
    findPop(lat, lon, results, cityName);

    cityName = "Canberra, Australia";
    lat = -35.293056;
    lon = 149.126944;
    findPop(lat, lon, results, cityName);

    cityName = "New Delhi, India";
    lat = 28.613895;
    lon = 77.209006;
    findPop(lat, lon, results, cityName);

    cityName = "Buenos Aires, Argentina";
    lat = -34.603333;
    lon = -58.381667;
    findPop(lat, lon, results, cityName);

    cityName = "Astana, Kazakhstan";
    lat = 51.166667;
    lon = 71.433333;
    findPop(lat, lon, results, cityName);

    cityName = "Algiers, Algeria";
    lat = 36.753889;
    lon = 3.058889;
    findPop(lat, lon, results, cityName);

    cityName = "Kinshasa, Democratic Republic of the Congo";
    lat = -4.325;
    lon = 15.322222;
    findPop(lat, lon, results, cityName);

    cityName = "Copenhagen, Denmark";
    lat = 55.676111;
    lon = 12.568333;
    findPop(lat, lon, results, cityName);

    lat = 24.633333;
    lon = 46.716667;
    cityName = "Riyadh, Saudi Arabia";
    findPop(lat, lon, results, cityName);

    lat = 19.433333;
    lon = -99.133333;
    cityName = "Mexico City, Mexico";
    findPop(lat, lon, results, cityName);

    lat = -6.175;
    lon = 106.8275;
    cityName = "Jakarta, Indonesia";
    findPop(lat, lon, results, cityName);

    lat = 15.500556;
    lon = 32.56;
    cityName = "Khartoum, Sudan";
    findPop(lat, lon, results, cityName);

    lat = 32.887222;
    lon = 13.191389;
    cityName = "Tripoli, Libya";
    findPop(lat, lon, results, cityName);

    lat = 35.689167;
    lon = 51.388889;
    cityName = "Tehran, Iran";
    findPop(lat, lon, results, cityName);

    lat = 47.920278;
    lon = 106.917222;
    cityName = "Ulaanbaatar, Mongolia";
    findPop(lat, lon, results, cityName);

    lat = -12.06;
    lon = -77.0375;
    cityName = "Lima, Peru";
    findPop(lat, lon, results, cityName);

    lat = 12.11;
    lon = 15.05;
    cityName = "N'Djamena, Chad";
    findPop(lat, lon, results, cityName);

    lat = 13.515;
    lon = 2.1175;
    cityName = "Niamey, Niger";
    findPop(lat, lon, results, cityName);

    lat = -8.838333;
    lon = 13.234444;
    cityName = "Luanda, Angola";
    findPop(lat, lon, results, cityName);

    lat = 12.639167;
    lon = -8.002778;
    cityName = "Bamako, Mali";
    findPop(lat, lon, results, cityName);

    lat = -25.746111;
    lon = 28.188056;
    cityName = "Pretoria, South Africa";
    findPop(lat, lon, results, cityName);

    lat = 4.711111;
    lon = -74.072222;
    cityName = "Bogota, Colombia";
    findPop(lat, lon, results, cityName);

    lat = 9.03;
    lon = 38.74;
    cityName = "Addis Ababa, Ethiopia";
    findPop(lat, lon, results, cityName);

    lat = -19.0475;
    lon = -65.26;
    cityName = "Sucre, Bolivia";
    findPop(lat, lon, results, cityName);

    lat = 18.08581;
    lon = -15.9785;
    cityName = "Nouakchott, Mauritania";
    findPop(lat, lon, results, cityName);

    lat = 30.044444;
    lon = 31.235833;
    cityName = "Cairo, Egypt";
    findPop(lat, lon, results, cityName);

    lat = -6.173056;
    lon = 35.741944;
    cityName = "Dodoma, Tanzania";
    findPop(lat, lon, results, cityName);

    lat = 9.066667;
    lon = 7.483333;
    cityName = "Abuja, Nigeria";
    findPop(lat, lon, results, cityName);

    lat = 10.480556;
    lon = -66.903611;
    cityName = "Caracas, Venezuela";
    findPop(lat, lon, results, cityName);

    lat = 33.693056;
    lon = 73.063889;
    cityName = "Islamabad, Pakistan";
    findPop(lat, lon, results, cityName);

    lat = -22.57;
    lon = 17.083611;
    cityName = "Windhoek, Namibia";
    findPop(lat, lon, results, cityName);

    lat = -25.966667;
    lon = 32.583333;
    cityName = "Maputo, Mozambique";
    findPop(lat, lon, results, cityName);

    lat = 39.93;
    lon = 32.85;
    cityName = "Ankara, Turkey";
    findPop(lat, lon, results, cityName);

    lat = -33.45;
    lon = -70.666667;
    cityName = "Santiago, Chile";
    findPop(lat, lon, results, cityName);

    lat = -15.416667;
    lon = 28.283333;
    cityName = "Lusaka, Zambia";
    findPop(lat, lon, results, cityName);

    lat = 19.7475;
    lon = 96.115;
    cityName = "Naypyidaw, Myanmar";
    findPop(lat, lon, results, cityName);

    lat = 34.525278;
    lon = 69.178333;
    cityName = "Kabul, Afghanistan";
    findPop(lat, lon, results, cityName);

    lat = 4.85;
    lon = 31.6;
    cityName = "Juba, South Sudan";
    findPop(lat, lon, results, cityName);

    lat = 48.856613;
    lon = 2.352222;
    cityName = "Paris, France";
    findPop(lat, lon, results, cityName);

    lat = 2.039167;
    lon = 45.341944;
    cityName = "Mogadishu, Somalia";
    findPop(lat, lon, results, cityName);

    lat = 4.373333;
    lon = 18.562778;
    cityName = "Bangui, Central African Republic";
    findPop(lat, lon, results, cityName);

    lat = 50.45;
    lon = 30.523333;
    cityName = "Kyiv, Ukraine";
    findPop(lat, lon, results, cityName);

    lat = -18.933333;
    lon = 47.516667;
    cityName = "Antananarivo, Madagascar";
    findPop(lat, lon, results, cityName);

    lat = -24.658056;
    lon = 25.912222;
    cityName = "Gaborone, Botswana";
    findPop(lat, lon, results, cityName);

    lat = -1.286389;
    lon = 36.817222;
    cityName = "Nairobi, Kenya";
    findPop(lat, lon, results, cityName);

    lat = 15.348333;
    lon = 44.206389;
    cityName = "Sanaa, Yemen";
    findPop(lat, lon, results, cityName);

    lat = 13.7525;
    lon = 100.494167;
    cityName = "Bangkok, Thailand";
    findPop(lat, lon, results, cityName);

    lat = 40.416667;
    lon = -3.7025;
    cityName = "Madrid, Spain";
    findPop(lat, lon, results, cityName);

    lat = 37.9375;
    lon = 58.38;
    cityName = "Ashgabat, Turkmenistan";
    findPop(lat, lon, results, cityName);

    lat = 3.866667;
    lon = 11.516667;
    cityName = "Yaounde, Cameroon";
    findPop(lat, lon, results, cityName);

    lat = -9.478889;
    lon = 147.149444;
    cityName = "Port Moresby, Papua New Guinea";
    findPop(lat, lon, results, cityName);

    lat = 59.329444;
    lon = 18.068611;
    cityName = "Stockholm, Sweden";
    findPop(lat, lon, results, cityName);

    lat = 41.311111;
    lon = 69.279722;
    cityName = "Tashkent, Uzbekistan";
    findPop(lat, lon, results, cityName);

    lat = 34.020882;
    lon = -6.84165;
    cityName = "Rabat, Morocco";
    findPop(lat, lon, results, cityName);

    lat = 33.315278;
    lon = 44.366111;
    cityName = "Baghdad, Iraq";
    findPop(lat, lon, results, cityName);

    lat = -25.3;
    lon = -57.633333;
    cityName = "Asuncion, Paraguay";
    findPop(lat, lon, results, cityName);

    lat = -17.829167;
    lon = 31.052222;
    cityName = "Harare, Zimbabwe";
    findPop(lat, lon, results, cityName);

    lat = 59.913333;
    lon = 10.738889;
    cityName = "Oslo, Norway";
    findPop(lat, lon, results, cityName);

    lat = 35.689722;
    lon = 139.692222;
    cityName = "Tokyo, Japan";
    findPop(lat, lon, results, cityName);

    lat = 52.52;
    lon = 13.405;
    cityName = "Berlin, Germany";
    findPop(lat, lon, results, cityName);

    lat = -4.269444;
    lon = 15.271389;
    cityName = "Brazzaville Republic of the Congo";
    findPop(lat, lon, results, cityName);

    lat = 60.170833;
    lon = 24.9375;
    cityName = "Helsinki, Finland";
    findPop(lat, lon, results, cityName);

    lat = 21.028333;
    lon = 105.854167;
    cityName = "Hanoi, Vietnam";
    findPop(lat, lon, results, cityName);

    lat = 3.147778;
    lon = 101.695278;
    cityName = "Kuala Lumpur, Malaysia";
    findPop(lat, lon, results, cityName);

    lat = 6.816111;
    lon = -5.274167;
    cityName = "Yamoussoukro, Ivory Coast";
    findPop(lat, lon, results, cityName);

    lat = 52.23;
    lon = 21.011111;
    cityName = "Warsaw, Poland";
    findPop(lat, lon, results, cityName);

    lat = 23.588889;
    cityName = "Muscat, Oman";
    lon = 58.408333;
    findPop(lat, lon, results, cityName);

    lat = 41.893333;
    lon = 12.482778;
    cityName = "Rome, Italy";
    findPop(lat, lon, results, cityName);

    lat = 14.5958;
    lon = 120.9772;
    cityName = "Manila, Philippines";
    findPop(lat, lon, results, cityName);

    lat = -0.22;
    lon = -78.5125;
    cityName = "Quito, Ecuador";
    findPop(lat, lon, results, cityName);

    lat = 12.368611;
    lon = -1.5275;
    cityName = "Ouagadougou, Burkina Faso";
    findPop(lat, lon, results, cityName);

    lat = -41.288889;
    lon = 174.777222;
    cityName = "Wellington, New Zealand";
    findPop(lat, lon, results, cityName);

    lat = 0.390278;
    lon = 9.454167;
    cityName = "Libreville, Gabon";
    findPop(lat, lon, results, cityName);

    lat = 9.509167;
    lon = -13.712222;
    cityName = "Conakry, Guinea";
    findPop(lat, lon, results, cityName);

    lat = 51.507222;
    lon = -0.1275;
    cityName = "London, UK";
    findPop(lat, lon, results, cityName);

    lat = 0.313611;
    lon = 32.581111;
    cityName = "Kampala, Uganda";
    findPop(lat, lon, results, cityName);

    lat = 5.55;
    lon = -0.2;
    cityName = "Accra, Ghana";
    findPop(lat, lon, results, cityName);

    lat = 44.4325;
    lon = 26.103889;
    cityName = "Bucharest, Romania";
    findPop(lat, lon, results, cityName);

    lat = 17.966667;
    lon = 102.6;
    cityName = "Vientiane, Laos";
    findPop(lat, lon, results, cityName);

    lat = 6.805833;
    lon = -58.150833;
    cityName = "Georgetown, Guyana";
    findPop(lat, lon, results, cityName);

    lat = 53.9;
    lon = 27.566667;
    cityName = "Minsk, Belarus";
    findPop(lat, lon, results, cityName);

    lat = 42.874722;
    lon = 74.612222;
    cityName = "Bishkek, Kyrgyzstan";
    findPop(lat, lon, results, cityName);

    lat = 14.692778;
    lon = -17.446667;
    cityName = "Dakar, Senegal";
    findPop(lat, lon, results, cityName);

    lat = 33.513056;
    lon = 36.291944;
    cityName = "Damascus, Syria";
    findPop(lat, lon, results, cityName);

    lat = 11.569444;
    lon = 104.921111;
    cityName = "Phnom Penh, Cambodia";
    findPop(lat, lon, results, cityName);

    lat = -34.883611;
    lon = -56.181944;
    cityName = "Montevideo, Uruguay";
    findPop(lat, lon, results, cityName);

    lat = 5.852222;
    lon = -55.203889;
    cityName = "Paramaribo, Suriname";
    findPop(lat, lon, results, cityName);

    lat = 36.806389;
    lon = 10.181667;
    cityName = "Tunis, Tunisia";
    findPop(lat, lon, results, cityName);

    lat = 23.763889;
    lon = 90.388889;
    cityName = "Dhaka, Bangladesh";
    findPop(lat, lon, results, cityName);

    lat = 27.7172;
    lon = 85.324;
    cityName = "Kathmandu, Nepal";
    findPop(lat, lon, results, cityName);

    lat = 38.536667;
    lon = 68.78;
    cityName = "Dushanbe, Tajikistan";
    findPop(lat, lon, results, cityName);

    lat = 37.984167;
    lon = 23.728056;
    cityName = "Athens, Greece";
    findPop(lat, lon, results, cityName);

    lat = 12.136389;
    lon = -86.251389;
    cityName = "Managua, Nicaragua";
    findPop(lat, lon, results, cityName);

    lat = 15.322778;
    lon = 38.925;
    cityName = "Asmara, Eritrea";
    findPop(lat, lon, results, cityName);

    lat = 39.019444;
    lon = 125.738056;
    cityName = "Pyongyang, North Korea";
    findPop(lat, lon, results, cityName);

    lat = -13.983333;
    lon = 33.783333;
    cityName = "Lilongwe, Malawi";
    findPop(lat, lon, results, cityName);

    lat = 6.497222;
    lon = 2.605;
    cityName = "Porto-Novo, Benin";
    findPop(lat, lon, results, cityName);

    lat = 14.1;
    lon = -87.216667;
    cityName = "Tegucigalpa, Honduras";
    findPop(lat, lon, results, cityName);

    lat = 6.313333;
    lon = -10.801389;
    cityName = "Monrovia, Liberia";
    findPop(lat, lon, results, cityName);

    lat = 42.7;
    lon = 23.33;
    cityName = "Sofia, Bulgaria";
    findPop(lat, lon, results, cityName);

    lat = 23.136667;
    lon = -82.358889;
    cityName = "Havana, Cuba";
    findPop(lat, lon, results, cityName);

    lat = 14.613333;
    lon = -90.535278;
    cityName = "Guatemala City, Guatemala";
    findPop(lat, lon, results, cityName);

    lat = 64.146667;
    lon = -21.94;
    cityName = "Reykjavik, Iceland";
    findPop(lat, lon, results, cityName);

    lat = 37.56;
    lon = 126.99;
    cityName = "Seoul, South Korea";
    findPop(lat, lon, results, cityName);

    lat = 47.4925;
    lon = 19.051389;
    cityName = "Budapest, Hungary";
    findPop(lat, lon, results, cityName);

    lat = 38.725267;
    lon = -9.150019;
    cityName = "Lisbon, Portugal";
    findPop(lat, lon, results, cityName);

    lat = 31.949722;
    lon = 35.932778;
    cityName = "Amman, Jordan";
    findPop(lat, lon, results, cityName);

    lat = 44.817778;
    lon = 20.456944;
    cityName = "Belgrade, Serbia";
    findPop(lat, lon, results, cityName);

    lat = 40.395278;
    lon = 49.882222;
    cityName = "Baku, Azerbaijan";
    findPop(lat, lon, results, cityName);

    lat = 48.2;
    lon = 16.366667;
    cityName = "Vienna, Austria";
    findPop(lat, lon, results, cityName);

    lat = 24.466667;
    lon = 54.366667;
    cityName = "Abu Dhabi, United Arab Emirates";
    findPop(lat, lon, results, cityName);

    lat = 50.0875;
    lon = 14.421389;
    cityName = "Prague, Czechia";
    findPop(lat, lon, results, cityName);

    lat = 8.983333;
    lon = -79.516667;
    cityName = "Panama City, Panama";
    findPop(lat, lon, results, cityName);

    lat = 8.484444;
    lon = -13.234444;
    cityName = "Freetown, Sierra Leone";
    findPop(lat, lon, results, cityName);

    lat = 53.35;
    lon = -6.260278;
    cityName = "Dublin, Ireland";
    findPop(lat, lon, results, cityName);

    lat = 41.7225;
    lon = 44.7925;
    cityName = "Tbilisi, Georgia";
    findPop(lat, lon, results, cityName);

    lat = 6.934444;
    lon = 79.842778;
    cityName = "Colombo, Sri Lanka";
    findPop(lat, lon, results, cityName);

    lat = 54.687222;
    lon = 25.28;
    cityName = "Vilnius, Lithuania";
    findPop(lat, lon, results, cityName);

    lat = 56.948889;
    lon = 24.106389;
    cityName = "Riga, Latvia";
    findPop(lat, lon, results, cityName);

    lat = 6.131944;
    lon = 1.222778;
    cityName = "Lome, Togo";
    findPop(lat, lon, results, cityName);

    lat = 45.816667;
    lon = 15.983333;
    cityName = "Zagreb, Croatia";
    findPop(lat, lon, results, cityName);

    lat = 43.856389;
    lon = 18.413056;
    cityName = "Sarajevo, Bosnia and Herzegovina";
    findPop(lat, lon, results, cityName);

    lat = 9.9325;
    lon = -84.08;
    cityName = "San Jose, Costa Rica";
    findPop(lat, lon, results, cityName);

    lat = 48.143889;
    lon = 17.109722;
    cityName = "Bratislava, Slovakia";
    findPop(lat, lon, results, cityName);

    lat = 18.466667;
    lon = -69.95;
    cityName = "Santo Domingo, Dominican Republic";
    findPop(lat, lon, results, cityName);

    lat = 59.437222;
    lon = 24.745278;
    cityName = "Tallinn, Estonia";
    findPop(lat, lon, results, cityName);

    lat = 52.372778;
    lon = 4.893611;
    cityName = "Amsterdam, Netherlands";
    findPop(lat, lon, results, cityName);

    lat = 46.948056;
    lon = 7.4475;
    cityName = "Bern, Switzerland";
    findPop(lat, lon, results, cityName);

    lat = 27.472222;
    lon = 89.636111;
    cityName = "Thimphu, Bhutan";
    findPop(lat, lon, results, cityName);

    lat = 25.066667;
    lon = 121.516667;
    cityName = "Taipei, Taiwan";
    findPop(lat, lon, results, cityName);

    lat = 11.85;
    lon = -15.566667;
    cityName = "Bissau, Guinea-Bissau";
    findPop(lat, lon, results, cityName);

    lat = 47.022778;
    lon = 28.835278;
    cityName = "Chisinau, Moldova";
    findPop(lat, lon, results, cityName);

    lat = 50.846667;
    lon = 4.3525;
    cityName = "Brussels, Belgium";
    findPop(lat, lon, results, cityName);

    lat = -29.31;
    lon = 27.48;
    cityName = "Maseru, Lesotho";
    findPop(lat, lon, results, cityName);

    lat = 40.181389;
    lon = 44.514444;
    cityName = "Yerevan, Armenia";
    findPop(lat, lon, results, cityName);

    lat = -9.431944;
    lon = 159.955556;
    cityName = "Honiara, Solomon Islands";
    findPop(lat, lon, results, cityName);

    lat = 41.328889;
    lon = 19.817778;
    cityName = "Tirana, Albania";
    findPop(lat, lon, results, cityName);

    lat = 3.752064;
    lon = 8.7737;
    cityName = "Malabo, Equitorial Guinea";
    findPop(lat, lon, results, cityName);

    lat = -3.428333;
    lon = 29.925;
    cityName = "Gitega, Burundi";
    findPop(lat, lon, results, cityName);

    lat = 18.533333;
    lon = -72.333333;
    cityName = "Port-au-Prince, Haiti";
    findPop(lat, lon, results, cityName);

    lat = -1.943889;
    lon = 30.059444;
    cityName = "Kigali, Rwanda";
    findPop(lat, lon, results, cityName);

    lat = 41.996111;
    lon = 21.431667;
    cityName = "Skopje, North Macedonia";
    findPop(lat, lon, results, cityName);

    lat = 11.588333;
    lon = 43.145;
    cityName = "Djibouti City, Djibouti";
    findPop(lat, lon, results, cityName);

    lat = 17.251389;
    lon = -88.766944;
    cityName = "Belmopan, Belize";
    findPop(lat, lon, results, cityName);

    lat = 13.698889;
    lon = -89.191389;
    cityName = "San Salvador, El Salvador";
    findPop(lat, lon, results, cityName);

    lat = 31.778889;
    lon = 35.225556;
    cityName = "Jerusalem, Israel";
    findPop(lat, lon, results, cityName);

    lat = 46.051389;
    lon = 14.506111;
    cityName = "Ljubljana, Slovenia";
    findPop(lat, lon, results, cityName);

    lat = -18.1416;
    lon = 178.4419;
    cityName = "Suva, Fiji";
    findPop(lat, lon, results, cityName);

    lat = 29.369722;
    lon = 47.978333;
    cityName = "Kuwait City, Kuwait";
    findPop(lat, lon, results, cityName);

    lat = -26.416667;
    lon = 31.166667;
    cityName = "Lobamba, Eswatini";
    findPop(lat, lon, results, cityName);

    lat = -8.553611;
    lon = 125.578333;
    cityName = "Dili, East Timor";
    findPop(lat, lon, results, cityName);

    lat = 25.078056;
    lon = -77.338611;
    cityName = "Nassau, Bahamas";
    findPop(lat, lon, results, cityName);

    lat = 42.441286;
    lon = 19.262892;
    cityName = "Podgorica, Montenegro";
    findPop(lat, lon, results, cityName);

    lat = -17.733333;
    lon = 168.316667;
    cityName = "Port Vila, Vanuatu";
    findPop(lat, lon, results, cityName);

    lat = 25.286667;
    lon = 51.533333;
    cityName = "Doha, Qatar";
    findPop(lat, lon, results, cityName);

    lat = 13.453056;
    lon = -16.5775;
    cityName = "Banjul, Gambia";
    findPop(lat, lon, results, cityName);

    lat = 17.971389;
    lon = -76.793056;
    cityName = "Kingston, Jamaica";
    findPop(lat, lon, results, cityName);

    lat = 42.663333;
    lon = 21.162222;
    cityName = "Pristina, Kosovo";
    findPop(lat, lon, results, cityName);

    lat = 33.886944;
    lon = 35.513056;
    cityName = "Beirut, Lebanon";
    findPop(lat, lon, results, cityName);

    lat = 35.1725;
    lon = 33.365;
    cityName = "Nicosia, Cyprus";
    findPop(lat, lon, results, cityName);

    lat = 31.778889;
    lon = 35.225556;
    cityName = "Jerusalem, Palestine";
    findPop(lat, lon, results, cityName);

    lat = 4.890278;
    lon = 114.942222;
    cityName = "Bandar Seri Begawan, Brunei";
    findPop(lat, lon, results, cityName);

    lat = 10.666667;
    lon = -61.516667;
    cityName = "Port of Spain, Trinidad and Tobago";
    findPop(lat, lon, results, cityName);

    lat = 14.918;
    lon = -23.509;
    cityName = "Praia, Cape Verde";
    findPop(lat, lon, results, cityName);

    lat = -13.833333;
    lon = -171.75;
    cityName = "Apia, Samoa";
    findPop(lat, lon, results, cityName);

    lat = 49.611667;
    lon = 6.131944;
    cityName = "Luxembourg City, Luxembourg";
    findPop(lat, lon, results, cityName);

    lat = -20.164444;
    lon = 57.504167;
    cityName = "Port Louis, Mauritius";
    findPop(lat, lon, results, cityName);

    lat = -11.699;
    lon = 43.256;
    cityName = "Moroni, Comoros";
    findPop(lat, lon, results, cityName);

    lat = 0.336111;
    lon = 6.730556;
    cityName = "Sao Tome, Sao Tome and Principe";
    findPop(lat, lon, results, cityName);

    lat = 1.433333;
    lon = 173;
    cityName = "South Tarawa, Kiribati";
    findPop(lat, lon, results, cityName);

    lat = 26.225;
    lon = 50.5775;
    cityName = "Manama, Bahrain";
    findPop(lat, lon, results, cityName);

    lat = 15.301389;
    lon = -61.388333;
    cityName = "Roseau, Dominica";
    findPop(lat, lon, results, cityName);

    lat = -21.133333;
    lon = -175.2;
    cityName = "Nuku'alofa, Tonga";
    findPop(lat, lon, results, cityName);

    lat = 1.283333;
    lon = 103.833333;
    cityName = "Singapore";
    findPop(lat, lon, results, cityName);

    lat = 6.917222;
    lon = 158.158889;
    cityName = "Palikir, Federated States of Micronesia";
    findPop(lat, lon, results, cityName);

    lat = 14.016667;
    lon = -60.983333;
    cityName = "Castries, Saint Lucia";
    findPop(lat, lon, results, cityName);

    lat = 42.5;
    lon = 1.5;
    cityName = "Andorra la Vella, Andorra";
    findPop(lat, lon, results, cityName);

    lat = 7.500556;
    lon = 134.624167;
    cityName = "Ngerulmud, Palau";
    findPop(lat, lon, results, cityName);

    lat = -4.6167;
    lon = 55.45;
    cityName = "Victoria, Seychelles";
    findPop(lat, lon, results, cityName);

    lat = 17.116667;
    lon = -61.85;
    cityName = "St. John's, Antigua and Barbuda";
    findPop(lat, lon, results, cityName);

    lat = 13.0975;
    lon = -59.616667;
    cityName = "Bridgetown, Barbados";
    findPop(lat, lon, results, cityName);

    lat = 13.157778;
    lon = -61.225;
    cityName = "Kingstown, Saint Vincent and the Grenadines";
    findPop(lat, lon, results, cityName);

    lat = 12.05;
    lon = -61.75;
    cityName = "St. George's, Grenada";
    findPop(lat, lon, results, cityName);

    lat = 35.898333;
    lon = 14.5125;
    cityName = "Valletta, Malta";
    findPop(lat, lon, results, cityName);

    lat = 4.175278;
    lon = 73.508889;
    cityName = "Male, Maldives";
    findPop(lat, lon, results, cityName);

    lat = 17.3;
    lon = -62.733333;
    cityName = "Basseterre, Saint Kitts and Nevis";
    findPop(lat, lon, results, cityName);

    lat = 7.1167;
    lon = 171.3667;
    cityName = "Delap-Uliga-Djarrit, Marshall Islands";
    findPop(lat, lon, results, cityName);

    lat = 47.141;
    lon = 9.521;
    cityName = "Vaduz, Liechtenstein";
    findPop(lat, lon, results, cityName);

    lat = 43.9346;
    lon = 12.4473;
    cityName = "City of San Marino, San Marino";
    findPop(lat, lon, results, cityName);

    lat = -8.516667;
    lon = 179.2;
    cityName = "Funafuti, Tuvalu";
    findPop(lat, lon, results, cityName);

    lat = -0.5477;
    lon = 166.920867;
    cityName = "Yaren District, Nauru";
    findPop(lat, lon, results, cityName);

    lat = 43.731111;
    lon = 7.42;
    cityName = "Monaco";
    findPop(lat, lon, results, cityName);

    lat = 41.9025;
    lon = 12.4525;
    cityName = "Vatican City";
    findPop(lat, lon, results, cityName);

    std::ofstream test;
    test.open("bigBoris.txt");
    for (const auto& resulte : results) {
        test << resulte << std::endl;
    }
    test.close();
}