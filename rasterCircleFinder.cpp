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
constexpr bool USE_2020_DATA = false;
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
                            // TODO: Deal with antimeridian and poles
                            maxEWExtent = maxEastWestExtent(kernel);
                            boundingBox = PixelBoundaries{
                                cenX - maxEWExtent, cenX + maxEWExtent,
                                cenY + kernel[kernelIndex(0, 2)] + 1,
                                cenY + kernel[kernelIndex(kernel.size() / KERNEL_WIDTH - 1, 3)]};
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

void findPop(double lat, double lon) {
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
    int radius = 2;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 5;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 10;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 20;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 50;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 100;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 200;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 500;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 1000;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 2000;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 5000;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
    radius = 10000;
    std::cout << "Pop within" << radius << " km of (" << lat << ", " << lon << "): " << ((long long)(popData.popWithinCircle(lat, lon, radius))) << std::endl;
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
    findLargestCitiesManyRadiuses();
    // double lat = 41.8816344;
    // double lon = -87.6288451;
    // findPop(lat, lon);
}