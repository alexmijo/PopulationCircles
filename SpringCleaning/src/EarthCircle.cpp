#pragma once
#include "RasterData.cpp"
#include "Utility.cpp"
#include <fstream>
#include <vector>

// The following commented out code is from RasterData.cpp

// struct Pixel {
//     // x from left to right, y from bottom to top
//     int x, y;
// };
// Returns the northernmost or southernmost row containing any pixel within <radius> km of the
    //  given row (<y>).
    // int boundingBoxEdge(const int y, const double radius, const bool isNorth) const {
    //     const double cenLat = lat(y);

    //     // Start at center
    //     int edge = y;
    //     int edgeOfMap = kNumRows - 1;
    //     // Added to edge to move one pixel in desired direction
    //     int incOrDec = isNorth ? -1 : 1;

    //     while (edge >= 0 && edge <= edgeOfMap) {
    //         double currLat;
    //         currLat = lat(edge);
    //         // Since this function only works for north and south, longitude never changes
    //         const double distanceFromCen = distance({currLat, 0}, {cenLat, 0});
    //         if (distanceFromCen > radius) {
    //             edge -= incOrDec; // Went too far, walk it back
    //             break;
    //         } else {
    //             edge += incOrDec;
    //         }
    //     }
    //     if (edge < 0) {
    //         edge = 0;
    //     } else if (edge > edgeOfMap) {
    //         edge = edgeOfMap;
    //     }
    //     return edge;
    // }
// // Defines a rectangular region of the raster data. Ranges are inclusive.
// class PixelRectangle {
//   public:
//     int west, east;
//     int south, north;

//     // TODO: Is Pixel small enough that pass by value is faster?
//     bool contains(const Pixel &p) const {
//         return west <= p.x && p.x <= east && north <= p.y && p.y <= south;
//     }
// };

// The following commented out code is from Utility.cpp

//

class EarthCircle {
  public:
    EarthCircle() = default;

    EarthCircle(const Pixel &center, const double radius, const RasterData &rasterData)
        : m_RasterData(rasterData), m_Radius(radius), m_Center(center) {
        makeRectangles();
    }

    bool contains(const Pixel &p) const {
        return distance(m_RasterData.pixelToLocation(m_Center), m_RasterData.pixelToLocation(p)) <=
               m_Radius;
    }

  private: 
    // TODO: See if this double counts part of 1 column when wrapping around all
    // TODO: Get rid of the cenX parameter.
    void makeRectangles() {
        const double cenLon = m_RasterData.lon(m_Center.x);
        const double cenLat = m_RasterData.lat(m_Center.y);
        const int northEdge = m_RasterData.boundingBoxEdge(m_Center.y, m_Radius, true);
        const int southEdge = m_RasterData.boundingBoxEdge(m_Center.y, m_Radius, false);
        const int maxPossibleLength = southEdge - northEdge + 1;
        // A faster way of doing a 2D array of dimensions maxPossibleSize x KERNEL_WIDTH
        // Each row consists of: {westX, eastX, northY, southY} describing a summation table
        //  rectangle (so KERNEL_WIDTH must be 4) relative to cenX and cenY
        std::vector<int> kernel;

        int y = northEdge;
        int horizontalOffset = 0; // From the verticle center line of the kernel
        while (y <= southEdge) {
            double currLon = m_RasterData.lon(m_Center.x + horizontalOffset);
            double currLat = m_RasterData.lat(y);
            if (distance({currLat, currLon}, {cenLat, cenLon}) > m_Radius) {
                if (horizontalOffset == 0) {
                    throw std::logic_error("Kernel making failed. Should never get here.");
                }
                horizontalOffset--;
            } else {
                // See how much farther out we can go at this y level and still be in the circle
                while (distance({currLat, currLon}, {cenLat, cenLon}) <= m_Radius) {
                    horizontalOffset++;
                    if (horizontalOffset > kNumCols / 2) {
                        // This rectangle wraps around the world
                        break;
                    }
                    currLon = m_RasterData.lon(m_Center.x + horizontalOffset);
                    currLat = m_RasterData.lat(y);
                }
                horizontalOffset--; // horizontalOffset is now maximally far (after this decrement)

                // Start a new rectangle at y
                const int west = -horizontalOffset - 1;
                const int east = horizontalOffset;
                const int north = y - m_Center.y - 1;
                int south;
                if (y == southEdge) {
                    // If we just started a new rectangle at southEdge, it must end there too
                    south = y - m_Center.y;
                    break; // huh?
                } else {
                    // Find where this new rectangle ends
                    while (true) {
                        y++;
                        if (y > southEdge) {
                            // Rectangle can't extend below south edge
                            south = y - m_Center.y - 1;
                            break;
                        }

                        // Check if the circle has widened
                        if (horizontalOffset < kNumCols / 2) {
                            // Rectangles that wrap around the whole world can't widen any more
                            currLon = m_RasterData.lon(m_Center.x + horizontalOffset + 1);
                            currLat = m_RasterData.lat(y);
                            if (distance({currLat, currLon}, {cenLat, cenLon}) <= m_Radius) {
                                // The circle has widened; the rectangle is done
                                kernel.emplace_back(y - m_Center.y - 1);
                                break;
                            }
                        }

                        currLon = m_RasterData.lon(m_Center.x + horizontalOffset);
                        currLat = m_RasterData.lat(y);
                        if (distance({currLat, currLon}, {cenLat, cenLon}) > m_Radius) {
                            // The y value can no longer be in the rectangle
                            kernel.emplace_back(y - m_Center.y - 1);
                            break;
                        }
                    }
                }
            }
        }
    }

    const Pixel m_Center;
    // TODO: double?
    // km
    const double m_Radius;
    std::vector<PixelRectangle> m_Rectangles;
    const RasterData &m_RasterData;
};