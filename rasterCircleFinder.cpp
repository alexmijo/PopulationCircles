// rasterCircleFinder.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Needs GDAL library
// This is an early version of the program. It's still kinda
// spaghetti code and has code I copied from random places on the internet unattributed (GeoTiff and most of distance())

#include <iostream>
#include <gdal.h>
#include <string>
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>

using namespace std;

const int KERNEL_WIDTH = 4;
// TODO: Make a class that allows the loaded summation table to include this information
// (as well as coordinate to indices and vice versa functions)
const int POP_NUM_ROWS = 2 * 60 * 180;
const int POP_NUM_COLS = 2 * 60 * 360;

class RasterData {
    //coords
    //getNumRows
    //getNumCols
    // Actually, the more I think about it, the more I think I can just move everything in here.
    // Actually, totally nevermind on the above thing, cause of the whole kernel issue. Instead I maybe should just
    //  use a struct.
    // Unless, I actually even put most of main's functionality in here (which honestly doesn't sound like a bad idea).
    //  By that I'm thinking like "most poplous circle in range" function (taking radius and range for center, also
    //  like a parameter about how much to print while doing it and how skippy to be). Yeah the more I think about it
    //  the more that seems like a good idea. I guess I should maybe do 6.031 first and see if that helps me decide.
};

// TODO: Make this less jank (final parameter)
// TODO: This gives very slightly different results than using the Geotiff class's coords function, fix pls
double coords(const int numRows, const int numCols, const int x, const int y, const int latOrLon) {
    if (latOrLon == 0) { // latOrLon is 0 if lattitude is wanted
        return -((((double)y) / ((double)numRows)) * 180.0 - 90.0);
    } else { // Otherwise, longitude is wanted
        return (((double)x) / ((double)numCols)) * 360.0 - 180.0;
    }
}

double distance(const double& lat1, const double& lat2, const double& lon1, const double& lon2)
{
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
        lambda = L + (1 - C) * f * sin_alpha * (sigma + C * sin_sigma * (cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m * cos_2sigma_m)));
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
    double delta_sigma = B * sin_sigma * (cos_2sigma_m + (B / 4.0) * (cos_sigma * (-1 * 2 * cos_2sigma_m_sq) - (B / 6.0) * cos_2sigma_m * (-3 + 4 * sin_sigma * sin_sigma) * (-3 + 4 * cos_2sigma_m_sq)));

    // Distance
    double s = b * A * (sigma - delta_sigma);
    return s;
}

// TODO: make direction an enum
// Only works for up (0) and down (1)
int boundingBoxEdge(const int x, const int y, const double radiusM, const int direction, const int numRows, const int numCols) {
    // Turn x and y (indices in the pop data) into geographic coordinates.
    const double cenLat = coords(numRows, numCols, x, y, 0);
    const double cenLon = coords(numRows, numCols, x, y, 1);

    int edge = y;
    int edgeOfMap = numRows - 1;
    // Start at center
    int incOrDec; // Added to edge to move one pixel in desired direction
    // 0 = up, 1 = down
    if (direction == 1) {
        incOrDec = 1; // Down: increasing
    }
    else {
        incOrDec = -1; // Up: decreasing
    }

    while (edge >= 0 && edge <= edgeOfMap) {
        double lat, lon;
        lat = coords(numRows, numCols, x, edge, 0);
        lon = coords(numRows, numCols, x, edge, 1);
        const double distanceFromCen = distance(lat, cenLat, lon, cenLon);
        if (distanceFromCen > radiusM) {
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
inline int kernelIndex(const int i, const int j) {
    return KERNEL_WIDTH * i + j;
}

// Makes the kernel for a specific lattitude. Assigns the kernel's length to the reference parameter kernelLength.
// Right now it just does one summation table rectangle for each row of the kernel
// TODO: make rectangles stretch across multiple rows where applicable
int* makeKernel(const int cenX, const int cenY, const double radiusM, const int numRows, const int numCols, int& kernelLength) {
    const double cenLat = coords(numRows, numCols, cenX, cenY, 0);
    const double cenLon = coords(numRows, numCols, cenX, cenY, 1);
    const int northEdge = boundingBoxEdge(cenX, cenY, radiusM, 0, numRows, numCols);
    const int southEdge = boundingBoxEdge(cenX, cenY, radiusM, 1, numRows, numCols);
    const int maxPossibleLength = southEdge - northEdge + 1;
    // A faster way of doing a 2D array of dimensions maxPossibleSize x KERNEL_WIDTH
    // Each row consists of: {westX, eastX, northY, southY} describing a summation table rectangle (so KERNEL_WIDTH
    // must be 4) relative to cenX and cenY
    int* tempKernel = new int[maxPossibleLength * KERNEL_WIDTH]; // Temp just cause we don't know length of real kernel yet

    int kernelRow = 0; // First index into kernel
    int y = northEdge;
    int horizontalOffset = 0; // From the verticle center line of the kernel
    while (y <= southEdge) {
        double lat = coords(numRows, numCols, cenX + horizontalOffset, y, 0);
        double lon = coords(numRows, numCols, cenX + horizontalOffset, y, 1);
        if (distance(lat, cenLat, lon, cenLon) > radiusM) {
            if (horizontalOffset == 0) {
                cout << "Something went wrong!1" << endl; // TODO: Probably some better way of logging/displaying this error
            }
            horizontalOffset--;
        }
        else {
            // See how much farther out we can go at this y level and still be in the circle
            while (distance(lat, cenLat, lon, cenLon) <= radiusM) {
                horizontalOffset++;
                if (horizontalOffset > numCols / 2) { // This rectangle wraps around the world
                    break;
                }
                lat = coords(numRows, numCols, cenX + horizontalOffset, y, 0);
                lon = coords(numRows, numCols, cenX + horizontalOffset, y, 1);
            }
            horizontalOffset--; // horizontalOffset is now maximally far (after this decrement)

            // Start a new rectangle at y
            tempKernel[kernelIndex(kernelRow, 0)] = -horizontalOffset - 1;
            tempKernel[kernelIndex(kernelRow, 1)] = horizontalOffset;
            tempKernel[kernelIndex(kernelRow, 2)] = y - cenY - 1;
            if (y == southEdge) { // If we just started a new rectangle at southEdge, it must end there too
                tempKernel[kernelIndex(kernelRow, 3)] = y - cenY;
                break;
            }
            else { // Find where this new rectangle ends
                while (true) {
                    y++;
                    if (y > southEdge) {
                        tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1; // Rectangle can't extend below south edge
                        break;
                    }

                    // Check if the circle has widened
                    if (horizontalOffset < numCols / 2) { // Rectangles that wrap around the whole world can't widen any more
                        lat = coords(numRows, numCols, cenX + horizontalOffset + 1, y, 0);
                        lon = coords(numRows, numCols, cenX + horizontalOffset + 1, y, 1);
                        if (distance(lat, cenLat, lon, cenLon) <= radiusM) {
                            tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1; // The circle has widened; the rectangle is done
                            kernelRow++;
                            break;
                        }
                    }

                    lat = coords(numRows, numCols, cenX + horizontalOffset, y, 0);
                    lon = coords(numRows, numCols, cenX + horizontalOffset, y, 1);
                    if (distance(lat, cenLat, lon, cenLon) > radiusM) {
                        tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1; // The y value can no longer be in the rectangle
                        kernelRow++;
                        break;
                    }
                }
            }
        }
    }

    kernelLength = kernelRow + 1; // Visible to the caller
    // Now that we know the length of kernel, we can construct it using the data in tempKernel and then return it
    int* kernel = new int[kernelLength * KERNEL_WIDTH];
    for (kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
        for (int kernelCol = 0; kernelCol < KERNEL_WIDTH; kernelCol++) {
            kernel[kernelIndex(kernelRow, kernelCol)] = tempKernel[kernelIndex(kernelRow, kernelCol)];
        }
    }
    delete[] tempKernel;
    return kernel;
}

// Returns the total population in the specified rectangle. West and north are 1 pixel outside of the rectangle,
// east and south are 1 pixel inside of the rectangle.
// TODO: Make this not use conditional logic by just padding the north and west sides of popSumTable with 0s
// TODO: If this makes the program slow, see if making it inline helps
double popWithinRectangle(const int west, const int east, const int north, const int south, double** popSumTable) {
    if (west == -1) {
        if (north == -1) {
            return popSumTable[south][east]; // West side of rectangle on the antimeridian and north side of rectangle on the north pole
        }
        return popSumTable[south][east] - popSumTable[north][east]; // West side of rectangle on the antimeridian
    } else if (north == -1) { // North side of rectangle on the north pole
        return popSumTable[south][east] - popSumTable[south][west];
    }
    return popSumTable[south][east] - popSumTable[north][east] - popSumTable[south][west] + popSumTable[north][west];
}

// Returns the population of a circle with a radius of <radiusKm> kilometers centered at the given latittude <cenLat>
// and longitude <cenLon>.
// Uses pop as the population data. dataTiff must be the Geotiff that pop is from. Need both to avoid constantly making
// new arrays I think. Could alternatively modify the Geotiff class to just include the pointer to the array I think.
// pop must have first dimension corresponding to longitude, second corresponding to reverse lattitude
// TODO: Make a coordinates class
// TODO: figure out where to put const around pop type
double popWithinKernel(const int cenX, const int cenY, int* kernel, const int kernelLength, double** popSumTable, const int numCols) {
    double totalPop = 0;

    for (int kernelRow = 0; kernelRow < kernelLength; kernelRow++) {
        // Sides of the rectangle
        const int west = kernel[kernelIndex(kernelRow, 0)] + cenX;
        const int east = kernel[kernelIndex(kernelRow, 1)] + cenX;
        const int north = kernel[kernelIndex(kernelRow, 2)] + cenY;
        const int south = kernel[kernelIndex(kernelRow, 3)] + cenY;

        if (kernel[kernelIndex(kernelRow, 1)] == numCols / 2) { // This rectangle encircles the entire latitude
            totalPop += popWithinRectangle(-1, numCols - 1, north, south, popSumTable);
        } else if (west < -1) { // Need to wrap around the antimeridian, creating two rectangles
            totalPop += popWithinRectangle(-1, east, north, south, popSumTable);
            totalPop += popWithinRectangle(numCols + west, numCols - 1, north, south, popSumTable);
        } else if (east >= numCols) { // Need to wrap around the antimeridian, creating two rectangles
            totalPop += popWithinRectangle(west, numCols - 1, north, south, popSumTable);
            totalPop += popWithinRectangle(-1, east - numCols, north, south, popSumTable);
        } else {
            totalPop += popWithinRectangle(west, east, north, south, popSumTable);
        }
    }

    return totalPop;
}

int main() {
    //-------------------------------Parameters---------------------------------------------
    // TODO: Add the ability to use multiple radiuses
    double radiusKm = 500;

    // TODO: Make it so that the program works even if populationMode is false
    const bool populationMode = true; // TODO: make this an enum. false means GDP PPP mode
    const bool smallestPopMode = false;
    
    const int printMod = 73; // Print lattitude when the Y index mod this is 0

    double upLat = 65;
    double downLat = -90;
    double leftLon = -180;
    double rightLon = 180;

    double smallest = 100000000000;
    double largest = 0;
    const int smallStep = 1; //1
    const int mediumStep = 4; //4
    const int largeStep = 16; //16
    const int xLStep = 64; //64
    const int xXLStep = 256; //256
    //-------------------------------Parameters-end-----------------------------------------

    fstream dataSumTableFile;
    dataSumTableFile.open("popSumTable.bin", ios::in | ios::binary);
    int numRows;
    int numCols;
    dataSumTableFile.read(reinterpret_cast<char *>(&numRows), sizeof(int));
    dataSumTableFile.read(reinterpret_cast<char *>(&numCols), sizeof(int));
    if (numRows != POP_NUM_ROWS || numCols != POP_NUM_COLS) {
        std::cout << "Pop summation table isn't expected dimensions" << std::endl;
        return 1;
    }

    double geoTransform[6];
    for (int i = 0; i < 6; i++) {
        dataSumTableFile.read(reinterpret_cast<char *>(&geoTransform[i]), sizeof(double));
    }

    double** dataSumTable = new double*[numRows];
    for(int r = 0; r < numRows; r++) {
        dataSumTable[r] = new double[numCols];
        for (int c = 0; c < numCols; c++) {
            // TODO: See if it's faster to read in entire rows at a time
            dataSumTableFile.read(reinterpret_cast<char *>(&dataSumTable[r][c]), sizeof(double));
        }
    }
    dataSumTableFile.close();

    cout << "Loaded " << (populationMode ? "population" : "GDP PPP") << " summation table." << endl;

    double radiusM = radiusKm * 1000;
    int step = smallStep;
    // Turn lat and lon into indices in the pop data.
    const int leftX = ((leftLon + 180.0) / 360.0) * numCols;
    const int rightX = ((rightLon + 180.0) / 360.0) * numCols;
    const int upY = ((-upLat + 90.0) / 180.0) * numRows;
    const int downY = ((-downLat + 90.0) / 180.0) * numRows;

    for (int cenY = upY; cenY < downY; cenY += step) {
        if (cenY % printMod == 0) {
            cout << "Current lattitude: " << coords(numRows, numCols, 100, cenY, 0) << endl;
        }

        int kernelLength;
        int* kernel = makeKernel(1000, cenY, radiusM, numRows, numCols, kernelLength); // Initializes kernelLength

        for (int cenX = leftX; cenX <= rightX; cenX += step) {
            double popWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength, dataSumTable, numCols);
            if (smallestPopMode) {
                if (popWithinNKilometers < smallest) {
                    cout << (populationMode ? "Population" : "GDP PPP") << " within " << radiusKm << " kilometers of (" \
                        << (coords(numRows, numCols, cenX, cenY, 0)) << ", " << coords(numRows, numCols, cenX, cenY, 1) << "): " \
                        << ((long long)popWithinNKilometers) << endl;
                    smallest = popWithinNKilometers;
                }

                if (popWithinNKilometers < smallest * 1.2) {
                    if (step > smallStep) {
                        cenX -= step;
                    }
                    step = smallStep;
                }
                else if (popWithinNKilometers < smallest * 1.4) {
                    if (step > mediumStep) {
                        cenX -= step;
                    }
                    step = mediumStep;
                }
                else if (popWithinNKilometers < smallest * 1.6) {
                    if (step > largeStep) {
                        cenX -= step;
                    }
                    step = largeStep;
                }
                else if (popWithinNKilometers < smallest * 1.8) {
                    if (step > xLStep) {
                        cenX -= step;
                    }
                    step = xLStep;
                }
                else {
                    step = xXLStep;
                }
            }
            else {
                if (popWithinNKilometers > largest) {
                    cout << (populationMode ? "Population" : "GDP PPP") << " within " << radiusKm << " kilometers of (" \
                        << coords(numRows, numCols, cenX, cenY, 0) << ", " << coords(numRows, numCols, cenX, cenY, 1) << "): " \
                        << ((long long)popWithinNKilometers) << endl;
                    largest = popWithinNKilometers;
                }

                if (popWithinNKilometers > largest * 0.8) {
                    if (step > smallStep) {
                        cenX -= step;
                    }
                    step = smallStep;
                }
                else if (popWithinNKilometers > largest * 0.6) {
                    if (step > mediumStep) {
                        cenX -= step;
                    }
                    step = mediumStep;
                }
                else if (popWithinNKilometers > largest * 0.4) {
                    if (step > largeStep) {
                        cenX -= step;
                    }
                    step = largeStep;
                }
                else if (popWithinNKilometers > largest * 0.2) {
                    if (step > xLStep) {
                        cenX -= step;
                    }
                    step = xLStep;
                }
                else {
                    step = xXLStep;
                }
            }
        }
        delete[] kernel;
        step = smallStep;
    }
}