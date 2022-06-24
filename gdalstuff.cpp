// gdalstuff.cpp : This file contains the 'main' function. Program execution begins and ends there.
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

using namespace std;

const int KERNEL_WIDTH = 4;

class Geotiff {

private: // NOTE: "private" keyword is redundant here.  
         // we place it here for emphasis. Because these
         // variables are declared outside of "public", 
         // they are private. 

    const char* filename;        // name of Geotiff
    GDALDataset* geotiffDataset; // Geotiff GDAL datset object. 
    double geotransform[6];      // 6-element geotranform array.
    int dimensions[3];           // X,Y, and Z dimensions. 
    int NROWS, NCOLS, NLEVELS;     // dimensions of data in Geotiff. 

public:

    // define constructor function to instantiate object
    // of this Geotiff class. 
    Geotiff(const char* tiffname) {
        filename = tiffname;
        GDALAllRegister();

        // set pointer to Geotiff dataset as class member.  
        geotiffDataset = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);

        // set the dimensions of the Geotiff 
        NROWS = GDALGetRasterYSize(geotiffDataset);
        NCOLS = GDALGetRasterXSize(geotiffDataset);
        NLEVELS = GDALGetRasterCount(geotiffDataset);

    }

    // define destructor function to close dataset, 
    // for when object goes out of scope or is removed
    // from memory. 
    ~Geotiff() {
        // close the Geotiff dataset, free memory for array.
        GDALClose(geotiffDataset);
        GDALDestroyDriverManager();
    }

    const char* GetFileName() {
        /*
         * function GetFileName()
         * This function returns the filename of the Geotiff.
         */
        return filename;
    }

    const char* GetProjection() {
        /* function const char* GetProjection():
         *  This function returns a character array (string)
         *  for the projection of the geotiff file. Note that
         *  the "->" notation is used. This is because the
         *  "geotiffDataset" class variable is a pointer
         *  to an object or structure, and not the object
         *  itself, so the "." dot notation is not used.
         */
        return geotiffDataset->GetProjectionRef();
    }

    double* GetGeoTransform() {
        /*
         * function double *GetGeoTransform()
         *  This function returns a pointer to a double that
         *  is the first element of a 6 element array that holds
         *  the geotransform of the geotiff.
         */
        geotiffDataset->GetGeoTransform(geotransform);
        return geotransform;
    }

    double coords(int x, int y, int idx) {
        double* geoTransform = GetGeoTransform();
        double lon = geoTransform[0] + x * geoTransform[1] + y * geoTransform[2];
        double lat = geoTransform[3] + x * geoTransform[4] + y * geoTransform[5];
        double coordinates[2];
        coordinates[0] = lat;
        coordinates[1] = lon;
        return coordinates[idx];
    };

    double GetNoDataValue() {
        /*
         * function GetNoDataValue():
         *  This function returns the NoDataValue for the Geotiff dataset.
         *  Returns the NoData as a double.
         */
        return (double)geotiffDataset->GetRasterBand(1)->GetNoDataValue();
    }

    int* GetDimensions() {
        /*
         * int *GetDimensions():
         *
         *  This function returns a pointer to an array of 3 integers
         *  holding the dimensions of the Geotiff. The array holds the
         *  dimensions in the following order:
         *   (1) number of columns (x size)
         *   (2) number of rows (y size)
         *   (3) number of bands (number of bands, z dimension)
         */
        dimensions[0] = NROWS;
        dimensions[1] = NCOLS;
        dimensions[2] = NLEVELS;
        return dimensions;
    }

    double** GetRasterBand(int z) {

        /*
         * function float** GetRasterBand(int z):
         * This function reads a band from a geotiff at a
         * specified vertical level (z value, 1 ...
         * n bands). To this end, the Geotiff's GDAL
         * data type is passed to a switch statement,
         * and the template function GetArray2D (see below)
         * is called with the appropriate C++ data type.
         * The GetArray2D function uses the passed-in C++
         * data type to properly read the band data from
         * the Geotiff, cast the data to float**, and return
         * it to this function. This function returns that
         * float** pointer.
         */

        double** bandLayer = new double* [NROWS];
        switch (GDALGetRasterDataType(geotiffDataset->GetRasterBand(z))) {
        case 0:
            return NULL; // GDT_Unknown, or unknown data type.
        case 1:
            // GDAL GDT_Byte (-128 to 127) - unsigned  char
            return GetArray2D<unsigned char>(z, bandLayer);
        case 2:
            // GDAL GDT_UInt16 - short
            return GetArray2D<unsigned short>(z, bandLayer);
        case 3:
            // GDT_Int16
            return GetArray2D<short>(z, bandLayer);
        case 4:
            // GDT_UInt32
            return GetArray2D<unsigned int>(z, bandLayer);
        case 5:
            // GDT_Int32
            return GetArray2D<int>(z, bandLayer);
        case 6:
            // GDT_Float32
            return GetArray2D<float>(z, bandLayer);
        case 7:
            // GDT_Float64
            return GetArray2D<double>(z, bandLayer);
        default:
            break;
        }
        return NULL;
    }

    template<typename T>
    double** GetArray2D(int layerIndex, double** bandLayer) {

        /*
         * function float** GetArray2D(int layerIndex):
         * This function returns a pointer (to a pointer)
         * for a float array that holds the band (array)
         * data from the geotiff, for a specified layer
         * index layerIndex (1,2,3... for GDAL, for Geotiffs
         * with more than one band or data layer, 3D that is).
         *
         * Note this is a template function that is meant
         * to take in a valid C++ data type (i.e. char,
         * short, int, float), for the Geotiff in question
         * such that the Geotiff band data may be properly
         * read-in as numbers. Then, this function casts
         * the data to a float data type automatically.
         */

         // get the raster data type (ENUM integer 1-12, 
         // see GDAL C/C++ documentation for more details)        
        GDALDataType bandType = GDALGetRasterDataType(
            geotiffDataset->GetRasterBand(layerIndex));

        // get number of bytes per pixel in Geotiff
        int nbytes = GDALGetDataTypeSizeBytes(bandType);

        // allocate pointer to memory block for one row (scanline) 
        // in 2D Geotiff array.  
        T* rowBuff = (T*)CPLMalloc(nbytes * NCOLS);

        for (int row = 0; row < NROWS; row++) {     // iterate through rows

          // read the scanline into the dynamically allocated row-buffer       
            CPLErr e = geotiffDataset->GetRasterBand(layerIndex)->RasterIO(
                GF_Read, 0, row, NCOLS, 1, rowBuff, NCOLS, 1, bandType, 0, 0);
            if (!(e == 0)) {
                cout << "Warning: Unable to read scanline in Geotiff!" << endl;
                exit(1);
            }

            bandLayer[row] = new double[NCOLS];
            for (int col = 0; col < NCOLS; col++) { // iterate through columns
                bandLayer[row][col] = (double)rowBuff[col];
            }
        }
        CPLFree(rowBuff);
        return bandLayer;
    }

};

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
int boundingBoxEdge(const int x, const int y, const double radiusM, const int direction, Geotiff& popTiff) {
    const int numRows = popTiff.GetDimensions()[0];
    // Turn x and y (indices in the pop data) into geographic coordinates.
    const double cenLat = popTiff.coords(x, y, 0);
    const double cenLon = popTiff.coords(x, y, 1);

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
        lat = popTiff.coords(x, edge, 0);
        lon = popTiff.coords(x, edge, 1);
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
int* makeKernel(const int cenX, const int cenY, const double radiusM, Geotiff& popTiff, int& kernelLength) {
    const int numCols = popTiff.GetDimensions()[1];
    const double cenLat = popTiff.coords(cenX, cenY, 0);
    const double cenLon = popTiff.coords(cenX, cenY, 1);
    const int northEdge = boundingBoxEdge(cenX, cenY, radiusM, 0, popTiff);
    const int southEdge = boundingBoxEdge(cenX, cenY, radiusM, 1, popTiff);
    const int maxPossibleLength = southEdge - northEdge + 1;
    // A faster way of doing a 2D array of dimensions maxPossibleSize x KERNEL_WIDTH
    // Each row consists of: {westX, eastX, northY, southY} describing a summation table rectangle (so KERNEL_WIDTH
    // must be 4) relative to cenX and cenY
    int* tempKernel = new int[maxPossibleLength * KERNEL_WIDTH]; // Temp just cause we don't know length of real kernel yet

    int kernelRow = 0; // First index into kernel
    int y = northEdge;
    int horizontalOffset = 0; // From the verticle center line of the kernel
    while (y <= southEdge) {
        double lat = popTiff.coords(cenX + horizontalOffset, y, 0);
        double lon = popTiff.coords(cenX + horizontalOffset, y, 1);
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
                lat = popTiff.coords(cenX + horizontalOffset, y, 0);
                lon = popTiff.coords(cenX + horizontalOffset, y, 1);
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
                        lat = popTiff.coords(cenX + horizontalOffset + 1, y, 0);
                        lon = popTiff.coords(cenX + horizontalOffset + 1, y, 1);
                        if (distance(lat, cenLat, lon, cenLon) <= radiusM) {
                            tempKernel[kernelIndex(kernelRow, 3)] = y - cenY - 1; // The circle has widened; the rectangle is done
                            kernelRow++;
                            break;
                        }
                    }

                    lat = popTiff.coords(cenX + horizontalOffset, y, 0);
                    lon = popTiff.coords(cenX + horizontalOffset, y, 1);
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

// Turn the passed in pop table into a summation table (mutates it)
void turnIntoSummationTable(double** pop, Geotiff& popTiff) {
    const int numRows = popTiff.GetDimensions()[0];
    const int numCols = popTiff.GetDimensions()[1];
    for (int x = 0; x < numCols; x++) {
        for (int y = 0; y < numRows; y++) {
            if (pop[y][x] < 0) {
                pop[y][x] = 0.0;
            }
            if (x == 0 && y == 0) {
                continue;
            }
            else if (x == 0) {
                pop[y][x] += pop[y - 1][x];
            }
            else if (y == 0) {
                pop[y][x] += pop[y][x - 1];
            }
            else {
                pop[y][x] += pop[y][x - 1] + pop[y - 1][x] - pop[y - 1][x - 1];
            }
        }
    }
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
// Uses pop as the population data. popTiff must be the Geotiff that pop is from. Need both to avoid constantly making
// new arrays I think. Could alternatively modify the Geotiff class to just include the pointer to the array I think.
// pop must have first dimension corresponding to longitude, second corresponding to reverse lattitude
// TODO: Make a coordinates class
// TODO: figure out where to put const around pop type
double popWithinKernel(const int cenX, const int cenY, int* kernel, const int kernelLength, double** popSumTable, Geotiff& popTiff) {
    const int numCols = popTiff.GetDimensions()[1];
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
    // Load population data
    const string popDataFilename = "GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif";
    const string gdpDataFilename = "/mnt/c/Users/Administrator/source/repos/gdalstuff/gdalstuff/gdpPPPdata.tif";
    Geotiff gdpTiff(gdpDataFilename.c_str());
    Geotiff popTiff(popDataFilename.c_str());
    const int numRows = popTiff.GetDimensions()[0];
    const int numCols = popTiff.GetDimensions()[1];
    if (numRows != gdpTiff.GetDimensions()[0] || numCols != gdpTiff.GetDimensions()[1]) {
        cout << "Population tiff and GDP tiff aren't same dimensions." << endl; // TODO: Probably some better way to do this
    }
    // The population data as a 2d array. First dimension is reverse lattitude, second dimension is longitude.
    double** pop = popTiff.GetRasterBand(1);

    cout << "Loaded population tiff and GDP tiff." << endl;

    turnIntoSummationTable(pop, popTiff); // Mutates pop

    cout << "Constructed population summation table." << endl;

    //-------------------------------Parameters---------------------------------------------
    // TODO: Add the ability to use multiple radiuses
    double radiusKm = 15000;

    const bool smallestPopMode = true;
    
    const int printMod = 11; // Print lattitude when the Y index mod this is 0

    double upLat = 90;
    double downLat = -90;
    double leftLon = -180;
    double rightLon = 180;

    double smallest = 100000000000;
    double largest = 0;
    const int smallStep = 16; //1
    const int mediumStep = 17; //4
    const int largeStep = 18; //16
    const int xLStep = 64; //64
    const int xXLStep = 256; //256
    //-------------------------------Parameters-end-----------------------------------------

    double radiusM = radiusKm * 1000;
    int step = smallStep;
    // Turn lat and lon into indices in the pop data.
    const int leftX = ((leftLon + 180.0) / 360.0) * numCols;
    const int rightX = ((rightLon + 180.0) / 360.0) * numCols;
    const int upY = ((-upLat + 90.0) / 180.0) * numRows;
    const int downY = ((-downLat + 90.0) / 180.0) * numRows;

    for (int cenY = upY; cenY < downY; cenY += step) {
        if (cenY % printMod == 0) {
            cout << "Current lattitude: " << popTiff.coords(100, cenY, 0) << endl;
        }

        int kernelLength;
        int* kernel = makeKernel(1000, cenY, radiusM, popTiff, kernelLength); // Initializes kernelLength

        for (int cenX = leftX; cenX <= rightX; cenX += step) {
            double popWithinNKilometers = popWithinKernel(cenX, cenY, kernel, kernelLength, pop, popTiff);
            if (smallestPopMode) {
                if (popWithinNKilometers < smallest) {
                    cout << "Population within " << radiusKm << " kilometers of (" << (popTiff.coords(cenX, cenY, 0)) << ", " << popTiff.coords(cenX, cenY, 1) << "): " \
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
                    cout << "Population within " << radiusKm << " kilometers of (" << popTiff.coords(cenX, cenY, 0) << ", " << popTiff.coords(cenX, cenY, 1) << "): " \
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