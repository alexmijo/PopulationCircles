// gdalstuff.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Needs GDAL library
// This is an early version of the program, the one used to make my first population circles map. It's still kinda
// spaghetti code and has code I copied from random places on the internet unattributed (GeoTiff and most of distance())
// Don't even have version control set up on this laptop yet.

#include <iostream>
#include <gdal.h>
#include "iostream"
#include "string"
#include "gdal_priv.h"
#include "cpl_conv.h"
//#include "gdalwarper.h"
#include "stdlib.h"
using namespace std;

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

    double* coords(int x, int y) {
        double* geoTransform = GetGeoTransform();
        double lon = geoTransform[0] + x * geoTransform[1] + y * geoTransform[2];
        double lat = geoTransform[3] + x * geoTransform[4] + y * geoTransform[5];
        double coordinates[2];
        coordinates[0] = lat;
        coordinates[1] = lon;
        return coordinates;
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

    float** GetRasterBand(int z) {

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

        float** bandLayer = new float* [NROWS];
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
    float** GetArray2D(int layerIndex, float** bandLayer) {

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

            bandLayer[row] = new float[NCOLS];
            for (int col = 0; col < NCOLS; col++) { // iterate through columns
                bandLayer[row][col] = (float)rowBuff[col];
            }
        }
        CPLFree(rowBuff);
        return bandLayer;
    }

};

double distance(double latp, double latc, double longp, double longc) {
    constexpr double req = 6378137.0;             //Radius at equator
    constexpr double flat = 1 / 298.257223563;    //flattening of earth
    constexpr double rpol = (1 - flat) * req;

    // I had to add this part
    if (latp == latc && longp == longc) {
        return 0.0;
    } else if (latp == 0.0 && latc == 0.0) {
        double lonDiff = abs(longc - longp);
        double fractionOfEquator = lonDiff / 360.0;
        double equatorLength = 2.0 * 3.14159265358979323 * req;
        return fractionOfEquator * equatorLength;
    }

    double sin_sigma, cos_sigma, sigma, sin_alpha, cos_sq_alpha, cos2sigma;
    double C, lam_pre;

    // convert to radians
    latp = M_PI * latp / 180.0;
    latc = M_PI * latc / 180.0;
    longp = M_PI * longp / 180.0;
    longc = M_PI * longc / 180.0;

    const double u1 = atan((1 - flat) * tan(latc));
    const double u2 = atan((1 - flat) * tan(latp));

    double lon = longp - longc;
    double lam = lon;
    double tol = pow(10., -12.); // iteration tolerance
    double diff = 1.;

    while (abs(diff) > tol) {
        sin_sigma = sqrt(pow((cos(u2) * sin(lam)), 2.) + pow(cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam), 2.));
        cos_sigma = sin(u1) * sin(u2) + cos(u1) * cos(u2) * cos(lam);
        sigma = atan(sin_sigma / cos_sigma);
        sin_alpha = (cos(u1) * cos(u2) * sin(lam)) / sin_sigma;
        cos_sq_alpha = 1 - pow(sin_alpha, 2.);
        cos2sigma = cos_sigma - ((2 * sin(u1) * sin(u2)) / cos_sq_alpha);
        C = (flat / 16) * cos_sq_alpha * (4 + flat * (4 - 3 * cos_sq_alpha));
        lam_pre = lam;
        lam = lon + (1 - C) * flat * sin_alpha * (sigma + C * sin_sigma * (cos2sigma + C * cos_sigma * (2 * pow(cos2sigma, 2.) - 1)));
        diff = abs(lam_pre - lam);
    }

    const double usq = cos_sq_alpha * ((pow(req, 2.) - pow(rpol, 2.)) / pow(rpol, 2.));
    const double A = 1 + (usq / 16384) * (4096 + usq * (-768 + usq * (320 - 175 * usq)));
    const double B = (usq / 1024) * (256 + usq * (-128 + usq * (74 - 47 * usq)));
    const double delta_sig = B * sin_sigma * (cos2sigma + 0.25 * B * (cos_sigma * (-1 + 2 * pow(cos2sigma, 2.)) -
        (1 / 6) * B * cos2sigma * (-3 + 4 * pow(sin_sigma, 2.)) *
        (-3 + 4 * pow(cos2sigma, 2.))));
    const double dis = rpol * A * (sigma - delta_sig);
    const double azi1 = atan2((cos(u2) * sin(lam)), (cos(u1) * sin(u2) - sin(u1) * cos(u2) * cos(lam)));

    return dis;
}

// TODO: make direction an enum
int boundingBoxEdge(const double cenLat, const double cenLon, const double radiusM, const int direction, Geotiff& popTiff) {
    const int numRows = popTiff.GetDimensions()[0];
    const int numCols = popTiff.GetDimensions()[1];
    // Turn lat and lon into indices in the pop data.
    const int x = ((cenLon + 180.0) / 360.0) * numCols;
    const int y = ((-cenLat + 90.0) / 180.0) * numRows;

    int edge;
    int edgeOfMap;
    // 0 = right, 1 = up, 2 = left, 3 = down
    // Start at center
    if (direction == 1 || direction == 3) {
        edge = y;
        edgeOfMap = numRows - 1;
    } else {
        edge = x;
        edgeOfMap = numCols - 1;
    }
    int incOrDec; // Added to edge to move one pixel in desired direction
    if (direction == 0 || direction == 3) {
        incOrDec = 1; // Right or down: increasing
    } else {
        incOrDec = -1; // Left or up: decreasing
    }

    while (edge >= 0 && edge <= edgeOfMap) {
        double lat, lon;
        if (direction == 1 || direction == 3) {
            lat = popTiff.coords(x, edge)[0];
            lon = popTiff.coords(x, edge)[1];
        } else {
            lat = popTiff.coords(edge, y)[0];
            lon = popTiff.coords(edge, y)[1];
        }
        if (distance(lat, cenLat, lon, cenLon) > radiusM) {
            edge -= incOrDec; // Went too far, walk it back
            break;
        }
        else {
            edge += incOrDec;
        }
    }
    return edge;
}

// Returns the population of a circle with a radius of <radiusKm> kilometers centered at the given latittude <cenLat>
// and longitude <cenLon>.
// Uses pop as the population data. popTiff must be the Geotiff that pop is from. Need both to avoid constantly making
// new arrays I think. Could alternatively modify the Geotiff class to just include the pointer to the array I think.
// pop must have first dimension corresponding to longitude, second corresponding to reverse lattitude
// TODO: Make a coordinates class
// TODO: figure out where to put const around pop type
float popWithinNKilometers(const double lat, const double lon, const double radiusKm, float** pop, Geotiff& popTiff) {
    const double radiusM = radiusKm * 1000; // Radius in meters
    const int numRows = popTiff.GetDimensions()[0]; // TODO: remove after parity
    const int numCols = popTiff.GetDimensions()[1];
    // Turn lat and lon into indices in the pop data.
    const int x = ((lon + 180.0) / 360.0) * numCols;
    const int y = ((-lat + 90.0) / 180.0) * numRows;

    const double cenLon = popTiff.coords(x, y)[1];
    const double cenLat = popTiff.coords(x, y)[0];

    // TODO: Major bug, can't find bounding box by only travelling left or right due to distortion
    //       up and down still work
    float thisPop = 0;

    int rightEdge = boundingBoxEdge(cenLat, cenLon, radiusM, 0, popTiff);
    int upEdge    = boundingBoxEdge(cenLat, cenLon, radiusM, 1, popTiff);
    int leftEdge  = boundingBoxEdge(cenLat, cenLon, radiusM, 2, popTiff);
    int downEdge  = boundingBoxEdge(cenLat, cenLon, radiusM, 3, popTiff);

    for (int X = leftEdge; X <= rightEdge; X++) {
        for (int Y = upEdge; Y <= downEdge; Y++) {
            double lattitude = popTiff.coords(X, Y)[0];
            double longitude = popTiff.coords(X, Y)[1];
            if (distance(lattitude, cenLat, longitude, cenLon) > radiusM) {
                continue;
            }
            // Indexing into pop requires Y first then X, unlike in Geotiff.coords() above
            float thisCellPop = pop[Y][X];
            if (thisCellPop >= 0.) {
                thisPop += thisCellPop;
            }
        }
    }
    return thisPop;
}

int main()
{
    // Load population data
    const string popDataFilename = "C:\\Users\\Administrator\\source\\repos\\gdalstuff\\gdalstuff\\GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0\\GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif";
    Geotiff popTiff(popDataFilename.c_str());

    // The population data as a 2d array. First dimension is longitude, second dimension is reverse lattitude.
    float** pop = popTiff.GetRasterBand(1);

    while (true) {
        double lat, lon, radiusKm;
        cout << "Enter radius in kilometers: ";
        cin >> radiusKm;
        cout << "Enter lattitude: ";
        cin >> lat;
        cout << "Enter longitude: ";
        cin >> lon;

        cout << "Population within " << radiusKm << " kilometers of (" << lat << ", " << lon << "): " \
             << popWithinNKilometers(lat, lon, radiusKm, pop, popTiff) << endl;
    }
}