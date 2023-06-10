#include <cpl_conv.h>
#include <fstream>
#include <gdal.h>
#include <gdal_priv.h>
#include <iostream>
#include <stdlib.h>
#include <string>

// TODO: Eliminate this class entirely, if it makes sense to do so
class Geotiff {

  private:                       // NOTE: "private" keyword is redundant here.
                                 // we place it here for emphasis. Because these
                                 // variables are declared outside of "public",
                                 // they are private.
    const char *filename;        // name of Geotiff
    GDALDataset *geotiffDataset; // Geotiff GDAL datset object.
    double geotransform[6];      // 6-element geotranform array.
    int dimensions[3];           // X,Y, and Z dimensions.
    int NROWS, NCOLS, NLEVELS;   // dimensions of data in Geotiff.

  public:
    // define constructor function to instantiate object
    // of this Geotiff class.
    Geotiff(const char *tiffname) {
        filename = tiffname;
        GDALAllRegister();

        // set pointer to Geotiff dataset as class member.
        geotiffDataset = (GDALDataset *)GDALOpen(filename, GA_ReadOnly);

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

    const char *GetFileName() {
        /*
         * function GetFileName()
         * This function returns the filename of the Geotiff.
         */
        return filename;
    }

    const char *GetProjection() {
        /* function const char* GetProjection():
         *  This function returns a character array (std::string)
         *  for the projection of the geotiff file. Note that
         *  the "->" notation is used. This is because the
         *  "geotiffDataset" class variable is a pointer
         *  to an object or structure, and not the object
         *  itself, so the "." dot notation is not used.
         */
        return geotiffDataset->GetProjectionRef();
    }

    double *GetGeoTransform() {
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
        double *geoTransform = GetGeoTransform();
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

    int *GetDimensions() {
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

    double **GetRasterBand(int z) {

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
         * the Geotiff, cast the data to double**, and return
         * it to this function. This function returns that
         * double** pointer.
         */

        double **bandLayer = new double *[NROWS];
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

    template <typename T> double **GetArray2D(int layerIndex, double **bandLayer) {

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
        GDALDataType bandType = GDALGetRasterDataType(geotiffDataset->GetRasterBand(layerIndex));

        // get number of bytes per pixel in Geotiff
        int nbytes = GDALGetDataTypeSizeBytes(bandType);

        // allocate pointer to memory block for one row (scanline)
        // in 2D Geotiff array.
        T *rowBuff = (T *)CPLMalloc(nbytes * NCOLS);

        for (int row = 0; row < NROWS; row++) { // iterate through rows

            // read the scanline into the dynamically allocated row-buffer
            CPLErr e = geotiffDataset->GetRasterBand(layerIndex)
                           ->RasterIO(GF_Read, 0, row, NCOLS, 1, rowBuff, NCOLS, 1, bandType, 0, 0);
            if (!(e == 0)) {
                std::cout << "Warning: Unable to read scanline in Geotiff!" << std::endl;
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

// Turn the passed in data table into a summation table (mutates it)
void turnIntoSummationTable(double **data, const int numRows, int const numCols) {
    for (int x = 0; x < numCols; x++) {
        for (int y = 0; y < numRows; y++) {
            if (data[y][x] < 0) {
                data[y][x] = 0.0;
            }
            if (x == 0 && y == 0) {
                continue;
            } else if (x == 0) {
                data[y][x] += data[y - 1][x];
            } else if (y == 0) {
                data[y][x] += data[y][x - 1];
            } else {
                data[y][x] += data[y][x - 1] + data[y - 1][x] - data[y - 1][x - 1];
            }
        }
    }
}

// TODO: Make an option for 2020 or 2015 data
int main() {
    const int POP_NUM_ROWS = 2 * 60 * 180;
    const int POP_NUM_COLS = 2 * 60 * 360;
    // const int GDP_NUM_ROWS = 2 * 60 * 180;
    // const int GDP_NUM_COLS = 2 * 60 * 360;

    // Load tiff data
    const std::string popTiffFilename = "GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif";
    // const std::string gdpTiffFilename = "gdpPPPdata.tif"; // This is too large to fit in github
    // by the way (11GB)
    Geotiff popTiff(popTiffFilename.c_str());
    // Geotiff gdpTiff(gdpTiffFilename.c_str());

    int numRows = popTiff.GetDimensions()[0];
    int numCols = popTiff.GetDimensions()[1];
    if (popTiff.GetDimensions()[0] != POP_NUM_ROWS || popTiff.GetDimensions()[1] != POP_NUM_COLS) {
        std::cout << "Pop tiff isn't expected dimensions" << std::endl;
        return 1;
    }
    double **popData = popTiff.GetRasterBand(1);
    turnIntoSummationTable(popData, POP_NUM_ROWS, POP_NUM_COLS); // Mutates popdata

    std::fstream popSumTableFile;
    popSumTableFile.open("popSumTable.bin", std::ios::out | std::ios::binary);
    popSumTableFile.write(reinterpret_cast<char *>(&numRows), sizeof(int));
    popSumTableFile.write(reinterpret_cast<char *>(&numCols), sizeof(int));
    for (int r = 0; r < POP_NUM_ROWS; r++) {
        for (int c = 0; c < POP_NUM_COLS; c++) {
            popSumTableFile.write(reinterpret_cast<char *>(&popData[r][c]), sizeof(double));
        }
    }
    popSumTableFile.close();
}