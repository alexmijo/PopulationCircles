import rasterio

printDimensions = True
printProgress = True

# World land raster
landTif = rasterio.open('5000WideLand.tif')
land = landTif.read(1)
landH, landW = land.shape
# World borders raster
bordersTif = rasterio.open('5000WideBorders.tif')
borders = bordersTif.read(1)

if printDimensions:
    print('Land Map Height, Width in pixels: ' + str(landH) + ', ' + str(landW))

# 0: Black (border)
# 1: Grey (land)
# 2: White (ocean)
# NA 3: Dark red (<= 35km of center)
# NA 4: Red (land in circle)
# NA 5: Light red ()

landNBorders = open('landNBorders.txt', 'w')

for latPix in range(0, landH):
    if latPix % 10 == 0 and printProgress:
        print(latPix)
    for lonPix in range(0, landW):
        if borders[latPix][lonPix] != 0:
            # Black
            landNBorders.write('0')
        elif land[latPix][lonPix] != 0:
            # Grey
            landNBorders.write('1')
        else:
            # White
            landNBorders.write('2')
    landNBorders.write('\n')

landNBorders.close()
landTif.close()
bordersTif.close()
