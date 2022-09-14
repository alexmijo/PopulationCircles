import rasterio
import json

for circleNum in range(916, 917):
    circleToDraw = str(circleNum // 10) + "." + str(circleNum % 10)
    print(circleToDraw)

    printDimensions = True
    printProgress = True

    # World land raster
    landTif = rasterio.open('5000WideLand.tif')
    output1 = landTif.read(1)
    output2 = landTif.read(1)
    output3 = landTif.read(1)
    num_rows, num_cols = output1.shape

    if printDimensions:
        print('Land Map Height, Width in pixels: ' + str(num_rows) + ', ' + str(num_cols))

    # cm = {
    #    'white': (255, 255, 255),
    #    0: (255, 0, 0),
    #    1: (255, 180, 180),
    #    2: (0, 128, 255),
    #    3: (152, 204, 255),
    #    4: (255, 128, 0),
    #    5: (255, 201, 146),
    #    6: (0, 97, 0),
    #    7: (126, 195, 126),
    #    8: (200, 200, 0),
    #    9: (255, 255, 165),
    #    10: (127, 0, 255),
    #    11: (206, 155, 255),
    #    12: (0, 255, 0),
    #    13: (169, 255, 169),
    #    14: (0, 0, 255),
    #    15: (169, 169, 255),
    #    16: (0, 200, 200),
    #    17: (179, 255, 255),
    #    'black': (0, 0, 0)}

    cm = {
        2: (255, 255, 255),  # White
        4: (255, 0, 0),      # Red
        5: (255, 180, 180),  # Light red
        0: (0, 0, 0),        # Black
        3: (97, 0, 0),       # Dark red
        1: (128, 128, 128)   # Grey
    }

    colorsJSONFilename = "ColorsJSONFiles2020/colorsJSON" + circleToDraw + ".txt"
    colorsJSON = open(colorsJSONFilename, 'r')
    colors = json.load(colorsJSON)
    colorsJSON.close()

    for latPix in range(0, num_rows):
        if latPix % 10 == 0 and printProgress:
            print(latPix)
        for lonPix in range(0, num_cols):
            output1[latPix][lonPix] = cm[colors[latPix][lonPix]][0]
            output2[latPix][lonPix] = cm[colors[latPix][lonPix]][1]
            output3[latPix][lonPix] = cm[colors[latPix][lonPix]][2]

    map = rasterio.open('/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps2020/'
                        + circleToDraw + "PercentCircle.tif", 'w', driver=landTif.driver,
                        height=landTif.height, width=landTif.width, count=3, dtype='uint8',
                        crs=landTif.crs, transform=landTif.transform)
    map.write(output1, indexes=1)
    map.write(output2, indexes=2)
    map.write(output3, indexes=3)
    map.close()
    landTif.close()
