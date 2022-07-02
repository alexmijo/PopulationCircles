import rasterio
import geopy.distance

circleToDraw = 1

circles = {1: (117, (23.629167, 90.2125)), 2: (177, (23.645833, 89.420833)),
           3: (239, (24.170833, 88.9875)), 4: (294, (24.529167, 87.795833)),
           5: (356, (24.695833, 87.9875)), 6: (425, (24.1875, 86.8125)),
           7: (490, (24.479167, 86.5125)), 8: (564, (24.479167, 85.804167)),
           9: (644, (24.379167, 85.020833)), 10: (705, (25.2125, 82.8125)),
           11: (757, (25.0875, 83.329167)), 12: (804, (23.4625, 82.620833)),
           13: (849, (23.670833, 82.229167)), 14: (901, (23.2625, 80.170833)),
           15: (950, (22.379167, 81.295833)), 16: (998, (25.179167, 79.795833)),
           17: (1046, (25.170833, 80.395833)), 18: (1105, (25.0125, 80.5125)),
           19: (1169, (23.2625, 79.904167)), 20: (1236, (23.9375, 79.129167)),
           21: (1308, (22.429167, 79.5125)), 22: (1415, (23.104167, 78.8125)),
           23: (1572, (23.020833, 81.129167)), 24: (1745, (22.6125, 80.7875)),
           25: (1838, (28.504167, 103.0125)), 26: (1884, (27.370833, 102.89583)),
           27: (1938, (27.3625, 102.2625)), 28: (1996, (27.620833, 101.72917))}

printDimensions = True
printProgress = True

# World land raster
landTif = rasterio.open('5000WideLand.tif')
land = landTif.read(1)
landH, landW = land.shape
# World borders raster
bordersTif = rasterio.open('5000WideBorders.tif')


if printDimensions:
    print('Land Map Height, Width in pixels: ' + str(landH) + ', ' + str(landW))

#cm = {
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
    'white': (255, 255, 255),
    0: (255, 0, 0),
    1: (255, 180, 180),
    2: (0, 230, 0),
    3: (186, 255, 186),
    4: (0, 0, 255),
    5: (180, 180, 255),
    'black': (0, 0, 0)}
darkestColors = {
    0: (97, 0, 0),
    2: (0, 94, 0),
    4: (0, 0, 97)}

shade1 = landTif.read(1)
shade2 = landTif.read(1)
shade3 = landTif.read(1)
borders = bordersTif.read(1)

for latPix in range(0, landH):
    if latPix % 10 == 0 and printProgress:
        print(latPix)
    for lonPix in range(0, landW):
        if borders[latPix][lonPix] != 0:
            shade1[latPix][lonPix] = cm['black'][0]
            shade2[latPix][lonPix] = cm['black'][1]
            shade3[latPix][lonPix] = cm['black'][2]
            continue
        
        center = circles[circleToDraw][1]
        distance = geopy.distance.distance((landTif.xy(latPix, lonPix)[1], landTif.xy(latPix, lonPix)[0]),
                                   center).km
        color = 0
        radius = circles[circleToDraw][0]

        if land[latPix][lonPix] != 0:
            shade1[latPix][lonPix] = 128
            shade2[latPix][lonPix] = 128
            shade3[latPix][lonPix] = 128
            if distance <= radius:
                shade1[latPix][lonPix] = cm[color][0]
                shade2[latPix][lonPix] = cm[color][1]
                shade3[latPix][lonPix] = cm[color][2]
        else:
            shade1[latPix][lonPix] = cm['white'][0]
            shade2[latPix][lonPix] = cm['white'][1]
            shade3[latPix][lonPix] = cm['white'][2]
            if distance > radius - 70 and distance <= radius:
                shade1[latPix][lonPix] = cm[color + 1][0]
                shade2[latPix][lonPix] = cm[color + 1][1]
                shade3[latPix][lonPix] = cm[color + 1][2]
        if distance <= 35:
            shade1[latPix][lonPix] = darkestColors[color][0]
            shade2[latPix][lonPix] = darkestColors[color][1]
            shade3[latPix][lonPix] = darkestColors[color][2]
        
colors = rasterio.open('/mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps/'
                       + str(circleToDraw) + "PercentCircle.tif", 'w', driver=landTif.driver,
                       height=landTif.height, width=landTif.width, count=3, dtype='uint8',
                       crs=landTif.crs, transform=landTif.transform)
colors.write(shade1, indexes=1)
colors.write(shade2, indexes=2)
colors.write(shade3, indexes=3)
colors.close()
landTif.close()
bordersTif.close()
