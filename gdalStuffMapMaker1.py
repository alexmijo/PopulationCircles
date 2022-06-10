import rasterio
import rasterio.plot
import numpy as np
import geopy.distance

rasterio.plot.get_plt() # import matplotlib.pyplot raise import error if matplotlib is not installed

centers = [(58.8333, 57.8333), (25.25, 64), (23.3333, 58.75), (27.25, 60), (28.9167, 64.4167), (22.25, 61.6667), (29.25, 64.9167), (36.3333, 65.0833), (40.75, 70.75),
           (41.2917, 77.5417), (25.083306, 103.417), (22.125, 101.167), (28.333306, 96.916694), (27.375, 95.5), (27.666694, 101.75), (22.75, 78.7083),
           (21.916694, 81.583306), (24.3583, 86.525)]
radii = [9000, 8500, 8000, 7500, 7000, 6500, 6000, 5500, 5000, 4500, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 500] # In kilometers

printDimensions = True
printProgress = True

# World land raster
landTif = rasterio.open('C:\Users\Administrator\Desktop\\thinkitworked.tif')
land = landTif.read(1)
[landH, landW] = land.shape
# World borders raster
bordersTif = rasterio.open('C:\Users\Administrator\Desktop\\testborders.tif')


if printDimensions:
    print('Land Map Height, Width in pixels: ' + str(landH) + ', ' + str(landW))
#    print('Population Map Height, Width in pixels: ' + str(popH) + ', ' + str(popW) + '\n')

cm = {
    'white': (255, 255, 255),
    0: (255, 0, 0),
    1: (255, 180, 180),
    2: (0, 128, 255),
    3: (152, 204, 255),
    4: (255, 128, 0),
    5: (255, 201, 146),
    6: (0, 97, 0),
    7: (126, 195, 126),
    8: (255, 255, 0),
    9: (255, 255, 165),
    10: (127, 0, 255),
    11: (206, 155, 255),
    12: (0, 255, 0),
    13: (169, 255, 169),
    14: (0, 0, 255),
    15: (169, 169, 255),
    16: (0, 255, 255),
    17: (179, 255, 255),
    'black': (0, 0, 0)}

shade1 = landTif.read(1)
shade2 = landTif.read(1)
shade3 = landTif.read(1)
borders = bordersTif.read(1)
print(shade1.shape)
rasterio.plot.show(landTif)



for latPix in range(0, landH):
    if latPix % 10 == 0 and printProgress:
        print(latPix)
    for lonPix in range(0, landW):
        if borders[latPix][lonPix] != 0:
            shade1[latPix][lonPix] = cm['black'][0]
            shade2[latPix][lonPix] = cm['black'][1]
            shade3[latPix][lonPix] = cm['black'][2]
            continue
        if land[latPix][lonPix] == 0:
            shade1[latPix][lonPix] = cm['white'][0]
            shade2[latPix][lonPix] = cm['white'][1]
            shade3[latPix][lonPix] = cm['white'][2]
        else:
            shade1[latPix][lonPix] = 128
            shade2[latPix][lonPix] = 128
            shade3[latPix][lonPix] = 128
        for i in range(18):
            color = (i % 9) * 2
            center = centers[i]
            radius = radii[i]
            try:
                distance = geopy.distance.vincenty((landTif.xy(latPix, lonPix)[1], landTif.xy(latPix, lonPix)[0]), center).km
            except:
                distance = 1000000000
            if distance < radius - 70 or distance > radius:
                continue
            if distance <= radius:
                if land[latPix][lonPix] == 0:
                    shade1[latPix][lonPix] = cm[color + 1][0]
                    shade2[latPix][lonPix] = cm[color + 1][1]
                    shade3[latPix][lonPix] = cm[color + 1][2]
                else:
                    shade1[latPix][lonPix] = cm[color][0]
                    shade2[latPix][lonPix] = cm[color][1]
                    shade3[latPix][lonPix] = cm[color][2]

        
colors = rasterio.open('C:\Users\Administrator\Desktop\outputrasterAllCircles2.tif', 'w', driver=landTif.driver, height=landTif.height, width=landTif.width, count=3, dtype='uint8', crs=landTif.crs, transform=landTif.transform)
colors.write(shade1, indexes=1)
colors.write(shade2, indexes=2)
colors.write(shade3, indexes=3)
colors.close()

colors = rasterio.open('C:\Users\Administrator\Desktop\outputrasterAllCircles2.tif')
rasterio.plot.show(colors)
