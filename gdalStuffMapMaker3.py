import rasterio
import rasterio.plot
import numpy as np
import geopy.distance

rasterio.plot.get_plt() # import matplotlib.pyplot raise import error if matplotlib is not installed

centers = [(27.6917, 63.7), (28.6833, 99.7083), (28.5083, 103.008)]
radii = [6165, 3281, 1838] # In kilometers


numCircles = len(centers)
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
        distances = []
        for i in range(numCircles):
            center = centers[i]
            try:
                distances.append(geopy.distance.vincenty((landTif.xy(latPix, lonPix)[1], landTif.xy(latPix, lonPix)[0]), center).km)
            except:
                distances.append(1000000000)
        if land[latPix][lonPix] != 0:
            shade1[latPix][lonPix] = 128
            shade2[latPix][lonPix] = 128
            shade3[latPix][lonPix] = 128
            for i in range(numCircles):
                color = (i % 9) * 2
                radius = radii[i]
                distance = distances[i]
                #if distance < radius - 60 or distance > radius:
                #    continue
                if distance <= radius:
                    if land[latPix][lonPix] == 0:
                        shade1[latPix][lonPix] = cm[color + 1][0]
                        shade2[latPix][lonPix] = cm[color + 1][1]
                        shade3[latPix][lonPix] = cm[color + 1][2]
                    else:
                        shade1[latPix][lonPix] = cm[color][0]
                        shade2[latPix][lonPix] = cm[color][1]
                        shade3[latPix][lonPix] = cm[color][2]
            for i in range(numCircles - 1, -1, -1):
                color = (i % 9) * 2
                radius = radii[i]
                distance = distances[i]
                if distance > radius - 70 and distance <= radius:
                    shade1[latPix][lonPix] = cm[color][0]
                    shade2[latPix][lonPix] = cm[color][1]
                    shade3[latPix][lonPix] = cm[color][2]
            for i in range(numCircles):
                color = (i % 9) * 2
                distance = distances[i]
                if distance <= 35:
                    shade1[latPix][lonPix] = darkestColors[color][0]
                    shade2[latPix][lonPix] = darkestColors[color][1]
                    shade3[latPix][lonPix] = darkestColors[color][2]
        else:
            shade1[latPix][lonPix] = cm['white'][0]
            shade2[latPix][lonPix] = cm['white'][1]
            shade3[latPix][lonPix] = cm['white'][2]
            for i in range(numCircles - 1, -1, -1):
                color = (i % 9) * 2
                radius = radii[i]
                distance = distances[i]
                if distance > radius - 70 and distance <= radius:
                    shade1[latPix][lonPix] = cm[color + 1][0]
                    shade2[latPix][lonPix] = cm[color + 1][1]
                    shade3[latPix][lonPix] = cm[color + 1][2]
            for i in range(numCircles):
                color = (i % 9) * 2
                distance = distances[i]
                if distance <= 35:
                    shade1[latPix][lonPix] = darkestColors[color][0]
                    shade2[latPix][lonPix] = darkestColors[color][1]
                    shade3[latPix][lonPix] = darkestColors[color][2]
            

        
colors = rasterio.open('C:\Users\Administrator\Desktop\outputrasterAllCircles5.tif', 'w', driver=landTif.driver, height=landTif.height, width=landTif.width, count=3, dtype='uint8', crs=landTif.crs, transform=landTif.transform)
colors.write(shade1, indexes=1)
colors.write(shade2, indexes=2)
colors.write(shade3, indexes=3)
colors.close()

colors = rasterio.open('C:\Users\Administrator\Desktop\outputrasterAllCircles5.tif')
rasterio.plot.show(colors)
