import os
for i in range(49, 132):
    cmd = "gdalwarp -s_srs EPSG:4326 -t_srs ESRI:54012 -r near -of GTiff /mnt/d/PopCirclesStuff/" + str(i) + ".tif /mnt/d/PopCirclesStuff/reprojected/" + str(i) + ".tif"
    os.system(cmd)