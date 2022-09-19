import os
for i in range(128, 838):
    cmd = "gdalwarp -s_srs EPSG:4326 -t_srs ESRI:54012 -r near -of GTiff /mnt/d/PopCirclesStuff/linuxPythonMadePercentMaps2020/" + str(i // 10) + "." + str(i % 10) + "PercentCircle.tif /mnt/d/PopCirclesStuff/linuxPythonMadePercentMaps2020/reprojected/" + str(i) + ".tif"
    os.system(cmd)