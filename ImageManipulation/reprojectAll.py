import os
for i in range(899, 901):
    cmd = "gdalwarp -s_srs EPSG:4326 -t_srs ESRI:54012 -r near -of GTiff /mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps2020/" + str(i // 10) + "." + str(i % 10) + "PercentCircle.tif /mnt/c/Users/Administrator/Desktop/linuxPythonMadePercentMaps2020/reprojected/" + str(i) + ".tif"
    os.system(cmd)