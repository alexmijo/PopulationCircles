rasterCircleFinder: rasterCircleFinder.cpp popSumTable.bin
	g++ -std=c++17 rasterCircleFinder.cpp -o rasterCircleFinder

popSumTable.bin: summationTableMaker GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif
	./summationTableMaker

summationTableMaker: summationTableMaker.cpp
	g++ -std=c++17 summationTableMaker.cpp -o summationTableMaker -I/usr/include/gdal -l gdal

clean:
	rm -f rasterCircleFinder
	rm -f popSumTable.bin
	rm -f summationTableMaker