rasterCircleFinder: rasterCircleFinder.cpp popSumTable2020.bin
	g++ -std=c++17 rasterCircleFinder.cpp -o rasterCircleFinder -I/usr/include/gdal -l gdal

popSumTable2020.bin: summationTableMaker NASA2020POPDATA.tif
	./summationTableMaker

summationTableMaker: summationTableMaker.cpp
	g++ -std=c++17 summationTableMaker.cpp -o summationTableMaker -I/usr/include/gdal -l gdal

clean:
	rm -f rasterCircleFinder
	rm -f popSumTable2020.bin
	rm -f summationTableMaker