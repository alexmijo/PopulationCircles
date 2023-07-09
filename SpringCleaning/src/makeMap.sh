#!/bin/bash

./MapMaker > thing.ppm
pnmtopng thing.ppm > /mnt/c/Users/Administrator/Desktop/SpringCleaningMaps/testmap.png
rm thing.ppm