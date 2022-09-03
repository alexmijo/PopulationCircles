# PopulationCircles
Note to any potential employers: This code kinda got spaghettified and I haven't yet brought it up to my normal standards for code quality (although it does still work well). If you'd like to see a better example of clean C++ code I've written, look at https://github.com/alexmijo/BlackjackSim

This program, among other things, can beat the previous best known results for finding the smallest possible circles on earth containing a given percentage of the
world’s population. Also wrote code in C++, Python and Java to render the maps.
This software is what I used to make these maps:
https://www.reddit.com/r/dataisbeautiful/comments/ventld/oc_the_largest_and_smallest_possible_circles/
https://www.reddit.com/r/MapPorn/comments/vdv3vw/the_smallest_possible_circle_containing_10_of_the/
https://www.reddit.com/r/MapPorn/comments/vc85v3/oc_the_smallest_possible_circles_containing_25_50/
https://www.reddit.com/r/dataisbeautiful/comments/vc77av/oc_the_smallest_possible_circles_containing_25_50/
https://www.reddit.com/r/dataisbeautiful/comments/vaszmp/oc_the_smallest_possible_circles_containing_25_50/
https://www.reddit.com/r/MapPorn/comments/vas3hu/the_worlds_most_populous_circles_of_radius_1000km/
https://www.reddit.com/r/dataisbeautiful/comments/v9ei3v/oc_the_worlds_most_populous_circles_of_radius/

A much earlier version of this software (even earlier than the first commit of this project on github) is what I used to make these maps (and it wouldn't take much
modification to add this functionality back in):
https://www.reddit.com/r/dataisbeautiful/comments/rwxajo/each_stripe_contains_02_of_the_world_population/
https://www.reddit.com/r/dataisbeautiful/comments/qotq6p/each_stripe_contains_1_of_the_worlds_population/
https://www.reddit.com/r/MapPorn/comments/ny0yxr/1_population_bands_diagonal_edition/
https://www.reddit.com/r/dataisbeautiful/comments/nxyayy/oc_each_vertical_band_contains_01_of_the_worlds/
https://www.reddit.com/r/dataisbeautiful/comments/nwi16u/oc_each_horizontal_band_contains_01_of_the_worlds/
https://www.reddit.com/r/dataisbeautiful/comments/nvnmkn/oc_each_horizontal_band_contains_1_of_the_earths/
https://www.reddit.com/r/dataisbeautiful/comments/nvnlot/oc_each_vertical_band_contains_1_of_the_earths/
https://www.reddit.com/r/MapPorn/comments/ib26u9/10_bands_of_equal_population_latitude_and/
https://www.reddit.com/r/MapPorn/comments/iaext2/world_split_into_10_latitudinal_bands_of_equal/

Here's some more info on the code and how it was able to be relatively fast, copied from a reddit comment I made in one of the above links:

See this Wikipedia page for prior work on the 50% population circle: https://en.wikipedia.org/wiki/Valeriepieris_circle

The original viral map (the first image on that Wikipedia page) was cool but it almost kinda annoyed me when I saw it on reddit several years ago since they didn't account for the distortion of the map projection (it looked like a circle on that image but wouldn't look like a circle on a globe) and I didn't know if that was the smallest that they could've made that circle. A Singaporean professor named Danny Quah apparently also had the same thoughts, and he found a circle (that would actually be a circle on a globe) of radius 3300km instead of ~4000km like in the original image; that's the second image on that wikipedia page. I achieved a better result than Quah for the 50% circle (3281km instead of 3300km) since I analyzed the population data at a <1km resolution instead of 100km resolution like he did (I'd guess we actually used the same population data since he also used 2015 data and there aren't many competing datasets for this sort of stuff). I was able to do this without the code taking 10,000 to 100,000,000 times longer (1002 to 1004 , depending on what exactly it means to be analyzing the data at a 100km resolution) by using this technique https://en.wikipedia.org/wiki/Summed-area_table and generating a single circular kernel for each latitude. Even with this considerable speedup, the population data was so high resolution (much higher resolution than this image) that I had to run the program overnight. I find it interesting that unlike both the original circle and Quah's circle, my 50% circle doesn't include any of the island of Java (it's better to be further north to get more of northeastern China, Korea and Japan it seems).
