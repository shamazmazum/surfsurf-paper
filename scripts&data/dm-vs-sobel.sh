#!/bin/sh

convert images/Fss_mean_dir.png \
        images/sobel.png -gravity Center -geometry 800x800+600-600 -composite \
        images/distance_map.png -gravity Center -geometry 800x800+600+300 -composite \
        images/dm_sobel.png
