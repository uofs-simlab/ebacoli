#!/usr/bin/env bash

# Objective:
# - run the driver,
# - plot the output,
# - store the results in a subdirectory

# Step 1: Run
echo "Running driver."
{ time ./driver-richards-frameoutput 2>&1 | tee driver.out ; } 2> time.out

# Step 2: Plot the results
# - may have to adjust wildcards
echo "Plotting output."
./plot_all_bsplines.py Bsplines000???
ffmpeg -i 'Bsplines%06d.png' -r 24 -vcodec libx264 Bsplines_movie.mp4 2>&1 > /dev/null

./plot_heatmap_from_bsplines.py Bsplines000???

# Step 3: Store results
echo "Enter a new dir to store results (Do not enter an existing dir!):"
read dirName

mkdir $dirName
mv Bsplines* $dirName
mv heatmap* $dirName
mv driver.out $dirName
mv time.out $dirName
