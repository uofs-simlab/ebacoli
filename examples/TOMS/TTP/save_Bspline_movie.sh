#!/usr/bin/env bash
#
# Save the most recently run movie, data, figs, and parameters
#

# File to parse for parameters
driverFile=driver-TT_monodomain-frameoutput.f95

# list of parameters in the source file (strings formatted for how they are defined in driverFile)
declare -a vars=("nin ="
		 "tstop ="
		 "ntout ="
		 "atol ="
		 "rtol =")

prompt for directory to store things
pr=1
while [ $pr -eq 1 ]
do
    echo "Save in which directory?"
    read dirName

    if [ -d "$dirName" ]; then
	echo "Directory exists! enter another name."
    else
	pr=0
    fi
done

mkdir $dirName
mkdir $dirName/data
mkdir $dirName/figs
mkdir $dirName/heatmaps
echo "# parameters in $driverFile:" > $dirName/parameters
for v in "${vars[@]}";
do
    grep "$v" $driverFile | grep -i -v format >> $dirName/parameters
done
mv Bsplines0*.png $dirName/figs
mv Bsplines?????? $dirName/data
mv heatmap_*.png $dirName/heatmaps
mv Bsplines_movie.mp4 $dirName/
