#!/bin/bash

# This script fully processes the repo directory with all
# SPOT 5 HRS stereo pairs, producing 1 DEM for each pair.


# Special instructions to use sparse_disp.
export PATH="${PATH}:/PATH/TO/STEREO/PIPELINE/bin"


# Constant paths.
dir_repo=`realpath "/PATH/TO/TREE/OF/STEREO/PAIRS/"`
seed_dem_path=`realpath "/PATH/TO/SEED/DEM/GEOTIFF"`
dir_out_path=`realpath "../"`
dir_out_dem_path=$dir_out_path/dem_v1
mkdir -p $dir_out_dem_path


# Constants for processing.
# Our CRS is EPSG:32642.
proj4="+proj=utm +zone=42 +datum=WGS84 +units=m +no_defs +type=crs"
mapproj_in_res=5
pllpar="--processes 1 --corr-tile-size 2000 --threads-singleprocess 8 --corr-timeout 7200 --corr-memory-limit-mb 32768 -t spot5maprpc --stereo-algorithm asp_mgm --cost-mode 3 --subpixel-mode 9 --corr-kernel 7 7 --subpixel-kernel 15 15"
point2dem_filter_par=""
dem_res=30

minZ=0
maxZ=7600


# The 4 lines below are used to process all scenes in the order they are found.
dir_repo_list=`ls $dir_repo`
dir_repo_scenes_n=`ls -l $dir_repo | wc -l`
echo "Found $(($dir_repo_scenes_n-1)) stereo pairs..."


# Processing loop.
for pair_name in $dir_repo_list;
do
 
  echo "Processing $pair_name..."
  cd $dir_repo/$pair_name
  front_dir=*HRS-1*
  back_dir=*HRS-2*

  add_spot_rpc $front_dir/SCENE01/METADATA.DIM --min-height $minZ --max-height $maxZ -o $front_dir/SCENE01/METADATA.DIM
  add_spot_rpc $back_dir/SCENE01/METADATA.DIM --min-height $minZ --max-height $maxZ -o $back_dir/SCENE01/METADATA.DIM

  mapproject -t rpc --t_srs "$proj4" --mpp $mapproj_in_res $seed_dem_path $front_dir/SCENE01/IMAGERY.TIF $front_dir/SCENE01/METADATA.DIM front_map_proj.tif
  mapproject -t rpc --t_srs "$proj4" --mpp $mapproj_in_res $seed_dem_path $back_dir/SCENE01/IMAGERY.TIF $back_dir/SCENE01/METADATA.DIM back_map_proj.tif
	
  parallel_stereo $pllpar front_map_proj.tif back_map_proj.tif $front_dir/SCENE01/METADATA.DIM $back_dir/SCENE01/METADATA.DIM tmp/out $seed_dem_path

  point2dem -r earth $point2dem_filter_par --tr $dem_res tmp/out-PC.tif
	mv tmp/out-DEM.tif $dir_out_dem_path/${pair_name}_DEM_${dem_res}m.tif

  cd ../..
	
	rm -r $dir_repo/$pair_name
	
	echo "=== Finished DEM production"
	echo "=== Currently in dir:"
	pwd
	
	
	dem_geoid --geoid EGM2008 $dir_out_dem_path/${pair_name}_DEM_${dem_res}m.tif
	
	rm $dir_out_dem_path/${pair_name}_DEM_${dem_res}m.tif
	rm *log-dem_geoid*.txt
	
	mv ${pair_name}_DEM_${dem_res}m-adj.tif $dir_out_dem_path
	
done
